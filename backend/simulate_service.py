"""Shared logic between the ranking and single-pair-detail endpoints: builds a
CombinationResultOut for one drug pair in one strain, at whatever level of
detail (cheap Bliss-only vs. full Bliss+FIC checkerboard) the caller needs.
"""

from __future__ import annotations

from .engine.drug_panel import Drug, Strain, load_resistance_mechanisms
from .engine.synergy import bliss_score, dose_response_curve, fic_index
from .reasoning.engine import RuleBasedReasoningEngine
from .resistance.overlay import DrugResistanceEffect, effective_multiplier, strain_resistance_profile
from .schemas import (
    BlissOut,
    CombinationResultOut,
    DrugOut,
    ExplanationOut,
    FICOut,
    ResistanceEffectOut,
)

_reasoning_engine = RuleBasedReasoningEngine()


def _drug_out(drug: Drug) -> DrugOut:
    return DrugOut(
        id=drug.id,
        display_name=drug.display_name,
        drug_class=drug.drug_class,
        target_gene_name=drug.target_gene_name,
        target_enzyme=drug.target_enzyme,
        pathway=drug.pathway,
        simulation_mode=drug.simulation_mode,
        notes=drug.notes,
    )


def _resistance_out(effect: DrugResistanceEffect | None) -> ResistanceEffectOut | None:
    if effect is None:
        return None
    mechanisms = load_resistance_mechanisms()
    names = [mechanisms[mid].name for mid in effect.contributing_mechanisms if mid in mechanisms]
    return ResistanceEffectOut(
        multiplier=effect.multiplier,
        contributing_mechanisms=list(effect.contributing_mechanisms),
        mechanism_names=names,
    )


def build_combination_result(
    drug_a: Drug,
    drug_b: Drug,
    strain: Strain,
    resistance_profile: dict[str, DrugResistanceEffect] | None = None,
    include_fic: bool = False,
) -> CombinationResultOut:
    profile = resistance_profile if resistance_profile is not None else strain_resistance_profile(strain)
    resistance_a = profile.get(drug_a.id)
    resistance_b = profile.get(drug_b.id)
    mult_a = effective_multiplier(drug_a.id, profile)
    mult_b = effective_multiplier(drug_b.id, profile)

    bliss = None
    fic = None
    if drug_a.simulation_mode == "direct" and drug_b.simulation_mode == "direct":
        curve_a = dose_response_curve(drug_a.id)
        curve_b = dose_response_curve(drug_b.id)
        if curve_a.ic50_sim is not None and curve_b.ic50_sim is not None:
            # Dose each drug at maximal simulated inhibition strength,
            # attenuated only by how much this strain's resistance blunts it
            # -- "give the strongest reasonable dose; how much of it actually
            # lands". Deliberately *not* each drug's own IC50_sim here: doing
            # so would normalize away potency differences (every unaffected
            # direct-mode drug alone would land at exactly 50% growth by
            # construction), making the ranking dominated by ties. IC50_sim
            # is still used below for the ΣFIC checkerboard, where a common
            # per-drug reference dose is exactly what a checkerboard needs.
            strength_a = min(1.0, mult_a)
            strength_b = min(1.0, mult_b)
            bliss = bliss_score(drug_a, strength_a, drug_b, strength_b)
            if include_fic:
                fic = fic_index(drug_a, drug_b, resistance_multiplier_a=mult_a, resistance_multiplier_b=mult_b)

    explanation = _reasoning_engine.explain_combination(
        drug_a, drug_b, strain, bliss, fic, resistance_a, resistance_b
    )

    return CombinationResultOut(
        drug_a=_drug_out(drug_a),
        drug_b=_drug_out(drug_b),
        resistance_a=_resistance_out(resistance_a),
        resistance_b=_resistance_out(resistance_b),
        bliss=(
            BlissOut(
                growth_fraction_a_alone=bliss.growth_fraction_a_alone,
                growth_fraction_b_alone=bliss.growth_fraction_b_alone,
                growth_fraction_combo=bliss.growth_fraction_combo,
                expected_independent_inhibition=bliss.expected_independent_inhibition,
                observed_inhibition=bliss.observed_inhibition,
                bliss_score=bliss.bliss_score,
                classification=bliss.classification,
            )
            if bliss
            else None
        ),
        fic=(
            FICOut(
                fic_a=fic.fic_a,
                fic_b=fic.fic_b,
                sigma_fic=fic.sigma_fic,
                classification=fic.classification,
            )
            if fic
            else None
        ),
        explanation=ExplanationOut(
            headline=explanation.headline,
            rationale=explanation.rationale,
            confidence=explanation.confidence,
            caveats=list(explanation.caveats),
        ),
    )


def ranking_sort_key(result: CombinationResultOut) -> tuple[int, float, float]:
    """Lower sorts first. Primary group: pairs with a numeric FBA growth
    prediction (lower combined growth = better = group 0), ordered by that
    growth fraction ascending -- see synergy.py docstring for why ranking by
    Bliss score alone would be misleading here. Both drugs dosed at maximal
    simulated strength mean many pairs tie at growth_fraction_combo == 0.0
    (any pair containing one resistance-unaffected direct-mode drug already
    fully suppresses growth on its own); Bliss score breaks those ties
    (higher/more-synergistic first) so the ordering among ties isn't
    arbitrary, even though it structurally skews non-positive (see above).
    Pairs with no FBA signal at all (group 1) sort after, in a stable
    (insertion) order.
    """
    if result.bliss is not None:
        return (0, result.bliss.growth_fraction_combo, -result.bliss.bliss_score)
    return (1, 0.0, 0.0)
