"""Layer 3 (reasoning) for the MVP: a deterministic, rule-based explanation
generator grounded in the actual FBA + resistance-overlay outputs for a given
combination -- not a template-only mock. It reads the same Bliss/FIC numbers
and resistance mechanism metadata the rest of the app computed and turns them
into a mechanistic sentence.

Why not the originally-scoped Llama 3.3 70B / Qwen 2.5 via Ollama or vLLM:
that needs either a GPU-hosted model server or outbound access to a model
weights host, neither of which this sandboxed build environment has (see
README "Honest constraints"). `ReasoningEngine` is the seam meant to carry
that upgrade -- swap `RuleBasedReasoningEngine` for an `LLMReasoningEngine`
that prompts a real model with the same structured inputs, and the rest of
the app (routers, frontend) doesn't need to change.

ChemBERTa/MolFormer (SMILES) and ESM-2 (protein) embeddings from the original
Layer 3 spec are not wired in for the same reason (they need model weights
downloaded from HuggingFace, which this environment cannot reach). The
`literature_synergy` table below is the honest fallback for what those
embeddings would otherwise help rank -- a small hand-curated set of published
synergy findings for pairs this MVP's FBA layer cannot score numerically.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field

from ..engine.drug_panel import Drug, ResistanceMechanism, Strain, load_resistance_mechanisms
from ..engine.synergy import BlissResult, FICResult
from ..resistance.overlay import DrugResistanceEffect

# Hand-curated, literature-referenced synergy notes for pairs where at least
# one drug is simulation_mode="non_metabolic" and so has no Bliss/FIC score.
# Keys are frozenset({drug_id_a, drug_id_b}). This is explicitly a small,
# static stand-in for the DrugComb / literature-mined training data the
# original Layer 3 spec called for -- see README.
LITERATURE_SYNERGY_NOTES: dict[frozenset[str], str] = {
    frozenset({"ciprofloxacin", "fosfomycin"}): (
        "Fosfomycin-induced cell wall/precursor stress and fluoroquinolone-induced "
        "DNA damage act on unrelated essential processes; combination use in complicated "
        "UTI/MDR uropathogens is reported in clinical literature, though evidence quality "
        "varies by study."
    ),
    frozenset({"colistin", "meropenem"}): (
        "Colistin permeabilizes the outer membrane, which is reported to improve carbapenem "
        "penetration in some MDR Enterobacterales and Pseudomonas studies; a mechanistically "
        "plausible pairing even where FBA cannot score the carbapenem-membrane interaction directly."
    ),
    frozenset({"gentamicin", "ampicillin"}): (
        "Classic cell-wall-active + aminoglycoside pairing: beta-lactam-induced cell wall "
        "disruption is reported to increase aminoglycoside uptake, a mechanism outside what a "
        "metabolic-flux model can represent but well documented for enterococci and some "
        "Enterobacterales."
    ),
    frozenset({"amikacin", "ampicillin"}): (
        "Same rationale as gentamicin + ampicillin: beta-lactam cell wall disruption is reported "
        "to aid aminoglycoside entry."
    ),
}


@dataclass(frozen=True)
class CombinationExplanation:
    headline: str
    rationale: str
    confidence: str  # "fba_grounded" | "fba_partial" | "literature_only" | "insufficient_data"
    caveats: tuple[str, ...] = field(default_factory=tuple)


class ReasoningEngine(ABC):
    """Interface a real LLM-backed engine (Ollama/vLLM + Llama 3.3 70B or
    Qwen 2.5, per the original Layer 3 design) would implement in place of
    RuleBasedReasoningEngine, taking the same structured inputs."""

    @abstractmethod
    def explain_combination(
        self,
        drug_a: Drug,
        drug_b: Drug,
        strain: Strain,
        bliss: BlissResult | None,
        fic: FICResult | None,
        resistance_a: DrugResistanceEffect | None,
        resistance_b: DrugResistanceEffect | None,
    ) -> CombinationExplanation:
        raise NotImplementedError


def _mechanism_names(effect: DrugResistanceEffect | None) -> list[str]:
    if effect is None:
        return []
    mechanisms = load_resistance_mechanisms()
    return [mechanisms[mid].name for mid in effect.contributing_mechanisms if mid in mechanisms]


def _resistance_clause(drug: Drug, effect: DrugResistanceEffect | None) -> str | None:
    if effect is None:
        return None
    names = ", ".join(_mechanism_names(effect))
    pct_active = round(effect.multiplier * 100)
    if effect.multiplier <= 0.15:
        return f"{drug.display_name} is largely neutralized in this strain by {names} (~{pct_active}% of nominal potency remains)."
    if effect.multiplier <= 0.6:
        return f"{drug.display_name} is partially attenuated in this strain by {names} (~{pct_active}% of nominal potency remains)."
    return f"{drug.display_name} retains most of its activity despite {names} (~{pct_active}% of nominal potency remains)."


class RuleBasedReasoningEngine(ReasoningEngine):
    def explain_combination(
        self,
        drug_a: Drug,
        drug_b: Drug,
        strain: Strain,
        bliss: BlissResult | None,
        fic: FICResult | None,
        resistance_a: DrugResistanceEffect | None,
        resistance_b: DrugResistanceEffect | None,
    ) -> CombinationExplanation:
        pair_key = frozenset({drug_a.id, drug_b.id})
        caveats: list[str] = []

        target_clause = (
            f"{drug_a.display_name} inhibits {drug_a.target_enzyme} ({drug_a.pathway}); "
            f"{drug_b.display_name} inhibits {drug_b.target_enzyme} ({drug_b.pathway})."
        )
        same_pathway = drug_a.pathway == drug_b.pathway and drug_a.simulation_mode == "direct" == drug_b.simulation_mode

        res_clauses = [
            c for c in (_resistance_clause(drug_a, resistance_a), _resistance_clause(drug_b, resistance_b)) if c
        ]
        resistance_text = " ".join(res_clauses)

        if bliss is not None and fic is not None:
            headline = f"{drug_a.display_name} + {drug_b.display_name}: {bliss.classification} (Bliss), ΣFIC={fic.sigma_fic:.2f} ({fic.classification})"
            rationale_parts = [
                target_clause,
                f"FBA-predicted growth remaining: {drug_a.display_name} alone "
                f"{bliss.growth_fraction_a_alone:.0%}, {drug_b.display_name} alone {bliss.growth_fraction_b_alone:.0%}, "
                f"combined {bliss.growth_fraction_combo:.0%} (Bliss-independent expectation was "
                f"{1 - bliss.expected_independent_inhibition:.0%} growth remaining).",
            ]
            if resistance_text:
                rationale_parts.append(resistance_text)
            if same_pathway:
                caveats.append(
                    "Both targets sit on the same linear metabolic pathway in this model; steady-state FBA "
                    "captures shared-bottleneck capacity but not kinetic/sequential-depletion effects, so it "
                    "can under-predict synergy well documented in vitro for pathway pairs like this one."
                )
            confidence = "fba_grounded"
            rationale = " ".join(rationale_parts)
        elif bliss is not None:
            headline = f"{drug_a.display_name} + {drug_b.display_name}: {bliss.classification} (Bliss); ΣFIC not available"
            rationale = " ".join([target_clause, resistance_text]).strip()
            confidence = "fba_partial"
            caveats.append("Checkerboard ΣFIC could not be computed (no IC50_sim crossing for one drug in this strain).")
        else:
            non_metabolic = [d.display_name for d in (drug_a, drug_b) if d.simulation_mode != "direct"]
            note = LITERATURE_SYNERGY_NOTES.get(pair_key)
            headline = f"{drug_a.display_name} + {drug_b.display_name}: not FBA-simulable"
            rationale_parts = [target_clause]
            if note:
                rationale_parts.append(note)
            if resistance_text:
                rationale_parts.append(resistance_text)
            rationale = " ".join(rationale_parts)
            confidence = "literature_only" if note else "insufficient_data"
            caveats.append(
                f"{' and '.join(non_metabolic)} target processes (ribosome, DNA gyrase, or a multi-target "
                "prodrug activation) not represented as reactions in this genome-scale reconstruction, so no "
                "growth-based score is available -- see README 'Honest constraints'."
            )

        return CombinationExplanation(
            headline=headline,
            rationale=rationale,
            confidence=confidence,
            caveats=tuple(caveats),
        )
