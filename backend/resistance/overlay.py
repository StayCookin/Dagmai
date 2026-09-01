"""Strain resistance overlay.

Maps a strain's annotated resistance mechanisms onto a per-drug "effective
multiplier" in (0, 1]: 1.0 means the drug reaches its target at full nominal
potency, values near 0 mean the mechanism has essentially neutralized it
(beta-lactamase hydrolysis, target bypass), and intermediate values represent
partial attenuation (efflux, reduced target affinity).

When a strain carries multiple mechanisms that each affect the same drug, we
take the minimum (most severe) multiplier rather than multiplying them
together -- e.g. an ESBL plus an efflux pump acting on ceftriaxone is
dominated by the beta-lactamase essentially destroying the drug; naive
multiplication of two already-small numbers would understate the mechanism
that alone is already sufficient to explain resistance.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..engine.drug_panel import ResistanceMechanism, Strain, load_resistance_mechanisms


@dataclass(frozen=True)
class DrugResistanceEffect:
    drug_id: str
    multiplier: float
    contributing_mechanisms: tuple[str, ...]


def strain_resistance_profile(strain: Strain) -> dict[str, DrugResistanceEffect]:
    """Return, for every drug affected by any of the strain's mechanisms, the
    effective multiplier and which mechanism(s) drove it."""
    all_mechanisms = load_resistance_mechanisms()
    per_drug: dict[str, list[tuple[float, str]]] = {}

    for mech_id in strain.mechanisms:
        mechanism: ResistanceMechanism | None = all_mechanisms.get(mech_id)
        if mechanism is None:
            continue
        for drug_id, multiplier in mechanism.affected_drugs.items():
            per_drug.setdefault(drug_id, []).append((multiplier, mechanism.id))

    profile: dict[str, DrugResistanceEffect] = {}
    for drug_id, entries in per_drug.items():
        min_multiplier = min(m for m, _ in entries)
        contributing = tuple(mid for m, mid in entries if m == min_multiplier)
        # Also keep any other mechanism that meaningfully contributes (within
        # 2x of the dominant one) so the reasoning layer can mention it.
        other = tuple(mid for m, mid in entries if mid not in contributing and m <= min_multiplier * 2)
        profile[drug_id] = DrugResistanceEffect(
            drug_id=drug_id,
            multiplier=min_multiplier,
            contributing_mechanisms=tuple(dict.fromkeys(contributing + other)),
        )
    return profile


def effective_multiplier(drug_id: str, profile: dict[str, DrugResistanceEffect]) -> float:
    effect = profile.get(drug_id)
    return 1.0 if effect is None else effect.multiplier
