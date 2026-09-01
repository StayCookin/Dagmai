"""Core FBA simulation: applies drug inhibition(s) as reaction capacity constraints
and reports predicted growth relative to the unperturbed model.

Modeling choice (see README "Layer 1 modeling choices" for the full derivation):
a drug with `simulation_mode="direct"` is represented as a capacity cap on its
target reaction(s), sized relative to the flux those reactions actually carry
near the wild-type optimum (via tight-bound FVA, cached in model_loader). A
naive cap relative to the model's generic +/-1000 default bound produces a
step function (no effect until the drug is at ~100% strength) because default
bounds are ~40,000x larger than the flux any single essential reaction here
actually needs -- capacity-relative scaling gives a graded, monotonic
dose-response instead, which is what makes Bliss/FIC scoring meaningful.

`strength` is a dimensionless 0..1 "effective inhibition" knob, not a real
drug concentration -- see README for why a true PK-anchored dose axis is out
of scope for this MVP.
"""

from __future__ import annotations

from dataclasses import dataclass

from .drug_panel import Drug, all_target_reactions
from .model_loader import CAPACITY_HEADROOM, FALLBACK_TIMEOUT_SECONDS, get_loaded_model


@dataclass(frozen=True)
class GrowthResult:
    growth: float
    growth_fraction: float  # relative to wild-type, clamped to [0, 1]


def _clamp_growth(value: float | None, wt_growth: float) -> float:
    if value is None or value != value:  # None or NaN (infeasible LP)
        return 0.0
    return max(0.0, value)


def _apply_drug(model, capacities: dict[str, float], drug: Drug, strength: float) -> None:
    """Mutate reaction bounds in-place for the duration of a `with model:` block.

    Tightens toward the drug's cap rather than overwriting outright, so that
    when two drugs in the same combination share a target reaction (e.g. two
    beta-lactams both constraining the PBP transpeptidation reactions), the
    more restrictive of the two caps wins instead of whichever drug happened
    to be applied last silently discarding the other's constraint.
    """
    if strength <= 0 or drug.simulation_mode != "direct":
        return
    strength = min(strength, 1.0)
    for rid in drug.reactions:
        reaction = model.reactions.get_by_id(rid)
        cur_lb, cur_ub = reaction.bounds
        cap = capacities[rid] * (1 - strength) * CAPACITY_HEADROOM
        new_lb = max(cur_lb, -cap) if cur_lb < 0 else 0.0
        new_ub = min(cur_ub, cap) if cur_ub > 0 else 0.0
        reaction.bounds = (new_lb, new_ub)


def _solve_with_fallback(model) -> float | None:
    """slim_optimize(), retrying once with presolve forced on if the fast
    (presolve-off, short-timeout) attempt hits its timeout instead of
    reaching optimality -- see model_loader.py for why. The vast majority of
    solves never take this branch."""
    raw = model.slim_optimize()
    if model.solver.status == "optimal":
        return raw
    original_timeout = model.solver.configuration.timeout
    model.solver.configuration.presolve = True
    model.solver.configuration.timeout = FALLBACK_TIMEOUT_SECONDS
    try:
        return model.slim_optimize()
    finally:
        model.solver.configuration.presolve = False
        model.solver.configuration.timeout = original_timeout


def simulate_growth(drug_strengths: list[tuple[Drug, float]]) -> GrowthResult:
    """Run FBA with one or more drugs applied simultaneously.

    `drug_strengths` is a list of (Drug, effective_strength) pairs. Drugs with
    simulation_mode != "direct" are accepted but contribute no constraint
    (their combined effect is out of FBA's scope; see the reasoning layer for
    how they're still surfaced to the user).
    """
    loaded = get_loaded_model(tuple(all_target_reactions()))
    model = loaded.model
    with model:
        for drug, strength in drug_strengths:
            _apply_drug(model, loaded.reaction_capacity, drug, strength)
        raw = _solve_with_fallback(model)
    growth = _clamp_growth(raw, loaded.wt_growth)
    fraction = 0.0 if loaded.wt_growth <= 0 else min(1.0, growth / loaded.wt_growth)
    return GrowthResult(growth=growth, growth_fraction=fraction)


def wild_type_growth() -> float:
    return get_loaded_model(tuple(all_target_reactions())).wt_growth


def has_capacity_data(reaction_id: str) -> bool:
    loaded = get_loaded_model(tuple(all_target_reactions()))
    return reaction_id in loaded.reaction_capacity
