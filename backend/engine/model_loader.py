"""Loads the genome-scale metabolic model and caches per-reaction flux capacities.

Model choice: iJO1366 (E. coli K-12 MG1515, 2583 reactions / 1367 genes), the
direct predecessor of iML1515. It ships bundled inside the `cobra` package, so
it loads with no network access. iML1515 itself is only distributed from the
BiGG Models repository (bigg.ucsd.edu), which is not reachable from this
sandboxed environment's network policy -- see README "Honest constraints" for
how to swap it in once that access exists. iJO1366 and iML1515 share the same
core reaction network and gene lineage, so target genes/reactions used here
carry over with at most cosmetic ID changes.
"""

from __future__ import annotations

import os
import threading
from dataclasses import dataclass
from functools import lru_cache

import cobra
from cobra.flux_analysis import flux_variability_analysis

# How close to the true optimum we require when measuring a target reaction's
# achievable flux range. Kept tight (99.9%) so the resulting "capacity"
# reflects what's actually needed for near-optimal growth rather than what's
# theoretically reachable if growth were allowed to degrade -- see the
# derivation notes in README.md ("Layer 1 modeling choices").
FVA_FRACTION_OF_OPTIMUM = 0.999

# Multiplicative slack applied on top of the measured capacity so that
# strength=0 (no drug) never nudges growth below the true wild-type optimum
# due to solver/FVA tolerance noise.
CAPACITY_HEADROOM = 1.15

_MODEL_PATH = os.path.join(os.path.dirname(cobra.__file__), "data", "iJO1366.xml.gz")

_lock = threading.Lock()


@dataclass(frozen=True)
class LoadedModel:
    model: cobra.Model
    wt_growth: float
    reaction_capacity: dict[str, float]


def _compute_capacities(model: cobra.Model, reaction_ids: list[str]) -> dict[str, float]:
    if not reaction_ids:
        return {}
    fva = flux_variability_analysis(
        model, reaction_list=reaction_ids, fraction_of_optimum=FVA_FRACTION_OF_OPTIMUM
    )
    capacities: dict[str, float] = {}
    for rid in reaction_ids:
        row = fva.loc[rid]
        capacities[rid] = max(abs(row["minimum"]), abs(row["maximum"]), 1e-6)
    return capacities


# GLPK's simplex occasionally warm-starts a new LP from a very different
# prior basis (after hundreds of sequential solves against the same shared
# model/solver, each with different tightened reaction bounds -- this shows
# up reliably once a reaction with a very small baseline flux, like
# colistin's lipid-A-pathway target, gets scanned across many strengths
# after ~150+ prior unrelated solves) and can then take a very long time --
# tens of seconds to minutes -- to pivot back to optimality on what is
# otherwise an easy problem. Enabling presolve avoids this by giving GLPK a
# fresh simplified basis instead of warm-starting, but costs ~2.5s on every
# single solve if left on permanently, which is too slow for the common
# (non-degenerate) case. So: a short timeout with presolve OFF is the fast
# path for ordinary solves, and `solve_with_fallback` in fba.py retries with
# presolve ON only on the rare solve that actually hits it. Reproduced and
# diagnosed empirically; see README "Layer 1 modeling choices" -> performance
# note.
FAST_TIMEOUT_SECONDS = 3
FALLBACK_TIMEOUT_SECONDS = 15


@lru_cache(maxsize=1)
def _load_base_model() -> cobra.Model:
    model = cobra.io.read_sbml_model(_MODEL_PATH)
    model.solver.configuration.timeout = FAST_TIMEOUT_SECONDS
    model.solver.configuration.presolve = False
    return model


@lru_cache(maxsize=4)
def _get_loaded_model_cached(target_reaction_ids: tuple[str, ...]) -> LoadedModel:
    model = _load_base_model()
    wt_growth = model.slim_optimize()
    capacities = _compute_capacities(model, sorted(set(target_reaction_ids)))
    return LoadedModel(model=model, wt_growth=wt_growth, reaction_capacity=capacities)


def get_loaded_model(target_reaction_ids: tuple[str, ...]) -> LoadedModel:
    """Return the cached base model plus capacities for the given reactions.

    `target_reaction_ids` should be the union of every reaction referenced by
    the drug panel; FVA-derived capacities are expensive (each reaction costs
    two extra LP solves) so they are computed once per distinct reaction set
    and cached -- `simulate_growth` calls this on every single FBA run, so
    without caching here every growth prediction would silently redo FVA from
    scratch.
    """
    with _lock:
        return _get_loaded_model_cached(tuple(sorted(set(target_reaction_ids))))
