"""Dose-response scanning, Bliss independence scoring, and a checkerboard-style
FIC index for pairs of drugs in the panel.

All three pieces of math operate on `simulate_growth` from fba.py, so they
inherit its scope limits: only pairs where both drugs are
`simulation_mode="direct"` produce numeric scores here. See README for how
non-metabolic drugs (fluoroquinolones, aminoglycosides, nitrofurantoin) are
still surfaced through the reasoning layer instead.

Empirical finding worth knowing before reading scores off this module: for
two *independently* essential single-bottleneck reactions (which is most of
this panel), steady-state FBA's combined-growth optimum consistently lands
at approximately min(growth_alone_A, growth_alone_B) -- the LP is bounded by
whichever constraint is tighter, and unless the two reactions are actually
metabolically coupled (shared cofactor pools, forced flux rerouting), adding
a second cap rarely reduces growth further than the first cap already did.
That makes Bliss scores skew additive-to-antagonistic by construction here,
even for pairs with strong published wet-lab synergy (trimethoprim +
sulfamethoxazole scores as antagonistic in this model, despite being the
textbook synergy example). This isn't a bug -- it's a real, useful
finding about what static capacity-constraint FBA can and can't show, and is
exactly the "Bliss/FIC as prioritisation signal, not efficacy prediction"
caveat from the README made concrete. See README "Honest constraints" for
the full discussion and why the ranking endpoint therefore sorts primarily
by predicted combined growth suppression rather than by Bliss score.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

from .drug_panel import Drug
from .fba import GrowthResult, simulate_growth

# Number of points used when scanning strength 0..1 to locate IC50_sim
# (the strength at which growth crosses 50% of wild-type). 21 points gives
# ~0.05 resolution, which is enough for a monotonic capacity-constraint curve.
DOSE_SCAN_POINTS = 21

# Target growth-fraction threshold used to define both IC50_sim and the FIC
# checkerboard's inhibitory endpoint. 0.5 mirrors the "50% growth inhibition"
# convention used for in vitro IC50s; a true MIC (no visible growth) would sit
# further down this curve, but 0.5 gives a numerically stable crossing point
# for as many of the panel's drugs as possible.
INHIBITORY_THRESHOLD = 0.5

# Checkerboard grid resolution (fraction of each drug's own IC50_sim tested).
# Extended slightly past 1.0 so a drug that needs *more* than its reference
# potency to reach the threshold in a resistant strain still has a chance to
# be found, rather than silently reporting "not achievable".
_FIC_GRID_FRACTIONS = [round(0.2 * i, 2) for i in range(0, 11)]  # 0.0 .. 2.0


@dataclass(frozen=True)
class DoseResponsePoint:
    strength: float
    growth_fraction: float


@dataclass(frozen=True)
class DoseResponseCurve:
    drug_id: str
    points: tuple[DoseResponsePoint, ...]
    ic50_sim: float | None  # strength at which growth_fraction crosses 0.5, or None


def _growth_fraction_alone(drug: Drug, strength: float) -> float:
    return simulate_growth([(drug, strength)]).growth_fraction


@lru_cache(maxsize=None)
def dose_response_curve(drug_id: str) -> DoseResponseCurve:
    from .drug_panel import load_drugs

    drug = load_drugs()[drug_id]
    if drug.simulation_mode != "direct":
        return DoseResponseCurve(drug_id=drug_id, points=(), ic50_sim=None)

    points = []
    for i in range(DOSE_SCAN_POINTS):
        strength = i / (DOSE_SCAN_POINTS - 1)
        points.append(DoseResponsePoint(strength=strength, growth_fraction=_growth_fraction_alone(drug, strength)))

    ic50 = _interpolate_crossing(points, INHIBITORY_THRESHOLD)
    return DoseResponseCurve(drug_id=drug_id, points=tuple(points), ic50_sim=ic50)


def _interpolate_crossing(points: list[DoseResponsePoint], threshold: float) -> float | None:
    for prev, curr in zip(points, points[1:]):
        if prev.growth_fraction >= threshold >= curr.growth_fraction and prev.growth_fraction != curr.growth_fraction:
            span = prev.growth_fraction - curr.growth_fraction
            frac = (prev.growth_fraction - threshold) / span
            return prev.strength + frac * (curr.strength - prev.strength)
    return None


@dataclass(frozen=True)
class BlissResult:
    growth_fraction_a_alone: float
    growth_fraction_b_alone: float
    growth_fraction_combo: float
    expected_independent_inhibition: float
    observed_inhibition: float
    bliss_score: float  # observed - expected; >0 synergy, ~0 additive, <0 antagonism

    @property
    def classification(self) -> str:
        if self.bliss_score > 0.1:
            return "synergy"
        if self.bliss_score < -0.1:
            return "antagonism"
        return "additive"


def bliss_score(drug_a: Drug, strength_a: float, drug_b: Drug, strength_b: float) -> BlissResult | None:
    if drug_a.simulation_mode != "direct" or drug_b.simulation_mode != "direct":
        return None
    fa = _growth_fraction_alone(drug_a, strength_a)
    fb = _growth_fraction_alone(drug_b, strength_b)
    fab = simulate_growth([(drug_a, strength_a), (drug_b, strength_b)]).growth_fraction

    ia, ib, iab = 1 - fa, 1 - fb, 1 - fab
    expected = ia + ib - ia * ib
    score = iab - expected
    return BlissResult(
        growth_fraction_a_alone=fa,
        growth_fraction_b_alone=fb,
        growth_fraction_combo=fab,
        expected_independent_inhibition=expected,
        observed_inhibition=iab,
        bliss_score=score,
    )


@dataclass(frozen=True)
class FICResult:
    fic_a: float
    fic_b: float
    sigma_fic: float
    classification: str  # "synergy" | "additive" | "indifferent" | "antagonism"
    frac_a_used: float
    frac_b_used: float


def _classify_fic(sigma: float) -> str:
    if sigma <= 0.5:
        return "synergy"
    if sigma <= 1.0:
        return "additive"
    if sigma <= 4.0:
        return "indifferent"
    return "antagonism"


def fic_index(
    drug_a: Drug,
    drug_b: Drug,
    resistance_multiplier_a: float = 1.0,
    resistance_multiplier_b: float = 1.0,
) -> FICResult | None:
    """Checkerboard-style ΣFIC search.

    Each drug's own IC50_sim (susceptible-reference potency) is scaled by the
    strain's resistance multiplier for that drug, then a grid of dose
    fractions (0..2x that strain-adjusted reference) is searched for the
    combination that reaches the inhibitory threshold using the smallest
    total fraction. Returns None if either drug is non-metabolic or has no
    IC50_sim, or if no grid point reaches the threshold (i.e. the strain is
    resistant enough that neither drug, nor their combination within the
    tested range, is predicted to work).
    """
    curve_a = dose_response_curve(drug_a.id)
    curve_b = dose_response_curve(drug_b.id)
    if curve_a.ic50_sim is None or curve_b.ic50_sim is None:
        return None

    best: FICResult | None = None
    for frac_a in _FIC_GRID_FRACTIONS:
        nominal_a = frac_a * curve_a.ic50_sim
        effective_a = min(1.0, nominal_a * resistance_multiplier_a)
        for frac_b in _FIC_GRID_FRACTIONS:
            nominal_b = frac_b * curve_b.ic50_sim
            effective_b = min(1.0, nominal_b * resistance_multiplier_b)
            growth = simulate_growth([(drug_a, effective_a), (drug_b, effective_b)]).growth_fraction
            if growth > INHIBITORY_THRESHOLD:
                continue
            sigma = frac_a + frac_b
            if best is None or sigma < best.sigma_fic:
                best = FICResult(
                    fic_a=frac_a,
                    fic_b=frac_b,
                    sigma_fic=sigma,
                    classification=_classify_fic(sigma),
                    frac_a_used=frac_a,
                    frac_b_used=frac_b,
                )
    return best
