import pytest

from backend.engine.drug_panel import load_drugs
from backend.engine.synergy import bliss_score, dose_response_curve, fic_index


@pytest.fixture(scope="module")
def drugs():
    return load_drugs()


def test_direct_drug_has_ic50_sim():
    curve = dose_response_curve("sulfamethoxazole")
    assert curve.ic50_sim is not None
    assert 0 < curve.ic50_sim < 1


def test_non_metabolic_drug_has_no_ic50_sim():
    curve = dose_response_curve("gentamicin")
    assert curve.ic50_sim is None
    assert curve.points == ()


def test_bliss_score_none_when_either_drug_non_metabolic(drugs):
    result = bliss_score(drugs["ciprofloxacin"], 1.0, drugs["fosfomycin"], 1.0)
    assert result is None


def test_bliss_score_computed_for_two_direct_drugs(drugs):
    result = bliss_score(drugs["fosfomycin"], 0.9, drugs["colistin"], 0.9)
    assert result is not None
    assert 0 <= result.growth_fraction_combo <= 1
    assert result.classification in {"synergy", "additive", "antagonism"}


def test_same_pathway_sequential_blockade_does_not_score_as_synergy(drugs):
    """Documented finding (see synergy.py module docstring / README): two
    reactions on the same linear pathway (trimethoprim/DHFR then
    sulfamethoxazole/DHPS, both folate biosynthesis) reduce combined growth to
    roughly the more restrictive single constraint under steady-state FBA, so
    this combo should NOT show a positive (synergistic) Bliss score here even
    though it is the textbook clinical synergy pair -- that gap is the point,
    not a bug, and this test pins the behavior so it doesn't silently change."""
    tri, sul = drugs["trimethoprim"], drugs["sulfamethoxazole"]
    tri_curve, sul_curve = dose_response_curve("trimethoprim"), dose_response_curve("sulfamethoxazole")
    result = bliss_score(tri, tri_curve.ic50_sim, sul, sul_curve.ic50_sim)
    assert result.bliss_score <= 0.05


def test_fic_index_none_when_either_drug_non_metabolic(drugs):
    assert fic_index(drugs["ciprofloxacin"], drugs["fosfomycin"]) is None


def test_fic_index_returns_sane_classification(drugs):
    result = fic_index(drugs["fosfomycin"], drugs["colistin"])
    assert result is not None
    assert result.classification in {"synergy", "additive", "indifferent", "antagonism"}
    assert result.sigma_fic > 0
