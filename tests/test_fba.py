import pytest

from backend.engine.drug_panel import load_drugs
from backend.engine.fba import simulate_growth, wild_type_growth


@pytest.fixture(scope="module")
def drugs():
    return load_drugs()


def test_wild_type_growth_is_positive():
    assert wild_type_growth() > 0.5


def test_no_drug_gives_full_growth(drugs):
    result = simulate_growth([])
    assert result.growth_fraction == pytest.approx(1.0, abs=1e-6)


def test_essential_direct_drug_at_full_strength_kills_growth(drugs):
    fosfomycin = drugs["fosfomycin"]
    result = simulate_growth([(fosfomycin, 1.0)])
    assert result.growth_fraction < 0.01


def test_dose_response_is_monotonic_non_increasing(drugs):
    sulfamethoxazole = drugs["sulfamethoxazole"]
    fractions = [simulate_growth([(sulfamethoxazole, s / 10)]).growth_fraction for s in range(11)]
    for a, b in zip(fractions, fractions[1:]):
        assert a >= b - 1e-9


def test_non_metabolic_drug_has_no_effect(drugs):
    ciprofloxacin = drugs["ciprofloxacin"]
    assert ciprofloxacin.simulation_mode == "non_metabolic"
    result = simulate_growth([(ciprofloxacin, 1.0)])
    assert result.growth_fraction == pytest.approx(1.0, abs=1e-6)


def test_gyrase_ribosome_rna_pol_genes_absent_from_model():
    """Confirms the modeling assumption documented in drug_panel.json / README:
    the reconstruction genuinely has no gene entries for these targets, which
    is *why* the corresponding drugs are marked non_metabolic rather than an
    arbitrary choice."""
    from backend.engine.model_loader import _load_base_model

    model = _load_base_model()
    for gene_id in ["b2231", "b3699", "b3342", "b3987"]:
        assert gene_id not in model.genes
