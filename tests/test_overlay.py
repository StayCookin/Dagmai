from backend.engine.drug_panel import load_strains
from backend.resistance.overlay import effective_multiplier, strain_resistance_profile


def test_wild_type_strain_has_empty_profile():
    strain = load_strains()["wt_susceptible"]
    profile = strain_resistance_profile(strain)
    assert profile == {}


def test_esbl_strain_cripples_ceftriaxone_but_spares_meropenem():
    strain = load_strains()["esbl_ecoli_uti"]
    profile = strain_resistance_profile(strain)
    assert effective_multiplier("ceftriaxone", profile) < 0.1
    assert effective_multiplier("meropenem", profile) > 0.85


def test_unaffected_drug_defaults_to_full_multiplier():
    strain = load_strains()["esbl_ecoli_uti"]
    profile = strain_resistance_profile(strain)
    # Colistin is not a substrate of any mechanism this strain carries.
    assert effective_multiplier("colistin", profile) == 1.0


def test_multiple_mechanisms_take_the_more_severe_multiplier():
    strain = load_strains()["pan_resistant_kpc"]
    profile = strain_resistance_profile(strain)
    # KPC alone would give ceftriaxone ~0.04; ESBL alone also ~0.04-0.05.
    # Either way it should stay near-fully neutralized, not compounded below
    # either individual mechanism's own floor.
    mult = effective_multiplier("ceftriaxone", profile)
    assert 0.0 <= mult <= 0.06
