from backend.engine.drug_panel import load_resistance_mechanisms, load_strains
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


def test_resistance_mechanisms_document_their_own_confidence_and_carry_no_fabricated_citations():
    """Every mechanism must say what's literature-grounded (the ordinal
    direction) vs illustrative (the exact decimal) about its multipliers --
    see data/resistance_mechanisms.json's top-level description. This also
    guards against a regression back to specific CARD/ARO accession numbers
    that were removed because they were never checked against a live CARD
    instance and would misrepresent hand-picked values as verified lookups."""
    for mechanism in load_resistance_mechanisms().values():
        assert mechanism.magnitude_confidence
        assert not hasattr(mechanism, "aro_accession")
