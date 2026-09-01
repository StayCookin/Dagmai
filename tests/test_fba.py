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


def test_dna_replication_transcription_translation_machinery_absent_from_model():
    """Confirms the modeling assumption behind marking ciprofloxacin,
    levofloxacin, gentamicin and amikacin non_metabolic: this is a genome-scale
    *metabolic* reconstruction, and DNA gyrase / RNA polymerase / the ribosome
    are not metabolic enzymes with a defined stoichiometric reaction, so no
    published GEM (this one included) represents them.

    Deliberately does NOT hardcode specific gene locus tags (e.g. "gyrA is
    b2231") and check those are absent -- an earlier version of this codebase
    did exactly that for a different drug (nfsA/nfsB, see git history) and got
    the locus tag wrong, silently "confirming absence" of the wrong gene. A
    keyword search over every reaction's name is robust to that failure mode:
    it doesn't matter what gyrA's real b-number is, because if DNA gyrase
    supercoiling were represented as a model reaction at all, it would show up
    under some name containing one of these terms, and none do.
    """
    from backend.engine.model_loader import _load_base_model

    model = _load_base_model()
    keywords = ["dna ", "rna polymerase", "ribosom", "translation", "replicat", "transcript", "gyrase", "topoisomerase"]
    matches = [r.name for r in model.reactions if r.name and any(k in r.name.lower() for k in keywords)]
    assert matches == []


def test_ijo1366_gene_and_reaction_counts_match_published_model():
    """Sanity check that the correct, well-documented model is loaded (Orth
    et al. 2011: iJO1366 has 1366 genes -- 1367 including the spontaneous
    reaction placeholder gene -- 2583 reactions, 1805 metabolites). If this
    ever fails, something is wrong with the bundled model file itself, not
    with this codebase's gene/reaction mapping choices."""
    from backend.engine.model_loader import _load_base_model

    model = _load_base_model()
    assert len(model.genes) == 1367
    assert len(model.reactions) == 2583
    assert len(model.metabolites) == 1805


def test_nfsb_present_but_only_via_unrelated_moonlighting_reaction():
    """Locks in a verified finding behind nitrofurantoin's non_metabolic
    classification: nfsB *is* a gene in this model (b0578), so "the gene is
    just absent" would be the wrong reason to give. What's actually true is
    narrower and was checked directly against the loaded model: nfsB's only
    reactions here (DHPTDNR/DHPTDNRN, dihydropteridine reductase) are an
    unrelated documented NfsB moonlighting activity, not nitrofuran
    reduction -- and no nitrofuran-activation reaction exists anywhere in the
    model under any gene."""
    from backend.engine.model_loader import _load_base_model

    model = _load_base_model()
    assert "b0578" in model.genes
    reaction_ids = {r.id for r in model.genes.get_by_id("b0578").reactions}
    assert reaction_ids == {"DHPTDNR", "DHPTDNRN"}
    # Broader "nitro"/"furan" substrings are too loose here -- they also match
    # unrelated reactions like "Nitrous oxide exchange" and
    # "...tetrahydrofuran synthesis", so check for the compound term instead.
    assert not any("nitrofuran" in (r.name or "").lower() or "nitroreductase" in (r.name or "").lower() for r in model.reactions)


def test_colistin_target_reaction_is_essential_and_single_gene():
    """UHGADA/b0096 (the corrected LpxC mapping -- see drug_panel.json notes
    for the earlier b0180/fabZ mix-up this replaced) should behave like a
    clean drug target: one gene, one reaction, essential for growth."""
    from backend.engine.model_loader import _load_base_model
    from cobra.manipulation.delete import knock_out_model_genes

    model = _load_base_model()
    gene = model.genes.get_by_id("b0096")
    assert {r.id for r in gene.reactions} == {"UHGADA"}
    with model:
        knock_out_model_genes(model, [gene])
        assert model.slim_optimize() == pytest.approx(0.0, abs=1e-6)
