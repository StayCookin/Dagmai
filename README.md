# AMR Drug Combination Simulator (MVP)

Prioritises antimicrobial drug combinations against multidrug-resistant (MDR)
uropathogens by running genome-scale flux balance analysis (FBA) on a
resistance-annotated strain, rather than by guessing from drug class alone.

**Scope, as designed:** one organism (*E. coli*), a fixed 15-drug panel, five
example MDR/ESBL strain genotypes. This mirrors the scoped-down MVP plan in
the original brief -- "one organism, ESBL E. coli, a fixed panel of ~15
drugs... a working pipeline in weeks rather than a stalled general-purpose
platform."

**What this tool is for:** narrowing a large combination space down to a
short, mechanistically-reasoned shortlist worth testing in the lab. It is
**not** a clinical efficacy predictor, and nothing it outputs should inform
patient treatment decisions.

## Quickstart

```bash
python3 -m venv venv && source venv/bin/activate
pip install -r requirements.txt
uvicorn backend.main:app --reload
```

Open `http://127.0.0.1:8000/` for the UI, or `http://127.0.0.1:8000/docs`
for the interactive API docs. First request in a fresh process is slower
(~5-10s) because per-drug dose-response curves and reaction-capacity data get
computed once and cached; every request after that is fast.

Run the test suite:

```bash
pytest tests/ -v
```

## Architecture

```
data/
  drug_panel.json            15-drug panel: target gene, target reaction(s), simulation_mode
  resistance_mechanisms.json 10 hand-curated resistance mechanisms (real gene families,
                              illustrative not literature-sourced multipliers -- see Layer 2)
  strains/demo_strains.json  5 example strain genotypes (mechanism id lists)

backend/
  engine/
    model_loader.py   loads iJO1366 (bundled with cobrapy), computes & caches
                       FVA-derived reaction capacities
    drug_panel.py      loads/indexes data/*.json into typed dataclasses
    fba.py             applies drug(s) as reaction-capacity constraints, runs FBA
    synergy.py          dose-response scan -> IC50_sim, Bliss independence score,
                        checkerboard ΣFIC index
  resistance/
    overlay.py          strain genotype -> per-drug effective-potency multiplier
  reasoning/
    engine.py            rule-based mechanistic explanation generator (Layer 3 MVP)
  routers/               FastAPI endpoints (drugs, strains, simulate)
  schemas.py             Pydantic request/response models
  db.py                   SQLAlchemy persistence (SQLite by default)
  main.py                 FastAPI app, mounts frontend/ as static files

frontend/                 vanilla HTML/CSS/JS single-page UI (no build step)
tests/                    pytest suite
```

### Request flow

1. User picks a strain (a set of resistance mechanism ids) and clicks "rank".
2. `POST /api/simulate/rank` iterates every unordered pair in the 15-drug
   panel (105 pairs). For each pair with two `simulation_mode="direct"`
   drugs, it doses each drug at maximal simulated inhibition strength
   attenuated by the strain's resistance multiplier for that drug (i.e.
   "give the strongest reasonable dose; how much of it actually lands"),
   runs FBA with both constraints applied simultaneously, and computes a
   Bliss independence score. Results are sorted by predicted combined growth
   (ascending -- most suppression first), persisted to SQLite, and returned.
3. Clicking a row calls `GET /api/simulate/pair`, which additionally runs the
   more expensive checkerboard ΣFIC search (a small grid of FBA calls) and
   the rule-based reasoning engine, producing a mechanistic explanation.

## Layer 1: FBA modeling choices (the substance, and its limits)

### Model: iJO1366, not iML1515

The original design called for iML1515. That model is distributed only from
the BiGG Models repository (`bigg.ucsd.edu`); this deployment's network
policy does not reach that host (confirmed: `CONNECT bigg.ucsd.edu:443` is
rejected by the egress proxy). **iJO1366** -- iML1515's direct predecessor,
same *E. coli* K-12 MG1655 lineage, 2583 reactions / 1367 genes -- ships
bundled inside the `cobra` package itself, so it loads with zero network
access. The two models share the same core reaction network; almost every
gene/reaction id used in `data/drug_panel.json` carries over unchanged.
Swapping in iML1515 later, once network access exists, is a one-line change
in `backend/engine/model_loader.py` (`_MODEL_PATH`).

### Representing a drug as a reaction constraint

A `simulation_mode="direct"` drug in the panel targets a specific enzyme;
`data/drug_panel.json` lists the reaction id(s) that enzyme catalyses in
iJO1366. Applying the drug means capping those reactions' flux bounds.

**Why capacity-relative, not bound-relative, scaling.** The naive approach
-- cap a reaction's bound at `(1 - strength) * default_bound` (default bounds
are ±1000) -- produces a step function: growth is unaffected until strength
is within a percent or two of 1.0, because 1000 is ~40,000x larger than the
flux any single essential reaction here actually carries at the optimal
growth rate. That gives a useless, all-or-nothing dose-response curve.
Instead, each target reaction's capacity is anchored to **flux variability
analysis (FVA) at 99.9% of the wild-type growth optimum** (see
`model_loader.py`, `FVA_FRACTION_OF_OPTIMUM`), computed once at startup and
cached. A drug at `strength=s` caps the reaction at
`FVA_capacity * (1 - s) * 1.15` (the 1.15 is headroom so `strength=0` never
nudges growth below the true wild-type optimum from solver tolerance noise).
This produces graded, monotonic dose-response curves for every direct-mode
drug in the panel (verified in `tests/test_fba.py`).

**Composing two drugs on the same reaction.** When a combination's two drugs
both constrain the same reaction (e.g. two beta-lactams both capping the PBP
transpeptidation reactions), the tighter of the two caps applies -- not
whichever drug happened to be applied second. See `fba.py::_apply_drug`.

**`strength` is dimensionless**, not a real drug concentration or MIC. It
is deliberately not anchored to pharmacokinetic data (urinary concentration,
protein binding, half-life), which is out of scope for a metabolic-flux
model. `IC50_sim` (the strength at which growth crosses 50% of wild-type,
found by scanning strength 0..1) stands in for "the reference dose that just
crosses an inhibitory threshold," used to make Bliss and FIC scores
comparable across drugs.

### Bliss independence score and ΣFIC index

`bliss_score()` computes each drug's growth-inhibiting fraction alone and in
combination, then compares the observed combined inhibition to the Bliss
*independence* expectation (`iA + iB - iA*iB`). Positive = synergy,
~0 = additive, negative = antagonism.

`fic_index()` runs a checkerboard-style grid search over fractions of each
drug's own (strain-adjusted) IC50_sim, finds the combination reaching 50%
growth inhibition with the smallest total dose fraction, and classifies it
using the standard clinical ΣFIC bands (≤0.5 synergy, ≤1 additive, ≤4
indifferent, >4 antagonism).

### An honest, empirically-observed limitation: FBA under-detects synergy

Run the numbers on **trimethoprim + sulfamethoxazole** -- co-trimoxazole,
the textbook sequential-blockade synergy pair (they inhibit consecutive
steps of the same folate biosynthesis pathway) -- and this model scores it
as **antagonistic** under Bliss, not synergistic (`bliss_score ≈ -0.25`,
pinned as a regression test in `tests/test_synergy.py`). This isn't a bug.
Two reactions on the same linear pathway share a single bottleneck: in a
steady-state LP, capping both rarely reduces the optimum further than
capping the more restrictive one alone, since the model can't represent the
kinetic, time-resolved depletion dynamics that make sequential blockade
synergistic in real cells. The same pattern shows up even for drugs on
*different* pathways (fosfomycin + colistin also scores non-synergistic),
because FBA's growth optimum is generally governed by whichever single
constraint is tightest, not by how many constraints are simultaneously
active, unless the pathways are metabolically coupled.

**Consequence for how this tool should be read:** Bliss/ΣFIC numbers are
useful context, not the primary ranking signal. The `/api/simulate/rank`
endpoint sorts by **predicted combined growth suppression** (does this pair,
given this strain's resistance profile, actually knock growth down?), which
correctly surfaces "the strain's remaining active drug(s), paired with
something that doesn't get destroyed by the same mechanism" -- a real and
useful screening signal -- without overclaiming synergy detection that
static FBA structurally cannot deliver. See `backend/engine/synergy.py`'s
module docstring for the full derivation.

## Layer 2: resistance overlay

`data/resistance_mechanisms.json` hand-curates 10 mechanisms (ESBL, KPC
carbapenemase, AcrAB-TolC efflux overexpression, GyrA target mutation,
DfrA/Sul1 target bypass, MCR-1, ArmA, FosA, NfsA loss), each naming a real,
well-established resistance gene family (`blaCTX-M-15`, `mcr-1`, `dfrA1`,
etc. -- not fabricated) and its effect on specific panel drugs, expressed as
an "effective potency multiplier" in (0, 1]. A strain
(`data/strains/demo_strains.json`) is just a list of mechanism ids;
`backend/resistance/overlay.py` resolves them into per-drug multipliers,
taking the most severe multiplier when several mechanisms hit the same drug.

**Two different kinds of confidence, kept explicit per mechanism.** The
*ordinal direction* of every multiplier is standard, textbook-level clinical
microbiology and is high confidence: CTX-M-15 reliably inactivates
penicillins/cephalosporins and spares carbapenems; KPC-2 does not spare
carbapenems; dfrA/sul1 are target-bypass mechanisms that near-completely
defeat their drug; and so on. The *exact decimal* within that direction
(0.03 vs 0.05, say) is an illustrative placeholder, not sourced from any
specific fold-MIC-change study -- this environment has no network access to
verify one against the literature. Each mechanism's `magnitude_confidence`
field in the JSON says explicitly which parts are which, rather than
presenting the whole multiplier as equally well-founded.

An earlier version of this file also carried specific numeric CARD/ARO
accession IDs (e.g. `"ARO:3001864"`) next to each mechanism, styled to look
like a live database lookup. They weren't -- this environment has no network
access to CARD either (confirmed: `CONNECT card.mcmaster.ca:443` is
rejected), so those IDs could not be checked and have been removed rather
than left in place implying a verification that never happened.

**What's not implemented (originally scoped as live AMRFinderPlus/RGI
annotation of an uploaded genome):** there's no FASTA upload, no BLAST
database, no live CARD lookup. That pipeline needs a local CARD/RGI install
and strain assembly data neither available nor sensible to fake in this
build. The mechanism list here is a fixed, hand-picked reference set instead
of database-driven strain-specific annotation -- the natural extension point
is swapping `strain_resistance_profile()`'s input (a list of mechanism ids)
for the output of an actual RGI run.

## Layer 3: reasoning (rule-based MVP, not an LLM)

The original design called for Llama 3.3 70B or Qwen 2.5 via Ollama/vLLM,
reading FBA output plus CARD annotations to generate mechanistic rationale,
fine-tuned or few-shot prompted on DrugComb / literature-mined synergy pairs.
That needs either a GPU-hosted model server or the ability to download model
weights from HuggingFace; this sandboxed build environment has neither
(outbound access is restricted to package registries and a small allowlist
-- `bigg.ucsd.edu` and equivalent model-weight hosts are not reachable).

`backend/reasoning/engine.py` implements `RuleBasedReasoningEngine` instead:
deterministic template logic that reads the *same* structured inputs an LLM
would (Bliss/FIC results, resistance mechanism metadata, drug target
annotations) and produces a grounded mechanistic sentence -- not a static
mock, since the output genuinely changes with the underlying numbers. A
small hand-curated `LITERATURE_SYNERGY_NOTES` table stands in for what
ChemBERTa/ESM-2 embeddings plus DrugComb fine-tuning would otherwise surface
for pairs with a non-metabolic drug (fluoroquinolones, aminoglycosides,
nitrofurantoin) that FBA cannot score.

`ReasoningEngine` is an abstract base class specifically so a real
LLM-backed implementation can be swapped in later (`LLMReasoningEngine`,
prompting Llama/Qwen with the same structured inputs) without touching
routers, schemas, or the frontend.

## Layer 4: app

FastAPI + SQLAlchemy (SQLite by default; set `DATABASE_URL` to a Postgres
DSN to switch, e.g. `postgresql+psycopg://user:pass@host/db` -- nothing
outside `backend/db.py` is SQLite-specific) + a plain HTML/CSS/JS frontend
(no build tooling, served directly as static files by FastAPI). A React
frontend is the natural next step for a production UI; vanilla JS was chosen
here to ship a working MVP without adding Node build complexity.

Endpoints:
- `GET /api/drugs`, `GET /api/drugs/{id}`, `GET /api/drugs/{id}/dose-response`
- `GET /api/drugs/meta/resistance-mechanisms`
- `GET /api/strains`, `GET /api/strains/{id}`
- `POST /api/simulate/rank` `{strain_id, drug_ids?}` -> ranked combinations (Bliss only, fast)
- `GET /api/simulate/pair?strain_id&drug_a&drug_b` -> full detail (Bliss + ΣFIC + reasoning)

## The panel: what's FBA-simulable and what isn't

7 of the 15 drugs are `simulation_mode="direct"` (fosfomycin, trimethoprim,
sulfamethoxazole, 5 beta-lactams sharing the PBP reaction set, and colistin
via the lipid A pathway). The other 8 (nitrofurantoin, ciprofloxacin,
levofloxacin, gentamicin, amikacin) target processes -- DNA gyrase/
topoisomerase IV, the ribosome, multi-target prodrug activation -- that
genuinely have no representation in this genome-scale *metabolic*
reconstruction, because none of those processes are metabolic reactions with
defined stoichiometry; no published bacterial GEM (this one included)
represents DNA replication, transcription, or translation machinery as model
reactions, only as a lumped biomass precursor demand.

This is checked structurally, not by trusting a specific gene locus tag:
`tests/test_fba.py::test_dna_replication_transcription_translation_machinery_absent_from_model`
searches every reaction name in the loaded model for DNA/RNA-polymerase/
ribosome/translation/replication/gyrase/topoisomerase keywords and asserts
zero matches. An earlier version of this codebase instead hardcoded specific
b-numbers (`gyrA = b2231`, etc.) and asserted those were absent -- which
turned out to be the wrong kind of check: a memorized locus tag can simply be
wrong, and a version of exactly that mistake was caught for a different drug
during a correctness pass (see "A verification pass caught two real errors"
below). The keyword search doesn't depend on getting any specific locus tag
right.

These 8 drugs are marked `non_metabolic` rather than silently given a fake
score, and the UI/reasoning layer clearly labels them as literature-only or
insufficient-data rather than presenting a number that doesn't exist.

## A verification pass caught two real errors

This codebase's gene/reaction mappings were built from training-data recall
of *E. coli* locus tags, not a live database lookup (no network access to
NCBI/EcoCyc/UniProt from this environment -- see below). Recall of an exact
b-number is more fallible than recall of a well-known biological fact, so
before trusting any drug-target mapping used for a real computation, each one
was cross-checked against the loaded model itself: does the gene the mapping
names actually catalyse the reaction the enzyme is known for? This caught
two real errors, both now fixed and pinned with regression tests:

1. **Nitrofurantoin's target genes were wrong.** The panel originally listed
   `b1636` as `nfsA`. It isn't -- `b1636` is `pdxY` (pyridoxal kinase), a
   completely unrelated gene; its single reaction (`PYDXK`) is nowhere near
   nitroreductase chemistry. Re-checked via UniProt accession cross-reference
   instead of b-number memory: `nfsB` (P38489) *is* in the model, as `b0578`
   -- but the only reactions it carries here (`DHPTDNR`/`DHPTDNRN`,
   dihydropteridine reductase) are a documented NfsB moonlighting activity
   unrelated to nitrofuran reduction. No nitrofuran-activation reaction
   exists anywhere in the model under any gene (confirmed by keyword search).
   `nfsA`'s real locus tag couldn't be independently confirmed from this
   environment, so `data/drug_panel.json` now leaves it as a named-but-unpinned
   gene rather than guessing again. Net effect on results: none --
   nitrofurantoin was and remains `non_metabolic`, so this only affected the
   informational gene-name field, not any computation.

2. **Colistin's target reaction was wrong, and this one *did* affect
   results.** The panel used `b0180` with reactions `OGMEACPD`/`OPMEACPD` as
   an "LpxC" proxy for lipid A pathway dependency. `b0180`'s actual reactions
   are `3-hydroxyacyl-[ACP] dehydratase` across many acyl-chain lengths --
   classic **FabZ** (fatty acid synthesis) substrate promiscuity, not LpxC's
   deacetylase chemistry. Likely cause: pattern-matching on similar-sounding
   reaction names ("hydroxyacyl...dehydratase" vs LpxC's real
   "hydroxymyristoyl...deacetylase") without checking the actual enzyme
   class. Found the real LpxC reaction by searching the model for lipid-A
   pathway keywords: `UHGADA` ("UDP-3-O-acetylglucosamine deacetylase"),
   gene `b0096` -- single-gene, single-reaction, confirmed essential (full
   knockout drives growth to exactly zero), and a much better-conditioned
   target than the original choice, whose tiny baseline flux (~2e-6, several
   orders of magnitude smaller than every other target reaction in the
   panel) turned out to also be the trigger for a GLPK solver pathology (see
   `model_loader.py`'s `FAST_TIMEOUT_SECONDS` comment) -- fixing the science
   here incidentally fixed a performance bug too. `data/drug_panel.json` and
   colistin's dose-response numbers throughout the app reflect the corrected
   reaction as of this pass.

**What's still unverified.** This environment has no network access to
NCBI/EcoCyc/UniProt/CARD (confirmed: `CONNECT` to each is rejected by the
egress proxy's organization policy), so every gene/reaction mapping in
`data/drug_panel.json` was checked the same way -- against the loaded
model's own reaction chemistry, not an external database -- and every one
that catalyses its textbook-expected reaction (fosfomycin/MurA/`UAGCVT`,
trimethoprim/FolA/`DHFR`, sulfamethoxazole/FolP/`DHPS2`, the PBP family
genes and their murein-crosslinking/transglycosylation reactions,
colistin/LpxC/`UHGADA` after the fix above) is reported with that
corroboration explicitly in mind. The non-metabolic drugs' gene names
(`gyrA/gyrB`, `16S rRNA`/ribosomal proteins) describe well-known
*functional* targets with high confidence; their specific locus tags, where
still shown, have not been independently re-verified and should not be
treated as citations.

## Known limitations (beyond the FBA-synergy finding above)

- **Fungal targets are entirely out of scope for this MVP** (AGORA/iMM904
  models, the original Layer 1 spec's fungal option, are not wired in).
- **No PK/PD layer.** Urinary concentration, protein binding, and half-life
  are not modeled; `strength`/IC50_sim are dimensionless model quantities,
  not real doses.
- **No biofilm, membrane-permeability, or efflux-kinetics modeling** beyond
  the coarse multiplicative resistance-overlay effect.
- **Beta-lactam subclasses share one reaction set.** Genome-scale metabolic
  models don't resolve PBP-binding-affinity differences between e.g.
  ceftriaxone and meropenem; they're differentiated here almost entirely
  through the resistance overlay (which beta-lactamases hydrolyse which
  drug), which is the more clinically important axis for MDR anyway, but is
  a real simplification worth knowing about.
- **Resistance mechanism data is hand-curated, not database-driven** (see
  Layer 2 above).

## What this tool defensibly claims to do

Narrow a large combination space to a short, mechanistically-explainable
shortlist worth testing in the lab -- not predict clinical efficacy. Every
number the UI shows traces back to a real FBA computation or is explicitly
labeled as a literature note or non-simulable, rather than presenting a
confident score for something the model can't actually see.
