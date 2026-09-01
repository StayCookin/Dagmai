"""Loads and indexes the fixed drug panel and resistance mechanism reference data."""

from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from functools import lru_cache

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data")
DRUG_PANEL_PATH = os.path.join(DATA_DIR, "drug_panel.json")
RESISTANCE_PATH = os.path.join(DATA_DIR, "resistance_mechanisms.json")
STRAINS_PATH = os.path.join(DATA_DIR, "strains", "demo_strains.json")


@dataclass(frozen=True)
class Drug:
    id: str
    display_name: str
    drug_class: str
    target_gene: str
    target_gene_name: str
    target_enzyme: str
    pathway: str
    reactions: tuple[str, ...]
    simulation_mode: str  # "direct" | "non_metabolic"
    notes: str


@dataclass(frozen=True)
class ResistanceMechanism:
    id: str
    name: str
    gene_marker: str
    effect_type: str  # "remove_inhibition" | "scale_efflux" | "shift_constraint"
    affected_drugs: dict[str, float]  # drug_id -> effective strength multiplier
    magnitude_confidence: str  # what's literature-grounded vs illustrative about the multipliers
    description: str


@dataclass(frozen=True)
class Strain:
    id: str
    display_name: str
    organism: str
    source_note: str
    mechanisms: tuple[str, ...] = field(default_factory=tuple)


@lru_cache(maxsize=1)
def load_drugs() -> dict[str, Drug]:
    with open(DRUG_PANEL_PATH) as f:
        raw = json.load(f)
    drugs = {}
    for d in raw["drugs"]:
        drugs[d["id"]] = Drug(
            id=d["id"],
            display_name=d["display_name"],
            drug_class=d["class"],
            target_gene=d["target_gene"],
            target_gene_name=d["target_gene_name"],
            target_enzyme=d["target_enzyme"],
            pathway=d["pathway"],
            reactions=tuple(d["reactions"]),
            simulation_mode=d["simulation_mode"],
            notes=d["notes"],
        )
    return drugs


@lru_cache(maxsize=1)
def load_resistance_mechanisms() -> dict[str, ResistanceMechanism]:
    with open(RESISTANCE_PATH) as f:
        raw = json.load(f)
    mechanisms = {}
    for m in raw["mechanisms"]:
        mechanisms[m["id"]] = ResistanceMechanism(
            id=m["id"],
            name=m["name"],
            gene_marker=m["gene_marker"],
            effect_type=m["effect_type"],
            affected_drugs=dict(m["affected_drugs"]),
            magnitude_confidence=m["magnitude_confidence"],
            description=m["description"],
        )
    return mechanisms


@lru_cache(maxsize=1)
def load_strains() -> dict[str, Strain]:
    with open(STRAINS_PATH) as f:
        raw = json.load(f)
    strains = {}
    for s in raw["strains"]:
        strains[s["id"]] = Strain(
            id=s["id"],
            display_name=s["display_name"],
            organism=s["organism"],
            source_note=s["source_note"],
            mechanisms=tuple(s["mechanisms"]),
        )
    return strains


def all_target_reactions() -> tuple[str, ...]:
    """Union of every reaction referenced by any 'direct' drug in the panel."""
    reactions: set[str] = set()
    for drug in load_drugs().values():
        reactions.update(drug.reactions)
    return tuple(sorted(reactions))
