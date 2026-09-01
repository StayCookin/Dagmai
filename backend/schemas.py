"""Pydantic request/response models for the API."""

from __future__ import annotations

from pydantic import BaseModel


class DrugOut(BaseModel):
    id: str
    display_name: str
    drug_class: str
    target_gene_name: str
    target_enzyme: str
    pathway: str
    simulation_mode: str
    notes: str


class StrainOut(BaseModel):
    id: str
    display_name: str
    organism: str
    source_note: str
    mechanisms: list[str]


class ResistanceMechanismOut(BaseModel):
    id: str
    name: str
    gene_marker: str
    effect_type: str
    magnitude_confidence: str
    description: str


class DoseResponsePointOut(BaseModel):
    strength: float
    growth_fraction: float


class DrugDoseResponseOut(BaseModel):
    drug_id: str
    simulation_mode: str
    points: list[DoseResponsePointOut]
    ic50_sim: float | None


class ResistanceEffectOut(BaseModel):
    multiplier: float
    contributing_mechanisms: list[str]
    mechanism_names: list[str]


class BlissOut(BaseModel):
    growth_fraction_a_alone: float
    growth_fraction_b_alone: float
    growth_fraction_combo: float
    expected_independent_inhibition: float
    observed_inhibition: float
    bliss_score: float
    classification: str


class FICOut(BaseModel):
    fic_a: float
    fic_b: float
    sigma_fic: float
    classification: str


class ExplanationOut(BaseModel):
    headline: str
    rationale: str
    confidence: str
    caveats: list[str]


class CombinationResultOut(BaseModel):
    drug_a: DrugOut
    drug_b: DrugOut
    resistance_a: ResistanceEffectOut | None
    resistance_b: ResistanceEffectOut | None
    bliss: BlissOut | None
    fic: FICOut | None
    explanation: ExplanationOut


class RankRequest(BaseModel):
    strain_id: str
    drug_ids: list[str] | None = None  # None = use full panel


class RankResponseOut(BaseModel):
    strain: StrainOut
    results: list[CombinationResultOut]
