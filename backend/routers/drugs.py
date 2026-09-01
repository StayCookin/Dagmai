from fastapi import APIRouter, HTTPException

from ..engine.drug_panel import load_drugs, load_resistance_mechanisms
from ..engine.synergy import dose_response_curve
from ..schemas import DoseResponsePointOut, DrugDoseResponseOut, DrugOut, ResistanceMechanismOut

router = APIRouter(prefix="/api/drugs", tags=["drugs"])


def _drug_out(drug) -> DrugOut:
    return DrugOut(
        id=drug.id,
        display_name=drug.display_name,
        drug_class=drug.drug_class,
        target_gene_name=drug.target_gene_name,
        target_enzyme=drug.target_enzyme,
        pathway=drug.pathway,
        simulation_mode=drug.simulation_mode,
        notes=drug.notes,
    )


@router.get("", response_model=list[DrugOut])
def list_drugs() -> list[DrugOut]:
    return [_drug_out(d) for d in load_drugs().values()]


@router.get("/{drug_id}", response_model=DrugOut)
def get_drug(drug_id: str) -> DrugOut:
    drugs = load_drugs()
    if drug_id not in drugs:
        raise HTTPException(status_code=404, detail=f"Unknown drug '{drug_id}'")
    return _drug_out(drugs[drug_id])


@router.get("/{drug_id}/dose-response", response_model=DrugDoseResponseOut)
def get_dose_response(drug_id: str) -> DrugDoseResponseOut:
    drugs = load_drugs()
    if drug_id not in drugs:
        raise HTTPException(status_code=404, detail=f"Unknown drug '{drug_id}'")
    curve = dose_response_curve(drug_id)
    return DrugDoseResponseOut(
        drug_id=drug_id,
        simulation_mode=drugs[drug_id].simulation_mode,
        points=[DoseResponsePointOut(strength=p.strength, growth_fraction=p.growth_fraction) for p in curve.points],
        ic50_sim=curve.ic50_sim,
    )


@router.get("/meta/resistance-mechanisms", response_model=list[ResistanceMechanismOut])
def list_resistance_mechanisms() -> list[ResistanceMechanismOut]:
    return [
        ResistanceMechanismOut(
            id=m.id,
            name=m.name,
            gene_marker=m.gene_marker,
            effect_type=m.effect_type,
            magnitude_confidence=m.magnitude_confidence,
            description=m.description,
        )
        for m in load_resistance_mechanisms().values()
    ]
