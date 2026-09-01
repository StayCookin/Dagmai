from itertools import combinations

from fastapi import APIRouter, HTTPException

from ..db import SimulationRun, get_session
from ..engine.drug_panel import load_drugs, load_strains
from ..resistance.overlay import strain_resistance_profile
from ..schemas import RankRequest, RankResponseOut, StrainOut
from ..simulate_service import build_combination_result, ranking_sort_key

router = APIRouter(prefix="/api/simulate", tags=["simulate"])


def _strain_out(strain) -> StrainOut:
    return StrainOut(
        id=strain.id,
        display_name=strain.display_name,
        organism=strain.organism,
        source_note=strain.source_note,
        mechanisms=list(strain.mechanisms),
    )


@router.post("/rank", response_model=RankResponseOut)
def rank_combinations(request: RankRequest) -> RankResponseOut:
    strains = load_strains()
    drugs = load_drugs()

    if request.strain_id not in strains:
        raise HTTPException(status_code=404, detail=f"Unknown strain '{request.strain_id}'")
    strain = strains[request.strain_id]

    if request.drug_ids:
        unknown = [d for d in request.drug_ids if d not in drugs]
        if unknown:
            raise HTTPException(status_code=404, detail=f"Unknown drug id(s): {unknown}")
        drug_ids = request.drug_ids
    else:
        drug_ids = list(drugs.keys())

    if len(drug_ids) < 2:
        raise HTTPException(status_code=400, detail="Need at least 2 drugs to form a combination")

    profile = strain_resistance_profile(strain)

    results = [
        build_combination_result(drugs[a], drugs[b], strain, resistance_profile=profile, include_fic=False)
        for a, b in combinations(drug_ids, 2)
    ]
    results.sort(key=ranking_sort_key)

    session = get_session()
    try:
        run = SimulationRun(
            strain_id=strain.id,
            requested_drug_ids=drug_ids,
            results=[r.model_dump() for r in results],
        )
        session.add(run)
        session.commit()
    finally:
        session.close()

    return RankResponseOut(strain=_strain_out(strain), results=results)


@router.get("/pair")
def simulate_pair(strain_id: str, drug_a: str, drug_b: str):
    strains = load_strains()
    drugs = load_drugs()

    if strain_id not in strains:
        raise HTTPException(status_code=404, detail=f"Unknown strain '{strain_id}'")
    if drug_a not in drugs:
        raise HTTPException(status_code=404, detail=f"Unknown drug '{drug_a}'")
    if drug_b not in drugs:
        raise HTTPException(status_code=404, detail=f"Unknown drug '{drug_b}'")
    if drug_a == drug_b:
        raise HTTPException(status_code=400, detail="drug_a and drug_b must be different")

    strain = strains[strain_id]
    result = build_combination_result(drugs[drug_a], drugs[drug_b], strain, include_fic=True)
    return result
