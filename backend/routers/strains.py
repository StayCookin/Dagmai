from fastapi import APIRouter, HTTPException

from ..engine.drug_panel import load_strains
from ..schemas import StrainOut

router = APIRouter(prefix="/api/strains", tags=["strains"])


def _strain_out(strain) -> StrainOut:
    return StrainOut(
        id=strain.id,
        display_name=strain.display_name,
        organism=strain.organism,
        source_note=strain.source_note,
        mechanisms=list(strain.mechanisms),
    )


@router.get("", response_model=list[StrainOut])
def list_strains() -> list[StrainOut]:
    return [_strain_out(s) for s in load_strains().values()]


@router.get("/{strain_id}", response_model=StrainOut)
def get_strain(strain_id: str) -> StrainOut:
    strains = load_strains()
    if strain_id not in strains:
        raise HTTPException(status_code=404, detail=f"Unknown strain '{strain_id}'")
    return _strain_out(strains[strain_id])
