import os

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles

from .db import init_db
from .routers import drugs, simulate, strains

app = FastAPI(
    title="AMR Drug Combination Simulator",
    description=(
        "Prioritises antimicrobial drug combinations against MDR uropathogens using "
        "genome-scale FBA (Layer 1), a strain resistance overlay (Layer 2), and a "
        "rule-based mechanistic reasoning layer (Layer 3 MVP). See README for scope "
        "and honest constraints."
    ),
    version="0.1.0",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(drugs.router)
app.include_router(strains.router)
app.include_router(simulate.router)


@app.on_event("startup")
def _startup() -> None:
    init_db()


_FRONTEND_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "frontend")
if os.path.isdir(_FRONTEND_DIR):
    app.mount("/", StaticFiles(directory=_FRONTEND_DIR, html=True), name="frontend")
