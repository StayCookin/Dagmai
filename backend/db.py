"""Persistence layer. SQLite by default (zero-setup for the MVP); set
DATABASE_URL to a Postgres DSN (e.g. postgresql+psycopg://user:pass@host/db)
to switch -- SQLAlchemy is the only thing that needs to know the difference,
nothing in routers/ or engine/ is SQLite-specific.
"""

from __future__ import annotations

import datetime
import os
import uuid

from sqlalchemy import JSON, DateTime, String, create_engine
from sqlalchemy.orm import DeclarativeBase, Mapped, Session, mapped_column, sessionmaker

DATABASE_URL = os.environ.get("DATABASE_URL", "sqlite:///./amr_combo.db")

_connect_args = {"check_same_thread": False} if DATABASE_URL.startswith("sqlite") else {}
engine = create_engine(DATABASE_URL, connect_args=_connect_args)
SessionLocal = sessionmaker(bind=engine, autoflush=False, autocommit=False)


class Base(DeclarativeBase):
    pass


class SimulationRun(Base):
    __tablename__ = "simulation_runs"

    id: Mapped[str] = mapped_column(String, primary_key=True, default=lambda: str(uuid.uuid4()))
    strain_id: Mapped[str] = mapped_column(String, index=True)
    requested_drug_ids: Mapped[list] = mapped_column(JSON)
    created_at: Mapped[datetime.datetime] = mapped_column(
        DateTime, default=lambda: datetime.datetime.now(datetime.UTC)
    )
    results: Mapped[list] = mapped_column(JSON)


def init_db() -> None:
    Base.metadata.create_all(bind=engine)


def get_session() -> Session:
    return SessionLocal()
