"""Read the canonical per-system/per-phase free-energy dataset layout."""

from __future__ import annotations

import csv
import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .schema import (
    REFERENCE_COLUMNS,
    SchemaError,
    finite_float,
    validate_phase_meta,
    validate_system_meta,
)


@dataclass(frozen=True)
class ReferencePoint:
    """One normalized free-energy reference point (K, GPa, eV/atom)."""

    T_K: float
    P_GPa: float
    G_eV_per_atom: float
    level: str
    source: str


@dataclass(frozen=True)
class Phase:
    """A validated benchmark phase and its provenance."""

    system: str
    name: str
    structure: Path
    G_table: tuple[ReferencePoint, ...]
    meta: dict[str, Any]
    system_meta: dict[str, Any]


def default_data_root() -> Path:
    configured = os.environ.get("FES_BENCH_DATA_ROOT")
    if configured:
        return Path(configured).expanduser().resolve()
    return Path(__file__).resolve().parents[2] / "data"


def _read_json(path: Path) -> dict[str, Any]:
    try:
        parsed = json.loads(path.read_text(encoding="utf-8"))
    except OSError as exc:
        raise SchemaError(f"cannot read {path}: {exc}") from exc
    except json.JSONDecodeError as exc:
        raise SchemaError(f"invalid JSON in {path}: {exc}") from exc
    if not isinstance(parsed, dict):
        raise SchemaError(f"{path} must contain a JSON object")
    return parsed


def _read_reference_table(path: Path) -> tuple[ReferencePoint, ...]:
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            fieldnames = tuple(reader.fieldnames or ())
            if fieldnames != REFERENCE_COLUMNS:
                raise SchemaError(
                    f"{path} must have exactly these columns in order: "
                    f"{', '.join(REFERENCE_COLUMNS)}"
                )
            rows: list[ReferencePoint] = []
            for index, row in enumerate(reader, start=2):
                level = (row["level"] or "").strip()
                source = (row["source"] or "").strip()
                if not level or not source:
                    raise SchemaError(f"{path}:{index} requires non-empty level and source")
                rows.append(
                    ReferencePoint(
                        T_K=finite_float(row["T_K"], f"{path}:{index}:T_K"),
                        P_GPa=finite_float(row["P_GPa"], f"{path}:{index}:P_GPa"),
                        G_eV_per_atom=finite_float(
                            row["G_eV_per_atom"], f"{path}:{index}:G_eV_per_atom"
                        ),
                        level=level,
                        source=source,
                    )
                )
    except OSError as exc:
        raise SchemaError(f"cannot read {path}: {exc}") from exc
    if not rows:
        raise SchemaError(f"{path} must contain at least one reference row")
    return tuple(rows)


def load(system: str, phase: str, data_root: str | Path | None = None) -> Phase:
    """Load one phase from a canonical data directory.

    Parameters have no implicit system-specific mapping: the names must match
    the directories and `system.json` phase list exactly.
    """

    if not system or Path(system).name != system:
        raise SchemaError("system must be a single non-empty directory name")
    if not phase or Path(phase).name != phase:
        raise SchemaError("phase must be a single non-empty directory name")

    root = Path(data_root).expanduser().resolve() if data_root is not None else default_data_root()
    system_dir = root / system
    phase_dir = system_dir / phase
    system_meta = _read_json(system_dir / "system.json")
    validate_system_meta(system_meta, system)
    if phase not in system_meta["phases"]:
        raise SchemaError(f"phase {phase!r} is not declared by {system_dir / 'system.json'}")

    structure = phase_dir / "structure.extxyz"
    if not structure.is_file():
        raise SchemaError(f"missing representative structure: {structure}")
    if structure.stat().st_size == 0:
        raise SchemaError(f"representative structure is empty: {structure}")

    meta = _read_json(phase_dir / "meta.json")
    validate_phase_meta(meta, str(phase_dir / "meta.json"))
    table = _read_reference_table(phase_dir / "reference_G.csv")
    return Phase(
        system=system,
        name=phase,
        structure=structure,
        G_table=table,
        meta=meta,
        system_meta=system_meta,
    )
