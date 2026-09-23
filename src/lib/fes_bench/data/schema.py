"""Schema constants and validation helpers for benchmark data."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any

REFERENCE_COLUMNS = ("T_K", "P_GPa", "G_eV_per_atom", "level", "source")
REQUIRED_META_FIELDS = (
    "functional",
    "dispersion",
    "supercell",
    "kpoints",
    "convergence",
    "size_error_eV_per_atom",
    "method",
    "doi",
    "notes",
)
REQUIRED_SYSTEM_FIELDS = ("phases", "type_map", "reference_grid", "truth_level")


class SchemaError(ValueError):
    """Raised when an on-disk benchmark artifact breaks the canonical schema."""


def require_mapping(value: Any, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise SchemaError(f"{label} must be a JSON object")
    return value


def require_fields(value: Mapping[str, Any], fields: Sequence[str], label: str) -> None:
    missing = [field for field in fields if field not in value]
    if missing:
        raise SchemaError(f"{label} is missing required fields: {', '.join(missing)}")


def finite_float(value: Any, label: str) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise SchemaError(f"{label} must be numeric, got {value!r}") from exc
    if not math.isfinite(parsed):
        raise SchemaError(f"{label} must be finite, got {value!r}")
    return parsed


def validate_system_meta(value: Any, system: str) -> Mapping[str, Any]:
    meta = require_mapping(value, f"data/{system}/system.json")
    require_fields(meta, REQUIRED_SYSTEM_FIELDS, f"data/{system}/system.json")
    phases = meta["phases"]
    if not isinstance(phases, list) or not all(isinstance(item, str) and item for item in phases):
        raise SchemaError(f"data/{system}/system.json field 'phases' must be a non-empty list of names")
    if not isinstance(meta["type_map"], list) or not meta["type_map"]:
        raise SchemaError(f"data/{system}/system.json field 'type_map' must be a non-empty list")
    if not isinstance(meta["reference_grid"], Mapping):
        raise SchemaError(f"data/{system}/system.json field 'reference_grid' must be an object")
    if not isinstance(meta["truth_level"], str) or not meta["truth_level"]:
        raise SchemaError(f"data/{system}/system.json field 'truth_level' must be a non-empty string")
    return meta


def validate_phase_meta(value: Any, label: str) -> Mapping[str, Any]:
    meta = require_mapping(value, label)
    require_fields(meta, REQUIRED_META_FIELDS, label)
    finite_float(meta["size_error_eV_per_atom"], f"{label}.size_error_eV_per_atom")
    if not isinstance(meta["method"], str) or not meta["method"]:
        raise SchemaError(f"{label}.method must be a non-empty string")
    return meta
