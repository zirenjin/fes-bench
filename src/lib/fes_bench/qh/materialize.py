"""Materialize an existing QH diagnostic into the benchmark F_QH schema.

This command is intentionally a provenance-preserving conversion, not a QH
calculation.  It is useful when a prior diagnostic run predates the canonical
``data/<system>/<phase>/fqh.csv`` schema.  A diagnostic that failed its
physical gate remains explicitly unusable downstream through
``phonon_report.json``; this module never infers or upgrades reliability.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path
from typing import Any

from fes_bench.config import ConfigError, load_mapping


FQH_FIELDS = (
    "T_K",
    "P_GPa",
    "V_min_A3_per_atom",
    "F_vib_eV_per_atom",
    "F_QH_eV_per_atom",
)


def _resolve(config_path: Path, value: str) -> Path:
    candidate = Path(value).expanduser()
    return candidate if candidate.is_absolute() else (config_path.parent / candidate).resolve()


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ConfigError(f"{path} contains no rows")
    return rows


def _pressure_by_temperature(path: Path) -> dict[float, float]:
    result: dict[float, float] = {}
    for row in _read_rows(path):
        try:
            temperature, pressure = float(row["T_K"]), float(row["P_GPa"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ConfigError(f"{path} must contain numeric T_K and P_GPa") from exc
        previous = result.setdefault(temperature, pressure)
        if not math.isclose(previous, pressure, rel_tol=0.0, abs_tol=1.0e-12):
            raise ConfigError(f"{path} has multiple pressures at T={temperature:g}; QH source is ambiguous")
    return result


def _volume_per_atom(extxyz: Path) -> float:
    lines = extxyz.read_text(encoding="utf-8").splitlines()
    if len(lines) < 2:
        raise ConfigError(f"{extxyz} is not an extxyz structure")
    try:
        natoms = int(lines[0].strip())
    except ValueError as exc:
        raise ConfigError(f"{extxyz} has an invalid atom count") from exc
    lattice = re.search(r'Lattice="([^"]+)"', lines[1])
    if lattice is None:
        raise ConfigError(f"{extxyz} lacks an extxyz Lattice field")
    try:
        values = [float(value) for value in lattice.group(1).split()]
    except ValueError as exc:
        raise ConfigError(f"{extxyz} has a non-numeric lattice") from exc
    if len(values) != 9 or natoms <= 0:
        raise ConfigError(f"{extxyz} must have 9 lattice components and positive atom count")
    a, b, c, d, e, f, g, h, i = values
    determinant = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
    volume = abs(determinant) / natoms
    if volume <= 0.0:
        raise ConfigError(f"{extxyz} has zero cell volume")
    return volume


def _source_phase_summary(summary: dict[str, Any], phase: str) -> dict[str, Any]:
    phases = summary.get("phases")
    if not isinstance(phases, dict) or not isinstance(phases.get(phase), dict):
        raise ConfigError(f"source summary does not describe phase {phase!r}")
    return dict(phases[phase])


def materialize(config_path: str | Path) -> dict[str, Any]:
    """Convert configured diagnostic files and return a compact manifest."""
    config_file = Path(config_path).expanduser().resolve()
    config = load_mapping(config_file)
    for key in ("data_root", "source_root", "system"):
        if not isinstance(config.get(key), str):
            raise ConfigError(f"{key} is required")
    phases = config.get("phases")
    reliability = config.get("qh_reliable")
    if not isinstance(phases, list) or not phases or not all(isinstance(phase, str) for phase in phases):
        raise ConfigError("phases must be a non-empty string list")
    if not isinstance(reliability, dict) or set(reliability) != set(phases) or not all(
        isinstance(reliability[phase], bool) for phase in phases
    ):
        raise ConfigError("qh_reliable must give an explicit boolean for every phase")

    data_system = _resolve(config_file, config["data_root"]) / config["system"]
    source_root = _resolve(config_file, config["source_root"])
    source_summary_path = source_root / str(config.get("source_summary", "qh_summary.json"))
    try:
        summary = json.loads(source_summary_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ConfigError(f"cannot read QH source summary {source_summary_path}: {exc}") from exc
    if not isinstance(summary, dict):
        raise ConfigError("QH source summary must be a JSON mapping")

    diagnostic_only = bool(config.get("diagnostic_only", True))
    source_pattern = str(config.get("source_pattern", "{phase}_fqh.csv"))
    manifest: dict[str, Any] = {
        "system": config["system"],
        "source_root": str(source_root),
        "diagnostic_only": diagnostic_only,
        "phases": {},
    }
    for phase in phases:
        phase_dir = data_system / phase
        reference_pressure = _pressure_by_temperature(phase_dir / "reference_G.csv")
        base_volume = _volume_per_atom(phase_dir / "structure.extxyz")
        source_path = source_root / source_pattern.format(phase=phase)
        source_rows = _read_rows(source_path)
        normalized: list[dict[str, str]] = []
        for row in source_rows:
            try:
                temperature = float(row["T_K"])
                fqh = float(row["F_QH_eV_per_atom"])
                static = float(row["E_static_eV_per_atom"])
                scale = float(row["volume_scale"])
            except (KeyError, TypeError, ValueError) as exc:
                raise ConfigError(f"{source_path} lacks a numeric prototype-QH field") from exc
            if temperature not in reference_pressure:
                raise ConfigError(f"{source_path} has T={temperature:g} not present in reference_G.csv")
            if scale <= 0.0 or not all(math.isfinite(value) for value in (fqh, static, scale)):
                raise ConfigError(f"{source_path} has an invalid QH row at T={temperature:g}")
            normalized.append(
                {
                    "T_K": f"{temperature:.12g}",
                    "P_GPa": f"{reference_pressure[temperature]:.12g}",
                    "V_min_A3_per_atom": f"{base_volume * scale:.16g}",
                    "F_vib_eV_per_atom": f"{fqh - static:.16g}",
                    "F_QH_eV_per_atom": f"{fqh:.16g}",
                }
            )
        target = phase_dir / "fqh.csv"
        with target.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=FQH_FIELDS)
            writer.writeheader()
            writer.writerows(normalized)
        phase_summary = _source_phase_summary(summary, phase)
        report = {
            "qh_reliable": reliability[phase],
            "diagnostic_only": diagnostic_only,
            "source": {
                "prototype_csv": f"{config['source_root'].rstrip('/')}/{source_pattern.format(phase=phase)}",
                "summary": f"{config['source_root'].rstrip('/')}/{config.get('source_summary', 'qh_summary.json')}",
                "checkpoint": summary.get("checkpoint"),
                "head": summary.get("head"),
                "supercell_matrix": summary.get("supercell_matrix"),
                "mesh": summary.get("mesh"),
                "volume_scales": summary.get("volume_scales"),
            },
            "n_temperatures": len(normalized),
            "reference_pressure_inherited": True,
            "base_volume_A3_per_atom": base_volume,
            "minimum_frequency_THz_by_volume": phase_summary.get("minimum_frequency_THz_by_volume"),
            "imaginary_mode_fraction": None,
            "limitations": [
                "Materialized from a three-volume diagnostic prototype; no Vinet fit is claimed.",
                "Imaginary-mode fraction and force-constant cache were not produced by the source run.",
            ],
        }
        if not reliability[phase]:
            report["limitations"].append("This curve failed the Hf QH reliability/accuracy gate and is excluded from production predictors.")
        else:
            report["limitations"].append("Reliability is unlocked only for this hcp head after the recorded imaginary-mode and smooth-deviation gate; absolute reference agreement remains a separate limitation.")
        (phase_dir / "phonon_report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        manifest["phases"][phase] = {"fqh_csv": str(target), "phonon_report": str(phase_dir / "phonon_report.json"), "n_temperatures": len(normalized)}
        print(f"QH diagnostic materialized: {config['system']}/{phase} rows={len(normalized)} reliable={reliability[phase]}")
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    args = parser.parse_args(argv)
    try:
        materialize(args.config)
    except (ConfigError, OSError, ValueError, KeyError) as exc:
        print(f"QH diagnostic materialization failed: {exc}")
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
