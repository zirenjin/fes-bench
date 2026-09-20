"""Normalize per-phase thermodynamic tables into the benchmark reference schema."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from fes_bench.config import ConfigError, load_mapping


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="YAML/JSON import configuration")
    return parser


def _resolve(config_path: Path, value: str) -> Path:
    path = Path(value).expanduser()
    return path if path.is_absolute() else (config_path.parent / path).resolve()


def _phase_rows(source: Path, *, delimiter: str, temperature_column: str, energy_column: str, pressure_gpa: float, energy_factor: float, level: str, provenance: str, skip_rows: int = 0) -> list[dict[str, str]]:
    with source.open("r", encoding="utf-8", newline="") as handle:
        for _ in range(skip_rows):
            next(handle, None)
        reader = csv.DictReader(handle, delimiter=delimiter)
        if reader.fieldnames is None or temperature_column not in reader.fieldnames or energy_column not in reader.fieldnames:
            raise ConfigError(f"{source} lacks configured columns {temperature_column!r}, {energy_column!r}")
        rows: list[dict[str, str]] = []
        for line, row in enumerate(reader, start=2):
            try:
                temperature = float(row[temperature_column])
                energy = float(row[energy_column]) * energy_factor
            except (TypeError, ValueError) as exc:
                raise ConfigError(f"{source}:{line} has non-numeric thermodynamic data") from exc
            rows.append(
                {
                    "T_K": f"{temperature:.12g}",
                    "P_GPa": f"{pressure_gpa:.12g}",
                    "G_eV_per_atom": f"{energy:.16g}",
                    "level": level,
                    "source": provenance,
                }
            )
    if not rows:
        raise ConfigError(f"{source} contains no thermodynamic rows")
    return rows


def run(config_path: str | Path) -> int:
    config_file = Path(config_path).expanduser().resolve()
    config = load_mapping(config_file)
    raw_root = config.get("data_root")
    system = config.get("system")
    phase_specs = config.get("phases")
    if not isinstance(raw_root, str) or not isinstance(system, str) or not isinstance(phase_specs, list) or not phase_specs:
        raise ConfigError("data_root, system, and a non-empty phases list are required")
    system_dir = _resolve(config_file, raw_root) / system
    system_dir.mkdir(parents=True, exist_ok=True)
    shared = config.get("shared_meta")
    if not isinstance(shared, dict):
        raise ConfigError("shared_meta must be a mapping")
    all_temperatures: set[float] = set()
    phase_names: list[str] = []
    for spec in phase_specs:
        if not isinstance(spec, dict):
            raise ConfigError("each phase specification must be a mapping")
        name = spec.get("name")
        source_value = spec.get("source")
        if not isinstance(name, str) or not isinstance(source_value, str):
            raise ConfigError("each phase requires name and source")
        rows = _phase_rows(
            _resolve(config_file, source_value),
            delimiter=str(spec.get("delimiter", ",")),
            temperature_column=str(spec.get("temperature_column", "temperature")),
            energy_column=str(spec.get("energy_column", "Gibbs_energy")),
            pressure_gpa=float(spec.get("pressure_gpa", 0.0)),
            energy_factor=float(spec.get("energy_to_eV_per_atom", 1.0)),
            level=str(spec.get("level", "reference")),
            provenance=str(spec.get("provenance", source_value)),
            skip_rows=int(spec.get("skip_rows", 0)),
        )
        phase_dir = system_dir / name
        phase_dir.mkdir(parents=True, exist_ok=True)
        with (phase_dir / "reference_G.csv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["T_K", "P_GPa", "G_eV_per_atom", "level", "source"])
            writer.writeheader()
            writer.writerows(rows)
        meta = dict(shared)
        meta.update({"interpolated": False, "raw_table": source_value, "unit_conversion": {"energy_to_eV_per_atom": float(spec.get("energy_to_eV_per_atom", 1.0)), "pressure_to_GPa": 1.0}})
        (phase_dir / "meta.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
        phase_names.append(name)
        all_temperatures.update(float(row["T_K"]) for row in rows)
        print(f"normalized: {system}/{name} rows={len(rows)}")
    system_meta = {
        "phases": phase_names,
        "type_map": config.get("type_map"),
        "reference_grid": {"T_K": sorted(all_temperatures), "P_GPa": config.get("pressures_GPa", [0.0])},
        "truth_level": config.get("truth_level"),
    }
    (system_dir / "system.json").write_text(json.dumps(system_meta, indent=2) + "\n", encoding="utf-8")
    print(f"system metadata: {system_dir / 'system.json'}")
    return 0


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        return run(args.config)
    except (ConfigError, OSError, ValueError) as exc:
        print(f"thermodynamic import failed: {exc}")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
