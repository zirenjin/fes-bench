"""Freeze a common exact T/P grid across all phases in one system."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from fes_bench.config import ConfigError, load_mapping


FIELDS = ("T_K", "P_GPa", "G_eV_per_atom", "level", "source")


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="YAML/JSON grid-alignment configuration")
    return parser


def _resolve(config_path: Path, value: str) -> Path:
    path = Path(value).expanduser()
    return path if path.is_absolute() else (config_path.parent / path).resolve()


def _read(path: Path) -> dict[tuple[float, float], dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != FIELDS:
            raise ConfigError(f"noncanonical reference table: {path}")
        values: dict[tuple[float, float], dict[str, str]] = {}
        for row in reader:
            key = (float(row["T_K"]), float(row["P_GPa"]))
            if key in values:
                raise ConfigError(f"duplicate T/P row in {path}: {key}")
            values[key] = row
    return values


def run(config_path: str | Path) -> int:
    config_file = Path(config_path).expanduser().resolve()
    config = load_mapping(config_file)
    if not isinstance(config.get("data_root"), str) or not isinstance(config.get("system"), str):
        raise ConfigError("data_root and system are required strings")
    root = _resolve(config_file, config["data_root"])
    system_dir = root / config["system"]
    system_meta_path = system_dir / "system.json"
    system_meta = json.loads(system_meta_path.read_text(encoding="utf-8"))
    phases = system_meta.get("phases")
    if not isinstance(phases, list) or not phases:
        raise ConfigError(f"invalid phase list in {system_meta_path}")
    tables = {phase: _read(system_dir / phase / "reference_G.csv") for phase in phases}
    common = set.intersection(*(set(table) for table in tables.values()))
    if not common:
        raise ConfigError(f"no shared T/P grid for {config['system']}")
    ordered = sorted(common)
    for phase, table in tables.items():
        original_count = len(table)
        path = system_dir / phase / "reference_G.csv"
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS)
            writer.writeheader()
            writer.writerows(table[key] for key in ordered)
        meta_path = system_dir / phase / "meta.json"
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        meta["interpolated"] = original_count != len(ordered)
        meta["grid_harmonization"] = {
            "original_points": original_count,
            "common_points": len(ordered),
            "method": "exact T/P intersection; no free-energy value interpolation",
        }
        meta_path.write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
        print(f"aligned: {config['system']}/{phase} {original_count} -> {len(ordered)}")
    system_meta["reference_grid"] = {
        "T_K": [key[0] for key in ordered],
        "P_GPa": sorted({key[1] for key in ordered}),
    }
    system_meta_path.write_text(json.dumps(system_meta, indent=2) + "\n", encoding="utf-8")
    print(f"common grid: points={len(ordered)} T_K=[{ordered[0][0]:g}, {ordered[-1][0]:g}]")
    return 0


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        return run(args.config)
    except (ConfigError, OSError, ValueError, json.JSONDecodeError) as exc:
        print(f"grid alignment failed: {exc}")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
