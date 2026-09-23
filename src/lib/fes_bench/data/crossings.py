"""Compute reference free-energy crossing topology on a common system grid."""

from __future__ import annotations

import argparse
import csv
import itertools
import json
from pathlib import Path

from fes_bench.config import ConfigError, load_mapping


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="YAML/JSON crossing configuration")
    return parser


def _resolve(config_path: Path, value: str) -> Path:
    path = Path(value).expanduser()
    return path if path.is_absolute() else (config_path.parent / path).resolve()


def _table(path: Path) -> dict[float, float]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
    pressures = {float(row["P_GPa"]) for row in rows}
    if len(pressures) != 1:
        raise ConfigError(f"{path} has multiple pressures; pressure-aware crossing grids are not yet implemented")
    return {float(row["T_K"]): float(row["G_eV_per_atom"]) for row in rows}


def _roots(temperatures: list[float], delta: list[float]) -> list[dict[str, float]]:
    roots: list[dict[str, float]] = []
    index = 0
    while index < len(temperatures) - 1:
        t0, t1 = temperatures[index], temperatures[index + 1]
        d0, d1 = delta[index], delta[index + 1]
        if d0 == 0.0:
            start = index
            while index + 1 < len(delta) and delta[index + 1] == 0.0:
                index += 1
            end = index
            slopes: list[float] = []
            if start > 0:
                slopes.append(abs(delta[start] - delta[start - 1]) / (temperatures[start] - temperatures[start - 1]) * 1000.0)
            if end + 1 < len(delta):
                slopes.append(abs(delta[end + 1] - delta[end]) / (temperatures[end + 1] - temperatures[end]) * 1000.0)
            roots.append({"T_K": (temperatures[start] + temperatures[end]) / 2.0, "slope_meV_per_atom_per_K": max(slopes, default=0.0), "quantized_plateau": end > start})
        elif d0 * d1 < 0.0:
            slope = (d1 - d0) / (t1 - t0)
            roots.append({"T_K": t0 - d0 / slope, "slope_meV_per_atom_per_K": abs(slope) * 1000.0})
        index += 1
    return roots
def run(config_path: str | Path) -> int:
    config_file = Path(config_path).expanduser().resolve()
    config = load_mapping(config_file)
    if not isinstance(config.get("data_root"), str) or not isinstance(config.get("system"), str):
        raise ConfigError("data_root and system are required strings")
    system_dir = _resolve(config_file, config["data_root"]) / config["system"]
    system_meta = json.loads((system_dir / "system.json").read_text(encoding="utf-8"))
    phases = system_meta["phases"]
    tables = {phase: _table(system_dir / phase / "reference_G.csv") for phase in phases}
    temperatures = sorted(set.intersection(*(set(table) for table in tables.values())))
    if len(temperatures) < 2:
        raise ConfigError(f"{system_dir} has fewer than two common temperatures")
    pairs: dict[str, dict[str, object]] = {}
    for left, right in itertools.combinations(phases, 2):
        delta = [tables[left][temperature] - tables[right][temperature] for temperature in temperatures]
        pairs[f"{left}_minus_{right}"] = {
            "left": left,
            "right": right,
            "crossings": _roots(temperatures, delta),
            "min_delta_eV_per_atom": min(delta),
            "max_delta_eV_per_atom": max(delta),
        }
    output = {
        "system": config["system"],
        "pressure_GPa": float(system_meta.get("reference_grid", {}).get("P_GPa", [0.0])[0]),
        "grid": {"n_points": len(temperatures), "T_min_K": temperatures[0], "T_max_K": temperatures[-1], "step_K": 1.0},
        "pairs": pairs,
    }
    output_path = system_dir / "reference_crossings.json"
    output_path.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    for name, values in pairs.items():
        roots = values["crossings"]
        print(f"{name}: crossings={len(roots)} {roots}")
    return 0


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        return run(args.config)
    except (ConfigError, OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
        print(f"crossing computation failed: {exc}")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
