"""Freeze reproducible temperature, phase, and system holdout splits."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Any

from fes_bench.config import ConfigError, load_mapping
from fes_bench.data.load import load


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    return parser


def _resolve(config_path: Path, value: str, default: str) -> Path:
    raw = Path(value if value else default).expanduser()
    return raw.resolve() if raw.is_absolute() else (config_path.parent / raw).resolve()


def _triples(
    data_root: Path, systems: list[str]
) -> tuple[
    dict[str, list[dict[str, Any]]],
    dict[str, list[str]],
    dict[tuple[str, str, int], float],
]:
    rows: dict[str, list[dict[str, Any]]] = defaultdict(list)
    phases_by_system: dict[str, list[str]] = {}
    temperatures: dict[tuple[str, str, int], float] = {}
    for system in systems:
        system_meta = json.loads((data_root / system / "system.json").read_text(encoding="utf-8"))
        phases_by_system[system] = list(system_meta["phases"])
        for phase in phases_by_system[system]:
            table = load(system, phase, data_root).G_table
            for index, point in enumerate(table):
                rows[system].append({"system": system, "phase": phase, "T_index": index})
                temperatures[(system, phase, index)] = point.T_K
    return dict(rows), phases_by_system, temperatures


def _temperature_threshold(
    data_root: Path,
    system: str,
    rows: list[dict[str, Any]],
    temperatures: dict[tuple[str, str, int], float],
) -> float:
    crossings = json.loads((data_root / system / "reference_crossings.json").read_text(encoding="utf-8"))
    roots = [item["T_K"] for pair in crossings["pairs"].values() for item in pair.get("crossings", [])]
    if roots:
        return float(min(roots) - 100.0)
    values = sorted(temperatures[(system, row["phase"], row["T_index"])] for row in rows)
    return float(values[int(0.7 * (len(values) - 1))])


def _git_commit(root: Path) -> str:
    try:
        return subprocess.check_output(
            ["git", "-C", str(root), "rev-parse", "HEAD"], text=True, stderr=subprocess.DEVNULL
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unavailable"


def _write(path: Path, payload: dict[str, Any]) -> str:
    canonical = json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
    digest = hashlib.sha256(canonical.encode("utf-8")).hexdigest()
    payload["sha256"] = digest
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return digest


def run(config_path: str | Path) -> int:
    config_file = Path(config_path).expanduser().resolve()
    config = load_mapping(config_file)
    data_root_setting = str(config.get("data_root", "../data"))
    data_root = _resolve(config_file, data_root_setting, "../data")
    output_root = _resolve(config_file, str(config.get("output_root", "../splits")), "../splits")
    systems = config.get("systems")
    if not isinstance(systems, list) or not systems or not all(isinstance(item, str) for item in systems):
        raise ConfigError("systems must be a non-empty list of canonical system names")
    rows, phases_by_system, temperatures = _triples(data_root, systems)
    common = {
        "schema_version": 1,
        "data_root": data_root_setting,
        "systems": systems,
        "generation_parameters": {"crossing_buffer_K": 100.0},
        "git_commit": _git_commit(config_file.parents[1]),
    }

    thresholds = {
        system: _temperature_threshold(data_root, system, entries, temperatures)
        for system, entries in rows.items()
    }
    train: list[dict[str, Any]] = []
    test: list[dict[str, Any]] = []
    for system, entries in rows.items():
        for row in entries:
            temperature = temperatures[(system, row["phase"], row["T_index"])]
            (train if temperature <= thresholds[system] else test).append(row)
    temp = {**common, "name": "temp_extrap", "threshold_T_K": thresholds, "train": train, "test": test}
    print(f"temp_extrap: train={len(train)} test={len(test)} sha256={_write(output_root / 'temp_extrap.json', temp)}")

    lopo_folds: dict[str, dict[str, list[dict[str, Any]]]] = {}
    for system, phases in phases_by_system.items():
        for held_phase in phases:
            lopo_folds[f"{system}:{held_phase}"] = {
                "train": [row for entries in rows.values() for row in entries if not (row["system"] == system and row["phase"] == held_phase)],
                "test": [row for row in rows[system] if row["phase"] == held_phase],
            }
    lopo = {**common, "name": "phase_lopo", "folds": lopo_folds}
    print(f"phase_lopo: folds={len(lopo_folds)} sha256={_write(output_root / 'phase_lopo.json', lopo)}")

    loso_folds: dict[str, dict[str, list[dict[str, Any]]]] = {}
    for held_system in systems:
        loso_folds[held_system] = {
            "train": [row for system, entries in rows.items() if system != held_system for row in entries],
            "test": rows[held_system],
        }
    loso = {**common, "name": "system_loso", "folds": loso_folds}
    print(f"system_loso: folds={len(loso_folds)} sha256={_write(output_root / 'system_loso.json', loso)}")
    return 0


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        return run(args.config)
    except (ConfigError, OSError, KeyError, ValueError) as exc:
        print(f"split generation failed: {exc}")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
