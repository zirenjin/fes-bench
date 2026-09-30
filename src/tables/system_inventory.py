"""System inventory and reference-scale summary — plan Table 1.

Reads the normalized inventory and reference-statistics raw runs only.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from common import csv_write, inventory, meta_write


def _head_policy(root: Path) -> dict[str, str]:
    """Read system energy heads from the checked-in head policy."""
    policy = root / "configs/models/head_policy.yaml"
    try:
        payload = json.loads(policy.read_text(encoding="utf-8"))
        return {str(system): str(values.get("energy_head", "not reported")) for system, values in payload.get("systems", {}).items()}
    except (OSError, json.JSONDecodeError):
        pass
    result: dict[str, str] = {}
    current = None
    for line in policy.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if stripped.endswith(": {"):
            current = stripped.split(":", 1)[0].strip('"')
        elif current and stripped.startswith("energy_head:"):
            result[current] = stripped.split(":", 1)[1].strip().strip('"')
        elif stripped == "}":
            current = None
    return result


def _frozen_splits(root: Path, system: str) -> str:
    names: list[str] = []
    for split_root in (root / "data/processed/splits", root / "data/processed/splits_v2"):
      version = "v2" if split_root.name.endswith("v2") else "v1"
      for path in sorted(split_root.glob("*.json")):
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        rows = payload.get("train", []) + payload.get("test", [])
        if not rows and isinstance(payload.get("folds"), dict):
            rows = [row for fold in payload["folds"].values() for subset in ("train", "test") for row in fold.get(subset, [])]
        if any(str(row.get("system")) == system for row in rows):
            names.append(f"{version}:{path.stem}")
    return ";".join(names) or "not in frozen v1 splits"


def _reference_contract(phase_records: list[dict]) -> tuple[str, str, str]:
    values = {key: set() for key in ("supercell", "convergence", "size_error_eV_per_atom")}
    for item in phase_records:
        meta = item.get("meta", {})
        for key in values:
            value = meta.get(key, "not reported")
            # A numeric zero in the source schema is explicitly documented as
            # a placeholder, not as a measured finite-size error.
            if key == "size_error_eV_per_atom" and value == 0:
                value = "not reported"
            values[key].add(str(value))
    def render(key: str) -> str:
        return ";".join(sorted(values[key])) if values[key] else "not reported"
    return render("supercell"), render("convergence"), render("size_error_eV_per_atom")


def display_functional(phase_records: list[dict]) -> str:
    """Render the functional from processed metadata, including dispersion."""

    values: set[str] = set()
    for item in phase_records:
        meta = item.get("meta", {})
        functional = str(meta.get("functional", "")).strip()
        dispersion = str(meta.get("dispersion", "")).strip()
        if functional and dispersion and dispersion.lower() not in {"none", "n/a", "na"}:
            values.add(f"{functional}-{dispersion}")
        elif functional:
            values.add(functional)
    return ";".join(sorted(values))


def qh_counts_from_phase_inventory(path: Path) -> dict[str, int]:
    """Derive per-system reliable-QH counts from Table 2, never independently."""

    counts: dict[str, int] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if str(row.get("qh_reliable", "")).strip().lower() == "true":
                system = str(row["system"])
                counts[system] = counts.get(system, 0) + 1
    return counts


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--output-root", type=Path, default=Path("result/tables"))
    args = parser.parse_args()
    root = args.repo_root.resolve()
    output = args.output_root if args.output_root.is_absolute() else root / args.output_root
    inv = inventory(root)
    amplitudes: dict[tuple[str, str], dict[str, float]] = {}
    inputs: list[Path] = [root / "result/experiments/data_prep/raw_inventory/inventory/seed_none/metrics.json"]
    inputs.extend(sorted((root / "data/processed/splits").glob("*.json")))
    inputs.extend(sorted((root / "data/processed/splits_v2").glob("*.json")))
    statistics_path = root / "result/experiments/reference_statistics/raw_runs/reference_statistics/seed_none/metrics.json"
    if statistics_path.exists():
        inputs.append(statistics_path)
        body = json.loads(statistics_path.read_text(encoding="utf-8"))
        for record in body["metrics"].get("records", []):
            key = (str(record.get("system")), str(record.get("pair")))
            amplitudes[key] = record
    phase_inventory = output / "phase_inventory.csv"
    if not phase_inventory.exists():
        raise FileNotFoundError("Generate phase_inventory.csv before system_inventory.csv so QH counts have one source of truth.")
    qh_counts = qh_counts_from_phase_inventory(phase_inventory)
    inputs.append(phase_inventory)
    rows = []
    missing = []
    heads = _head_policy(root)
    for system, data in sorted(inv["systems"].items()):
        phase_names = list(data.get("phases", []))
        phase_records = [inv["phases"].get(f"{system}:{phase}", {}) for phase in phase_names]
        reference_supercell, reference_convergence, reference_finite_size_error = _reference_contract(phase_records)
        points = [item.get("reference_points") for item in phase_records if item.get("reference_points") is not None]
        temperatures = data.get("reference_grid", {}).get("T_K", [])
        stds = [float(item["std_eV_per_atom"]) * 1000.0 for (item_system, _), item in amplitudes.items() if item_system == system]
        rows.append({"system": system, "phases": ";".join(phase_names), "n_phases": len(phase_names), "type_map": ";".join(data.get("type_map", [])), "truth_level": data.get("truth_level", ""), "functional": display_functional(phase_records), "T_min_K": min(temperatures) if temperatures else "", "T_max_K": max(temperatures) if temperatures else "", "grid_points": min(points) if points and len(set(points)) == 1 else ";".join(str(point) for point in points), "delta_G_std_min_meV_per_atom": min(stds) if stds else "", "delta_G_std_max_meV_per_atom": max(stds) if stds else "", "relative_numbers": len(phase_names) * (len(phase_names) - 1) // 2, "qh_reliable_count": qh_counts.get(system, 0), "doi": ";".join(sorted({str(item.get("meta", {}).get("doi", "")) for item in phase_records})), "reference_supercell": reference_supercell, "reference_convergence": reference_convergence, "reference_finite_size_error": reference_finite_size_error, "in_frozen_splits": _frozen_splits(root, system), "head": heads.get(system, "not reported")})
        if not stds:
            missing.append({"field": f"{system}.delta_G_std", "reason": "no reference-statistics raw run for this system"})
    fields = ["system", "phases", "n_phases", "type_map", "truth_level", "functional", "T_min_K", "T_max_K", "grid_points", "delta_G_std_min_meV_per_atom", "delta_G_std_max_meV_per_atom", "relative_numbers", "qh_reliable_count", "doi", "reference_supercell", "reference_convergence", "reference_finite_size_error", "in_frozen_splits", "head"]
    csv_path = output / "system_inventory.csv"
    csv_write(csv_path, fields, rows)
    meta_write(root, csv_path, inputs, missing)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
