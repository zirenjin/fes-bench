"""System inventory and reference-scale summary — plan Table 1.

Reads the normalized inventory and reference-statistics raw runs only.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from common import csv_write, inventory, meta_write, raw_runs


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--output-root", type=Path, default=Path("results/tables"))
    args = parser.parse_args()
    root = args.repo_root.resolve()
    output = args.output_root if args.output_root.is_absolute() else root / args.output_root
    inv = inventory(root)
    amplitudes: dict[tuple[str, str], dict[str, float]] = {}
    inputs: list[Path] = [root / "results/raw_runs/_inventory/inventory/seed_none/metrics.json"]
    for path, body in raw_runs(root, "reference_statistics"):
        inputs.append(path)
        for record in body["metrics"].get("records", []):
            key = (str(record.get("system")), str(record.get("pair")))
            amplitudes[key] = record
    rows = []
    missing = []
    for system, data in sorted(inv["systems"].items()):
        phase_names = list(data.get("phases", []))
        phase_records = [inv["phases"].get(f"{system}:{phase}", {}) for phase in phase_names]
        points = [item.get("reference_points") for item in phase_records if item.get("reference_points") is not None]
        temperatures = data.get("reference_grid", {}).get("T_K", [])
        stds = [float(item["std_eV_per_atom"]) * 1000.0 for (item_system, _), item in amplitudes.items() if item_system == system]
        rows.append({"system": system, "phases": ";".join(phase_names), "n_phases": len(phase_names), "type_map": ";".join(data.get("type_map", [])), "truth_level": data.get("truth_level", ""), "functional": ";".join(sorted({str(item.get("meta", {}).get("functional", "")) for item in phase_records})), "T_min_K": min(temperatures) if temperatures else "", "T_max_K": max(temperatures) if temperatures else "", "grid_points": min(points) if points and len(set(points)) == 1 else ";".join(str(point) for point in points), "delta_G_std_min_meV_per_atom": min(stds) if stds else "", "delta_G_std_max_meV_per_atom": max(stds) if stds else "", "relative_numbers": len(phase_names) * (len(phase_names) - 1) // 2, "qh_reliable_count": sum(item.get("qh_reliable") is True for item in phase_records), "doi": ";".join(sorted({str(item.get("meta", {}).get("doi", "")) for item in phase_records}))})
        if not stds:
            missing.append({"field": f"{system}.delta_G_std", "reason": "no reference-statistics raw run for this system"})
    fields = ["system", "phases", "n_phases", "type_map", "truth_level", "functional", "T_min_K", "T_max_K", "grid_points", "delta_G_std_min_meV_per_atom", "delta_G_std_max_meV_per_atom", "relative_numbers", "qh_reliable_count", "doi"]
    csv_path = output / "system_inventory.csv"
    csv_write(csv_path, fields, rows)
    meta_write(root, csv_path, inputs, missing)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
