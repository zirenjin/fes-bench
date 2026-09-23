"""Per-phase structure and QH inventory — plan Table 2."""

from __future__ import annotations

import argparse
from pathlib import Path

from common import csv_write, inventory, meta_write


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--output-root", type=Path, default=Path("result/tables"))
    args = parser.parse_args()
    root = args.repo_root.resolve()
    output = args.output_root if args.output_root.is_absolute() else root / args.output_root
    inv = inventory(root)
    rows, missing = [], []
    inputs = [root / "result/experiments/data_prep/raw_inventory/inventory/seed_none/metrics.json"]
    qh_metrics_paths = sorted(root.glob("result/experiments/quasi_harmonic/raw_runs/*/seed_none/metrics.json"))
    for qh_metrics_path in qh_metrics_paths:
        inputs.append(qh_metrics_path)
        import json
        qh_metrics = json.loads(qh_metrics_path.read_text(encoding="utf-8"))
        missing.extend(qh_metrics.get("provenance", {}).get("missing", []))
    for key, item in sorted(inv["phases"].items()):
        meta = item.get("meta", {})
        representative = meta.get("representative", {})
        sg = representative.get("space_group", {})
        edges = item.get("lattice_edge_lengths_A", [])
        qh = item.get("qh", {})
        rows.append({"system": item.get("system"), "phase": item.get("phase"), "space_group_before": sg.get("before", representative.get("space_group_before", "")), "space_group_after": sg.get("after", representative.get("space_group_after", "")), "primitive_atoms": representative.get("primitive_atoms", item.get("atom_count", "")), "representative_source": representative.get("source", representative.get("source_structure", "")), "relaxation_status": representative.get("status", ""), "qh_supercell": qh.get("supercell", ""), "shortest_edge_A": item.get("shortest_edge_A", ""), "q_mesh": qh.get("mesh", ""), "imaginary_fraction": qh.get("imaginary_fraction", ""), "qh_reliable": qh.get("qh_reliable", "")})
        if not qh:
            missing.append({"field": f"{key}.qh", "reason": "QH raw summary is not available in the repository"})
        if not edges:
            missing.append({"field": f"{key}.shortest_edge_A", "reason": "structure lattice is not present in normalized inventory"})
    fields = ["system", "phase", "space_group_before", "space_group_after", "primitive_atoms", "representative_source", "relaxation_status", "qh_supercell", "shortest_edge_A", "q_mesh", "imaginary_fraction", "qh_reliable"]
    csv_path = output / "phase_inventory.csv"
    csv_write(csv_path, fields, rows)
    meta_write(root, csv_path, inputs, missing)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
