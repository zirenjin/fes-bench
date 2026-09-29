"""Per-phase structure and QH inventory — plan Table 2."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from common import csv_write, inventory, meta_write


def qh_from_report(root: Path, system: str, phase: str, fallback: dict) -> tuple[dict, Path | None]:
    """Prefer the per-phase phonon report over a hand-entered inventory field."""

    report = root / "data/processed" / system / phase / "phonon_report.json"
    if not report.exists():
        return fallback, None
    payload = json.loads(report.read_text(encoding="utf-8"))
    source = payload.get("source", {}) if isinstance(payload.get("source"), dict) else {}
    return {
        "supercell": payload.get("supercell_matrix", source.get("supercell_matrix", "")),
        "mesh": payload.get("mesh", source.get("mesh", "")),
        "imaginary_fraction": payload.get("imaginary_mode_fraction", ""),
        "qh_reliable": payload.get("qh_reliable", ""),
    }, report


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
    for key, item in sorted(inv["phases"].items()):
        meta = item.get("meta", {})
        representative = meta.get("representative", {})
        sg = representative.get("space_group", {})
        edges = item.get("lattice_edge_lengths_A", [])
        qh, report = qh_from_report(root, str(item.get("system")), str(item.get("phase")), item.get("qh", {}))
        if report is not None:
            inputs.append(report)
        source = f"data/processed/{item.get('system')}/{item.get('phase')}/structure.extxyz" if representative.get("status") == "relaxed" else representative.get("source", representative.get("source_structure", ""))
        rows.append({"system": item.get("system"), "phase": item.get("phase"), "space_group_before": sg.get("before", representative.get("space_group_before", "")), "space_group_after": sg.get("after", representative.get("space_group_after", "")), "primitive_atoms": representative.get("primitive_atoms", item.get("atom_count", "")), "representative_source": source, "relaxation_status": representative.get("status", ""), "qh_supercell": qh.get("supercell", ""), "shortest_edge_A": item.get("shortest_edge_A", ""), "q_mesh": qh.get("mesh", ""), "imaginary_fraction": qh.get("imaginary_fraction", ""), "qh_reliable": qh.get("qh_reliable", "")})
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
