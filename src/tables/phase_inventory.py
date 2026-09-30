"""Per-phase structure and QH inventory — plan Table 2b.

The table is a read-only view of processed structures and QH diagnostics. It
never accepts hand-entered reliability flags: ``qh_reliable`` is derived from
the recorded imaginary-mode fraction and the configured 0.01 threshold.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path
from typing import Any

from common import csv_write, inventory, meta_write

IMAGINARY_THRESHOLD = 0.01


def _representative_origin(representative: dict[str, Any]) -> str:
    """Render AFLOW/MP/DaRUS provenance for the representative structure."""
    parts: list[str] = []
    if representative.get("prototype"):
        parts.append(f"AFLOW:{representative['prototype']}")
    if representative.get("material_id"):
        parts.append(f"MP:{representative['material_id']}")
    if representative.get("icsd"):
        parts.append(f"ICSD:{representative['icsd']}")
    if representative.get("darus_phase"):
        parts.append(f"DaRUS:{representative['darus_phase']}")
    source = representative.get("source") or representative.get("source_structure")
    if source and not parts:
        parts.append(f"file:{source}")
    return ";".join(str(item) for item in parts) or "not reported"


def _lattice(path: Path) -> tuple[list[float], int | None]:
    if not path.exists():
        return [], None
    lines = path.read_text(encoding="utf-8").splitlines()[:2]
    atoms = None
    try:
        atoms = int(lines[0].strip())
    except (IndexError, ValueError):
        pass
    match = re.search(r'Lattice="([^"]+)"', " ".join(lines))
    if not match:
        return [], atoms
    values = [float(value) for value in match.group(1).split()]
    vectors = [values[row * 3 : row * 3 + 3] for row in range(3)]
    return [math.sqrt(sum(component * component for component in vector)) for vector in vectors], atoms


def _candidate_qh(root: Path, system: str, phase: str) -> tuple[dict[str, Any], Path | None, dict[str, Any]]:
    # One canonical diagnostic source per domain.  SiO2 uses the adopted
    # Domains_SSE_PBE run at equilibrium volume and [2,2,2]; the superseded
    # Domains_Alloy run is archived and must not silently win by glob ordering.
    if system == "sio2":
        candidates = sorted(root.glob("result/experiments/quasi_harmonic_sio2_sse_pbe/raw_runs/*/qh_summary.json"))
    else:
        candidates = sorted(root.glob("result/experiments/quasi_harmonic_10a/raw_runs/*/qh_summary.json"))
    for path in candidates:
        try:
            summary = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        if str(summary.get("system")) != system:
            continue
        phase_record = summary.get("phases", {}).get(phase)
        if isinstance(phase_record, dict):
            return phase_record, path, summary
    return {}, None, {}


def _diagnostic(root: Path, system: str, phase: str, summary_path: Path | None, phase_record: dict[str, Any]) -> tuple[float | None, Path | None]:
    if summary_path is not None:
        frac_by_volume = phase_record.get("negative_mode_fraction_by_volume", {})
        if isinstance(frac_by_volume, dict):
            value = frac_by_volume.get("1.0", frac_by_volume.get("1"))
            if value is not None:
                return float(value), summary_path
    reports = [root / "data/processed" / system / phase / "phonon_report.json"]
    reports += sorted(root.glob(f"result/experiments/**/{system}*{phase}*imaginary_diagnosis.json"))
    for path in reports:
        if not path.exists():
            continue
        payload = json.loads(path.read_text(encoding="utf-8"))
        value = payload.get("imaginary_mode_fraction", payload.get("negative_mode_fraction"))
        if value is not None:
            return float(value), path
    return None, None


def _soft_mode_fields(root: Path, system: str, phase: str) -> tuple[Any, Any]:
    if system != "sio2" or phase != "quartz_beta":
        return "", ""
    path = root / "result/experiments/imaginary_modes/quartz_soft_mode.json"
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
        if isinstance(payload.get("heads"), list) and payload["heads"]:
            payload = payload["heads"][0]
        return payload.get("minimum_frequency_THz", ""), payload.get("minimum_frequency_qpoints", "")
    except (OSError, json.JSONDecodeError):
        return "", ""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--output-root", type=Path, default=Path("result/tables"))
    args = parser.parse_args()
    root = args.repo_root.resolve()
    output = args.output_root if args.output_root.is_absolute() else root / args.output_root
    inv = inventory(root)
    rows: list[dict[str, Any]] = []
    missing: list[dict[str, str]] = []
    inputs = [root / "result/experiments/data_prep/raw_inventory/inventory/seed_none/metrics.json"]
    for key, item in sorted(inv["phases"].items()):
        system, phase = str(item["system"]), str(item["phase"])
        phase_dir = root / "data/processed" / system / phase
        for source_path in (phase_dir / "meta.json", phase_dir / "structure.extxyz"):
            if source_path.exists():
                inputs.append(source_path)
        meta = item.get("meta", {})
        representative = meta.get("representative", {}) if isinstance(meta.get("representative"), dict) else {}
        sg = representative.get("space_group", {}) if isinstance(representative.get("space_group"), dict) else {}
        primitive_edges, atom_count = _lattice(phase_dir / "structure.extxyz")
        phase_record, summary_path, summary = _candidate_qh(root, system, phase)
        if summary_path is not None:
            inputs.append(summary_path)
        imaginary_fraction, diagnostic_path = _diagnostic(root, system, phase, summary_path, phase_record)
        minimum_frequency, minimum_qpoints = _soft_mode_fields(root, system, phase)
        soft_mode_path = root / "result/experiments/imaginary_modes/quartz_soft_mode.json"
        if system == "sio2" and phase == "quartz_beta" and soft_mode_path.exists():
            inputs.append(soft_mode_path)
        if diagnostic_path is not None:
            inputs.append(diagnostic_path)
        qh_matrix = phase_record.get("supercell_matrix", summary.get("supercell_matrix", ""))
        if not isinstance(qh_matrix, list):
            qh_matrix = [int(value) for value in str(qh_matrix).replace("x", " ").split()] if qh_matrix else []
        qh_edges = [edge * int(multiplier) for edge, multiplier in zip(primitive_edges, qh_matrix)]
        derived_reliable = "" if imaginary_fraction is None else bool(float(imaginary_fraction) <= IMAGINARY_THRESHOLD)
        source = f"data/processed/{system}/{phase}/structure.extxyz" if representative.get("status") == "relaxed" else representative.get("source", representative.get("source_structure", ""))
        rows.append({
            "system": system,
            "phase": phase,
            "space_group_before": sg.get("before", representative.get("space_group_before", "")),
            "space_group_after": sg.get("after", representative.get("space_group_after", "")),
            "primitive_cell_atoms": 1 if phase == "bcc" else representative.get("primitive_atoms", atom_count or ""),
            "representative_source": source,
            "representative_origin": _representative_origin(representative),
            "relaxation_status": representative.get("status", ""),
            "qh_supercell": "x".join(str(value) for value in qh_matrix) if qh_matrix else "",
            "primitive_shortest_edge_A": min(primitive_edges) if primitive_edges else "",
            "qh_supercell_shortest_edge_A": min(qh_edges) if qh_edges else "",
            "q_mesh": phase_record.get("mesh", summary.get("mesh", "")),
            "imaginary_fraction": imaginary_fraction if imaginary_fraction is not None else "",
            "imaginary_source": str(summary_path.relative_to(root)) if summary_path is not None else (str(diagnostic_path.relative_to(root)) if diagnostic_path is not None else ""),
            "imaginary_volume_scale": "1.0" if summary_path is not None else "not reported",
            "minimum_frequency_THz": minimum_frequency,
            "minimum_frequency_qpoints": json.dumps(minimum_qpoints, separators=(",", ":")) if isinstance(minimum_qpoints, list) else minimum_qpoints,
            "qh_reliable": derived_reliable,
        })
        if not primitive_edges:
            missing.append({"field": f"{key}.primitive_shortest_edge_A", "reason": "structure lattice is not present"})
        if imaginary_fraction is None:
            missing.append({"field": f"{key}.imaginary_fraction", "reason": "QH diagnostic did not record an imaginary-mode fraction"})
        if not qh_matrix:
            missing.append({"field": f"{key}.qh_supercell", "reason": "QH summary is not available"})
    fields = ["system", "phase", "space_group_before", "space_group_after", "primitive_cell_atoms", "representative_source", "representative_origin", "relaxation_status", "qh_supercell", "primitive_shortest_edge_A", "qh_supercell_shortest_edge_A", "q_mesh", "imaginary_fraction", "imaginary_source", "imaginary_volume_scale", "minimum_frequency_THz", "minimum_frequency_qpoints", "qh_reliable"]
    csv_path = output / "phase_inventory.csv"
    csv_write(csv_path, fields, rows)
    meta_write(root, csv_path, inputs, missing, {"imaginary_reliable_fraction_threshold": IMAGINARY_THRESHOLD, "qh_reliable_derived": "imaginary_fraction <= threshold"})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
