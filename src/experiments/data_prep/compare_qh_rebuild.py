"""Compare old and rebuilt representative/QH diagnostics — plan Table 2b."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_fqh(path: Path) -> dict[float, float]:
    with path.open(encoding="utf-8", newline="") as handle:
        return {float(row["T_K"]): float(row["F_QH_eV_per_atom"]) for row in csv.DictReader(handle)}


def _atom_count(path: Path) -> int | None:
    try:
        return int(path.read_text(encoding="utf-8").splitlines()[0].strip())
    except (OSError, ValueError, IndexError):
        return None


def _nearest_common(old: dict[float, float], new: dict[float, float]) -> list[float]:
    return sorted(set(old).intersection(new))


def _phase_row(root: Path, system: str, phase: str, old_fqh: Path, old_diag: Path | None, new_fqh: Path, new_summary: Path, old_meta: Path, new_meta: Path) -> dict[str, Any]:
    old_values, new_values = _read_fqh(old_fqh), _read_fqh(new_fqh)
    common = _nearest_common(old_values, new_values)
    differences = [abs(new_values[t] - old_values[t]) for t in common]
    old_meta_payload = json.loads(old_meta.read_text(encoding="utf-8"))
    new_meta_payload = json.loads(new_meta.read_text(encoding="utf-8"))
    old_energy = old_meta_payload.get("representative", {}).get("energy_eV_per_atom", {}).get("after")
    new_energy = new_meta_payload.get("formal_energy_eV_per_atom")
    if new_energy is None:
        new_energy = new_meta_payload.get("representative", {}).get("energy_eV_per_atom", {}).get("after")
    old_fraction = None
    if old_diag and old_diag.exists():
        old_diag_payload = json.loads(old_diag.read_text(encoding="utf-8"))
        old_fraction = old_diag_payload.get("negative_mode_fraction", old_diag_payload.get("diagnostic", {}).get("negative_mode_fraction"))
    new_summary_payload = json.loads(new_summary.read_text(encoding="utf-8"))
    phase_summary = new_summary_payload["phases"][phase]
    new_fraction = phase_summary.get("negative_mode_fraction_by_volume", {}).get("1.0")
    return {
        "system": system,
        "phase": phase,
        "old_space_group": old_meta_payload.get("representative", {}).get("space_group", {}).get("after"),
        "new_space_group": new_meta_payload.get("representative", {}).get("space_group", {}).get("after"),
        "old_atoms": old_meta_payload.get("representative", {}).get("source_atoms"),
        "new_atoms": new_meta_payload.get("representative", {}).get("source_atoms"),
        "old_E_DPA_eV_per_atom": old_energy,
        "new_E_DPA_eV_per_atom": new_energy,
        "delta_E_DPA_eV_per_atom": None if old_energy is None or new_energy is None else float(new_energy) - float(old_energy),
        "old_imaginary_fraction": old_fraction,
        "new_imaginary_fraction": new_fraction,
        "delta_imaginary_fraction": None if old_fraction is None or new_fraction is None else float(new_fraction) - float(old_fraction),
        "common_temperature_points": len(common),
        "mean_abs_delta_F_QH_eV_per_atom": sum(differences) / len(differences) if differences else None,
        "max_abs_delta_F_QH_eV_per_atom": max(differences) if differences else None,
        "old_fqh_sha256": sha256(old_fqh),
        "new_fqh_sha256": sha256(new_fqh),
        "old_checkpoint_sha256": "unknown_legacy",
        "old_head": old_diag_payload.get("head", "unknown_legacy") if old_diag and old_diag.exists() else "unknown_legacy",
        "new_checkpoint_sha256": new_summary_payload.get("checkpoint_sha256", ""),
        "new_head": new_summary_payload.get("head", ""),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--system", required=True)
    parser.add_argument("--new-root", required=True, help="result experiment directory containing raw_runs/<system>")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    root = args.repo_root.resolve()
    system = args.system
    new_dir = root / args.new_root / "raw_runs" / system
    new_summary = new_dir / "qh_summary.json"
    rows = []
    source_paths: list[Path] = [new_summary]
    for phase_dir in sorted((root / "data/processed" / system).iterdir()):
        phase = phase_dir.name
        old_fqh = phase_dir / "fqh.csv"
        if not old_fqh.exists():
            old_fqh = root / "result/experiments/quasi_harmonic/raw_runs" / system / "seed_none" / f"{phase}_fqh.csv"
        new_fqh = new_dir / f"{phase}_fqh.csv"
        old_meta = phase_dir / "meta.json"
        archived_meta = phase_dir.parent / "_archive_md_snapshot" / phase / "meta.json"
        if archived_meta.exists():
            old_meta = archived_meta
        new_meta = phase_dir / "meta.json"
        if not (old_fqh.exists() and new_fqh.exists() and old_meta.exists()):
            continue
        old_diag = root / "result/experiments/quasi_harmonic/raw_runs" / system / "seed_none" / f"{phase}_imaginary_diagnosis.json"
        if not old_diag.exists():
            old_diag = root / "result/experiments/quasi_harmonic_10a/old_diagnostics" / system / f"{phase}_imaginary_diagnosis.json"
        row = _phase_row(root, system, phase, old_fqh, old_diag if old_diag.exists() else None, new_fqh, new_summary, old_meta, new_meta)
        row["old_atoms"] = row["old_atoms"] or _atom_count(old_meta.parent / "structure.extxyz")
        row["new_atoms"] = row["new_atoms"] or _atom_count(phase_dir / "structure.extxyz")
        rows.append(row)
        source_paths.extend([old_fqh, new_fqh, old_meta, new_meta])
        if old_diag.exists():
            source_paths.append(old_diag)
    output = root / args.output
    output.mkdir(parents=True, exist_ok=True)
    findings = output / "findings.csv"
    fields = list(rows[0]) if rows else ["system", "phase"]
    with findings.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)
    inputs = [new_summary] + [Path(row_path) for row_path in []]
    meta = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "git_commit": subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip(),
        "system": system,
        "old_sources": "canonical data/processed/<system>/<phase>/fqh.csv and archived/current metadata",
        "new_source": str(new_summary.relative_to(root)),
        "inputs_sha256": {str(path.relative_to(root)): sha256(path) for path in sorted(set(source_paths)) if path.exists()},
        "missing": [] if rows else [{"field": "findings", "reason": "no matching old/new F_QH files"}],
    }
    (output / "findings.meta.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    metrics = {
        "provenance": {
            "git_commit": meta["git_commit"],
            "split": "qh",
            "split_sha256": None,
            "predictor_config": f"configs/qh/{system}_qh_10a.yaml" if system != "sio2" else "configs/qh/sio2_sse_pbe_10a.yaml",
            "predictor_config_sha256": sha256(root / (f"configs/qh/{system}_qh_10a.yaml" if system != "sio2" else "configs/qh/sio2_sse_pbe_10a.yaml")),
            "checkpoint_path": "<external-checkpoint>/DPA-3.1-3M.pt",
            "checkpoint_sha256": "86dd3a804d78ca5d203ebf98747e8f16dff9713ba8950097ceb760b161e19907",
            "evaluator_version": "compare_qh_rebuild.py v1",
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
            "raw_summary": str(new_summary.relative_to(root)),
        },
        "metrics": {"rows": rows, "old_reference": "pre-rebuild canonical QH outputs", "new_reference": str(new_summary.relative_to(root))},
    }
    (output / "metrics.json").write_text(json.dumps(metrics, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    (new_dir / "metrics.json").write_text(json.dumps(metrics, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
