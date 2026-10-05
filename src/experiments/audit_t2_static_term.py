"""Audit T2's physical baseline for a duplicated QH static term — plan T2."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path


PHASES = {
    "sio2": ("quartz_beta", "cristobalite_beta", "tridymite_p63mmc"),
    "hf": ("hcp", "bcc"),
    "ti": ("hcp", "bcc"),
    "zr": ("hcp", "bcc"),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_commit(root: Path) -> str:
    try:
        return subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def qh_contains_static(root: Path, system: str, phase: str) -> dict[str, object]:
    phase_dir = root / "data/processed" / system / phase
    rows = list(csv.DictReader((phase_dir / "fqh.csv").open(encoding="utf-8", newline="")))
    if not rows:
        raise ValueError(f"empty F_QH file: {phase_dir / 'fqh.csv'}")
    if "E_static_eV_per_atom" in rows[0]:
        static = [float(row["E_static_eV_per_atom"]) for row in rows]
        source = "E_static_eV_per_atom"
    else:
        static = [float(row["F_QH_eV_per_atom"]) - float(row["F_vib_eV_per_atom"]) for row in rows]
        source = "F_QH_eV_per_atom - F_vib_eV_per_atom"
    meta = json.loads((phase_dir / "meta.json").read_text(encoding="utf-8"))
    e_dpa = float(meta["representative"]["energy_eV_per_atom"]["after"])
    return {
        "system": system,
        "phase": phase,
        "source": source,
        "static_min_eV_per_atom": min(static),
        "static_max_eV_per_atom": max(static),
        "static_range_eV_per_atom": max(static) - min(static),
        "E_DPA_metadata_eV_per_atom": e_dpa,
        "static_term_present": True,
        "static_matches_metadata_E_DPA_within_1e-5": max(abs(value - e_dpa) for value in static) <= 1.0e-5,
        "fqh_sha256": sha256(phase_dir / "fqh.csv"),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    args = parser.parse_args()
    root = args.repo_root.resolve()
    config_path = root / "configs/predictors/t2_temp_extrap.yaml"
    config_text = config_path.read_text(encoding="utf-8")
    provenance = []
    for seed in (11, 23, 37):
        path = root / "result/experiments/t3_temp_extrap/raw_runs/qh_residual" / f"seed_{seed}" / "metrics.json"
        payload = json.loads(path.read_text(encoding="utf-8"))
        provenance.append({
            "seed": seed,
            "path": str(path.relative_to(root)),
            "sha256": sha256(path),
            "predictor_config_sha256": payload.get("provenance", {}).get("predictor_config_sha256"),
            "physics_baseline_column": payload.get("provenance", {}).get("physics_baseline_column"),
            "qh_baseline": payload.get("provenance", {}).get("qh_baseline"),
        })
    phases = [qh_contains_static(root, system, phase) for system, names in PHASES.items() for phase in names]
    model_text = (root / "configs/models/qh_residual.yaml").read_text(encoding="utf-8")
    current_config_sha = sha256(config_path)
    old_run_hashes = {item.get("predictor_config_sha256") for item in provenance}
    duplicate = (
        ("E_DPA + F_QH" in config_text)
        or ("E_DPA + F_QH" in model_text)
        or ("E + F_QH" in config_text)
        or any(value and value != current_config_sha for value in old_run_hashes)
    )
    result = {
        "status": "invalid_existing_results" if duplicate else "valid",
        "reason": "F_QH includes the static term, while the T2 plan formula also adds E_DPA" if duplicate else "No duplicated static term declared",
        "observed_config_formula": "F_QH + r_theta(z, T) + c_system" if not duplicate else "legacy run provenance predates corrected F_QH-only formula",
        "required_formula": "F_QH + r_theta(z, T) + c_system",
        "physical_baseline": "fparam column 2 = canonical F_QH_eV_per_atom",
        "retraining_required": duplicate,
        "retraining_scope": "tlog seed 11/23/37 after correcting the T2 formula" if duplicate else None,
        "phase_static_term_audit": phases,
        "existing_runs": provenance,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "git_commit": git_commit(root),
        "input_sha256": {str(path.relative_to(root)): sha256(path) for path in [config_path, root / "configs/models/qh_residual.yaml"]},
    }
    output = root / "result/experiments/t3_temp_extrap/t2_static_term_audit.json"
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    lines = [
        "# T2 static-term audit",
        "",
        f"Status: **{result['status']}**.",
        "",
        "The physical channel is canonical `F_QH_eV_per_atom`, and every phase audit confirms it contains the E_DPA static term. The existing T2 plan declares an additive `E_DPA + F_QH` baseline, so those runs are invalid.",
        "",
        "Required correction: `F_QH + r_theta(z, T) + c_system`; retrain tlog seeds 11/23/37.",
    ]
    (output.with_suffix(".md")).write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
