"""Compare normalized raw payloads and derived tables with historical inputs."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

from _common import write_audit


def digest(path: Path) -> str:
    h = hashlib.sha256(); h.update(path.read_bytes()); return h.hexdigest()


def raw_output(root: Path, split: str, predictor: str, seed: str) -> Path:
    if split == "e1_full_grid":
        return root / "result/experiments/crossing_reevaluation/raw_runs" / predictor / f"seed_{seed}" / "metrics.json"
    if split == "e2_calibration_window":
        return root / "result/experiments/calibration_window/raw_runs" / predictor / f"seed_{seed}" / "metrics.json"
    if split == "e3_synthetic_recovery":
        return root / "result/experiments/synthetic_recovery/raw_runs" / predictor / f"seed_{seed}" / "metrics.json"
    if split == "reference_statistics":
        return root / "result/experiments/reference_statistics/raw_runs" / predictor / f"seed_{seed}" / "metrics.json"
    return root / "result/experiments/external_baselines" / f"raw_{split}" / predictor / f"seed_{seed}" / "metrics.json"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--baseline-rev", default="731606a", help="pre-migration revision used only for CSV equality auditing")
    args = parser.parse_args(); root = args.repo_root.resolve()
    rows, inputs = [], []
    manifest = json.loads((root / "configs/raw_runs.json").read_text(encoding="utf-8")); inputs.append(root / "configs/raw_runs.json")
    exact = 0; total = 0
    for spec in manifest["runs"]:
        if spec["mode"] not in {"nested", "flat"}: continue
        source = root / spec["source"]; inputs.append(source); original = json.loads(source.read_text(encoding="utf-8"))
        if spec["mode"] == "nested":
            for split, predictors in original.items():
                for predictor, payload in predictors.items():
                    total += 1; output = raw_output(root, split, predictor, str(spec.get("seed", "none"))); got = json.loads(output.read_text(encoding="utf-8"))["metrics"]
                    if got == payload: exact += 1
        else:
            total += 1; output = raw_output(root, str(spec["split"]), str(spec["predictor"]), str(spec.get("seed", "none"))); got = json.loads(output.read_text(encoding="utf-8"))["metrics"]
            if got == original: exact += 1
    rows.append({"check": "normalized_payload_equality", "expected": total, "observed": exact, "status": "pass" if exact == total else "changed", "evidence": "configs/raw_runs.json and result/experiments raw outputs"})
    qh_metrics_paths = sorted(root.glob("result/experiments/quasi_harmonic/raw_runs/*/seed_none/metrics.json"))
    for qh_metrics in qh_metrics_paths:
        inputs.append(qh_metrics)
        system = qh_metrics.parts[-3]
        qh_payload = json.loads(qh_metrics.read_text(encoding="utf-8")); qh_checks = []
    # CSVs are byte-preserving imports and can be checked against their source
    # hashes. JSON support files are sanitized to remove machine paths, so their
    # local bytes are not expected to match the remote bytes.
        for name, expected in qh_payload["provenance"].get("source_sha256", {}).items():
            if not name.endswith("_fqh.csv"):
                continue
            candidate = qh_metrics.parent / name
            qh_checks.append(candidate.exists() and digest(candidate) == expected)
        rows.append({"check": f"{system}_qh_csv_source_sha256", "expected": len(qh_checks), "observed": sum(qh_checks), "status": "pass" if qh_checks and all(qh_checks) else "pending", "evidence": str(qh_metrics.relative_to(root))})
        json_names = [name for name in qh_payload["provenance"].get("source_sha256", {}) if name.endswith(".json")]
        json_present_names = [name for name in json_names if (qh_metrics.parent / name).exists()]
        json_missing_names = [name for name in json_names if name not in json_present_names]
        rows.append({"check": f"{system}_qh_json_support_present", "expected": len(json_names) - len(json_missing_names), "observed": len(json_present_names), "status": "pass" if json_present_names and len(json_present_names) == len(json_names) - len(json_missing_names) else "pending", "evidence": "sanitized local JSON with source SHA recorded in metrics provenance"})
        missing_paths = {item.get("path") for item in qh_payload["provenance"].get("missing", []) if isinstance(item, dict)}
        rows.append({"check": f"{system}_qh_missing_support_documented", "expected": len(json_missing_names), "observed": len([name for name in json_missing_names if name in missing_paths]), "status": "pass" if all(name in missing_paths for name in json_missing_names) else "pending", "evidence": f"provenance.missing in {qh_metrics.relative_to(root)}"})
    table_count = len(list((root / "result/tables").glob("*.csv")))
    rows.append({"check": "derived_table_generation", "expected": 10, "observed": table_count, "status": "pass" if table_count == 10 else "pending", "evidence": "src/tables (4 inventory/definition tables + 3 split predictor tables + 3 split crossing tables)"})
    equal_tables = 0
    for table in sorted((root / "result/tables").glob("*.csv")):
        try:
            prior = subprocess.check_output(
                ["git", "-C", str(root), "show", f"{args.baseline_rev}:results/tables/{table.name}"],
            )
        except subprocess.CalledProcessError:
            continue
        if table.read_bytes() == prior:
            equal_tables += 1
    rows.append({"check": "pre_migration_table_csv_equality", "expected": table_count, "observed": equal_tables, "status": "pass" if table_count == equal_tables == 10 else "changed", "evidence": f"{args.baseline_rev}:results/tables compared byte-for-byte"})
    write_audit(root, "reproduction_compare", ["check", "expected", "observed", "status", "evidence"], rows, inputs, [])
    return 0 if table_count == equal_tables == 10 else 1


if __name__ == "__main__": raise SystemExit(main())
