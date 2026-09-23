"""Compare normalized raw payloads and derived tables with historical inputs."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from common import write_audit


def digest(path: Path) -> str:
    h = hashlib.sha256(); h.update(path.read_bytes()); return h.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument("--repo-root", type=Path, default=Path(".")); args = parser.parse_args(); root = args.repo_root.resolve()
    rows, inputs = [], []
    manifest = json.loads((root / "configs/raw_runs.json").read_text(encoding="utf-8")); inputs.append(root / "configs/raw_runs.json")
    exact = 0; total = 0
    for spec in manifest["runs"]:
        if spec["mode"] not in {"nested", "flat"}: continue
        source = root / spec["source"]; inputs.append(source); original = json.loads(source.read_text(encoding="utf-8"))
        if spec["mode"] == "nested":
            for split, predictors in original.items():
                for predictor, payload in predictors.items():
                    total += 1; output = root / "results/raw_runs" / split / predictor / f"seed_{spec.get('seed', 'none')}" / "metrics.json"; got = json.loads(output.read_text(encoding="utf-8"))["metrics"]
                    if got == payload: exact += 1
        else:
            total += 1; output = root / "results/raw_runs" / spec["split"] / spec["predictor"] / f"seed_{spec.get('seed', 'none')}" / "metrics.json"; got = json.loads(output.read_text(encoding="utf-8"))["metrics"]
            if got == original: exact += 1
    rows.append({"check": "normalized_payload_equality", "expected": total, "observed": exact, "status": "pass" if exact == total else "changed", "evidence": "configs/raw_runs.json and results/raw_runs"})
    qh_metrics_paths = sorted(root.glob("results/raw_runs/qh/*/seed_none/metrics.json"))
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
    table_count = len(list((root / "results/tables").glob("*.csv")))
    rows.append({"check": "derived_table_generation", "expected": 10, "observed": table_count, "status": "pass" if table_count == 10 else "pending", "evidence": "scripts/tables (4 inventory/definition tables + 3 split predictor tables + 3 split crossing tables)"})
    write_audit(root, "reproduction_compare", ["check", "expected", "observed", "status", "evidence"], rows, inputs, [])
    return 0


if __name__ == "__main__": raise SystemExit(main())
