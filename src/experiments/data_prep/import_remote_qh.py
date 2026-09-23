"""Import completed QH artifacts from an isolated run without recomputation."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sanitize(value):
    if isinstance(value, dict): return {key: sanitize(item) for key, item in value.items()}
    if isinstance(value, list): return [sanitize(item) for item in value]
    if isinstance(value, str) and (value.startswith("/") or ":/" in value):
        prefix = ""
        if ":/" in value and not value.startswith("/"):
            prefix, value = value.split(":/", 1); prefix += ":"
        return f"{prefix}external/{Path(value).name}"
    return value


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True, help="Local mirror of the isolated QH output directory")
    parser.add_argument("--system", required=True)
    parser.add_argument("--output-root", type=Path, default=Path("result/experiments/quasi_harmonic/raw_runs"))
    parser.add_argument("--source-uri", default="remote://isolated-qh-run")
    args = parser.parse_args(); root = Path(".").resolve(); source = args.source_root.resolve(); output = args.output_root / args.system / "seed_none"; output.mkdir(parents=True, exist_ok=True)
    if not source.exists(): raise FileNotFoundError(source)
    files = sorted(source.glob("*.csv")) + sorted(source.glob("*.json")); missing = []
    copied = []
    for path in files:
        target = output / path.name
        if path.suffix == ".json": target.write_text(json.dumps(sanitize(json.loads(path.read_text(encoding="utf-8"))), indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        else: shutil.copy2(path, target)
        copied.append(path)
    summary = output / "qh_summary.json"
    if not summary.exists(): missing.append({"field": "qh_summary.json", "reason": "source mirror has no qh_summary.json"})
    for phase in json.loads(summary.read_text(encoding="utf-8")).get("phases", {}) if summary.exists() else {}:
        csv_name = f"{phase}_fqh.csv"
        if not (output / csv_name).exists(): missing.append({"field": csv_name, "reason": "source mirror has no phase F_QH CSV"})
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    metrics = {"provenance": {"git_commit": commit, "split": "qh", "split_sha256": None, "predictor_config": None, "predictor_config_sha256": None, "checkpoint_path": "external/checkpoints/unknown", "checkpoint_sha256": None, "evaluator_version": "fes-bench QH raw-run import v1", "timestamp_utc": datetime.now(timezone.utc).isoformat(), "source_uri": args.source_uri, "source_sha256": {path.name: sha256(path) for path in copied}, "missing": missing}, "metrics": {"files": [path.name for path in copied]}}
    (output / "metrics.json").write_text(json.dumps(metrics, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(output)
    return 0


if __name__ == "__main__": raise SystemExit(main())
