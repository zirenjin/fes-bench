"""Write provenance indexes at each result/experiments/<name>/metrics.json — migration support."""

from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_commit(root: Path) -> str:
    return subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip()


def main() -> int:
    root = Path(__file__).resolve().parents[3]
    experiments = root / "result" / "experiments"
    for directory in sorted(path for path in experiments.iterdir() if path.is_dir() and path.name != "legacy_support"):
        inputs = sorted(path for path in directory.rglob("metrics.json") if path != directory / "metrics.json")
        payload = {
            "provenance": {
                "git_commit": git_commit(root),
                "split_sha256": None,
                "predictor_config": None,
                "predictor_config_sha256": None,
                "checkpoint_path": None,
                "checkpoint_sha256": None,
                "evaluator_version": "fes-bench experiment-index v1",
                "timestamp_utc": datetime.now(timezone.utc).isoformat(),
                "inputs": [{"path": str(path.relative_to(root)), "sha256": sha256(path)} for path in inputs],
            },
            "metrics": {"kind": "provenance index", "raw_metrics": [str(path.relative_to(directory)) for path in inputs]},
        }
        (directory / "metrics.json").write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
