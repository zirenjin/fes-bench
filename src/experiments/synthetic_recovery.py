"""Summarize synthetic coefficient recovery — plan experiment E3."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from _common import write_audit


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument("--repo-root", type=Path, default=Path(".")); args = parser.parse_args(); root = args.repo_root.resolve()
    path = root / "result/experiments/synthetic_recovery/raw_runs/synthetic_recovery/seed_none/metrics.json"; payload = json.loads(path.read_text(encoding="utf-8"))["metrics"]
    rows = [{"predictor": name, "basis": item.get("basis"), "MAE_eV_per_atom": item.get("MAE_eV_per_atom"), "coefficient_relative_error": json.dumps(item.get("coefficient_relative_error"), ensure_ascii=False)} for name, item in payload.get("fits", {}).items()]
    write_audit(root, "synthetic_recovery", ["predictor", "basis", "MAE_eV_per_atom", "coefficient_relative_error"], rows, [path], [])
    return 0


if __name__ == "__main__": raise SystemExit(main())
