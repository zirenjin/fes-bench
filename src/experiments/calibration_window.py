"""Compare train-window and held-out calibration — plan experiment E2."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from _common import write_audit


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument("--repo-root", type=Path, default=Path(".")); args = parser.parse_args(); root = args.repo_root.resolve()
    rows, inputs = [], []
    for path in sorted((root / "result/experiments/calibration_window/raw_runs/calibration_window").glob("seed_*/metrics.json")):
        inputs.append(path); payload = json.loads(path.read_text(encoding="utf-8"))["metrics"]
        for item in payload.get("rows", []): rows.append(item)
    fields = ["basis", "seed", "formal_G_MAE", "heldout_G_MAE", "delta_G_MAE", "relative_delta_percent", "mean_pair_delta_G_MAE", "max_Tc_delta_K"]
    write_audit(root, "calibration_window", fields, rows, inputs, [])
    return 0


if __name__ == "__main__": raise SystemExit(main())
