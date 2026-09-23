"""Audit split freezing, calibration windows, and read-only split usage — plan experiment E10."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from _common import write_audit


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument("--repo-root", type=Path, default=Path(".")); args = parser.parse_args(); root = args.repo_root.resolve()
    rows, inputs, missing = [], [], []
    for path in sorted((root / "data/processed/splits").glob("*.json")):
        split = json.loads(path.read_text(encoding="utf-8")); inputs.append(path); folds = split.get("folds", {"all": split})
        for fold, body in sorted(folds.items()):
            train = {(item["system"], item["phase"], item["T_index"]) for item in body.get("train", [])}; test = {(item["system"], item["phase"], item["T_index"]) for item in body.get("test", [])}; rows.append({"split": path.stem, "fold": fold, "train_frames": len(train), "test_frames": len(test), "train_test_overlap": len(train & test), "git_commit": split.get("git_commit", ""), "sha256": split.get("sha256", "")})
    e2 = root / "result/experiments/calibration_window/raw_runs/calibration_window/seed_none/metrics.json"
    if e2.exists(): inputs.append(e2); calibration_rows = len(json.loads(e2.read_text(encoding="utf-8"))["metrics"].get("rows", []))
    else: calibration_rows = 0; missing.append({"field": "calibration_window", "reason": "normalized E2 raw run missing"})
    rows.append({"split": "e2_calibration_window", "fold": "all", "train_frames": "", "test_frames": "", "train_test_overlap": "", "git_commit": "", "sha256": f"rows={calibration_rows}"})
    write_audit(root, "leakage", ["split", "fold", "train_frames", "test_frames", "train_test_overlap", "git_commit", "sha256"], rows, inputs, missing)
    return 0


if __name__ == "__main__": raise SystemExit(main())
