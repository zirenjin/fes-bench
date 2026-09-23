"""Collect per-phase imaginary-mode diagnostics — plan experiment E7."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from common import write_audit


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument("--repo-root", type=Path, default=Path(".")); args = parser.parse_args(); root = args.repo_root.resolve()
    rows, inputs, missing = [], [], []
    candidates = sorted((root / "results").glob("phase2/**/*imaginary_diagnosis.json")) + sorted((root / "results/raw_runs/qh").glob("*/seed_none/*_imaginary_diagnosis.json"))
    for path in candidates:
        inputs.append(path); payload = json.loads(path.read_text(encoding="utf-8")); phase = path.stem.removesuffix("_imaginary_diagnosis")
        rows.append({"phase": phase, "minimum_frequency_THz": payload.get("minimum_frequency_THz", payload.get("min_frequency_THz", "")), "negative_mode_fraction": payload.get("negative_mode_fraction", payload.get("negative_fraction", "")), "negative_qpoint_fraction": payload.get("negative_qpoint_fraction", ""), "qh_reliable": payload.get("qh_reliable", False), "source": str(path.relative_to(root))})
    if not rows: missing.append({"field": "imaginary_modes", "reason": "no diagnosis JSON found"})
    write_audit(root, "imaginary_modes", ["phase", "minimum_frequency_THz", "negative_mode_fraction", "negative_qpoint_fraction", "qh_reliable", "source"], rows, inputs, missing)
    return 0


if __name__ == "__main__": raise SystemExit(main())
