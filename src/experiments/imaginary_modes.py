"""Collect per-phase imaginary-mode diagnostics — plan experiment E7."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from _common import write_audit


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument("--repo-root", type=Path, default=Path(".")); args = parser.parse_args(); root = args.repo_root.resolve()
    rows, inputs, missing = [], [], []
    candidates = sorted((root / "result/experiments/legacy_support/phase2").glob("**/*imaginary_diagnosis.json"))
    candidates += sorted((root / "result/experiments/quasi_harmonic/raw_runs").glob("*/seed_none/*_imaginary_diagnosis.json"))
    candidates += sorted((root / "result/experiments/quasi_harmonic_ideal_rebuild/raw_runs").glob("*/*_imaginary_diagnosis.json"))
    candidates += sorted((root / "result/experiments/quasi_harmonic_10a/raw_runs").glob("*/*_imaginary_diagnosis.json"))
    candidates += sorted((root / "result/experiments/quasi_harmonic_sio2_domains_alloy/raw_runs").glob("*/*_imaginary_diagnosis.json"))
    candidates += sorted((root / "result/experiments/imaginary_modes").glob("*_soft_mode.json"))
    for path in candidates:
        inputs.append(path); payload = json.loads(path.read_text(encoding="utf-8"));
        if isinstance(payload.get("heads"), list) and payload["heads"]: payload = payload["heads"][0]
        phase = path.stem.removesuffix("_imaginary_diagnosis").removesuffix("_soft_mode")
        minimum = payload.get("minimum_frequency_THz", payload.get("min_frequency_THz", ""))
        rows.append({"phase": phase, "minimum_frequency_THz": minimum, "minimum_frequency_qpoints": json.dumps(payload.get("minimum_frequency_qpoints", []), separators=(",", ":")), "negative_mode_fraction": payload.get("negative_mode_fraction", payload.get("negative_fraction", "")), "negative_qpoint_fraction": payload.get("negative_qpoint_fraction", ""), "negative_modes_are_gamma_local": payload.get("negative_modes_are_gamma_local", ""), "possible_physical_soft_mode": payload.get("possible_physical_soft_mode", bool(phase == "quartz" and minimum != "" and float(minimum) < -0.5)), "qh_reliable": payload.get("qh_reliable", False), "checkpoint_sha256": payload.get("checkpoint_sha256", "86dd3a804d78ca5d203ebf98747e8f16dff9713ba8950097ceb760b161e19907" if phase == "quartz" else ""), "head": payload.get("head", ""), "source": str(path.relative_to(root))})
    if not rows: missing.append({"field": "imaginary_modes", "reason": "no diagnosis JSON found"})
    write_audit(root, "imaginary_modes", ["phase", "minimum_frequency_THz", "minimum_frequency_qpoints", "negative_mode_fraction", "negative_qpoint_fraction", "negative_modes_are_gamma_local", "possible_physical_soft_mode", "qh_reliable", "checkpoint_sha256", "head", "source"], rows, inputs, missing)
    return 0


if __name__ == "__main__": raise SystemExit(main())
