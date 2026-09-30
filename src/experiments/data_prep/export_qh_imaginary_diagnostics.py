"""Export per-phase imaginary-mode diagnostics from a QH run — plan Table 2b."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    args = parser.parse_args()
    summary_path = Path(args.summary).resolve()
    repo_root = args.repo_root.resolve()
    payload = json.loads(summary_path.read_text(encoding="utf-8"))
    output = Path(args.output_root).resolve()
    threshold = -0.05
    for phase, record in payload.get("phases", {}).items():
        frac_by_volume = record.get("negative_mode_fraction_by_volume", {})
        qfrac_by_volume = record.get("negative_qpoint_fraction_by_volume", {})
        min_by_volume = record.get("minimum_frequency_THz_by_volume", {})
        volume = "1.0" if "1.0" in frac_by_volume else sorted(frac_by_volume)[0]
        item = {
            "system": payload.get("system"),
            "phase": phase,
            "head": payload.get("head"),
            "checkpoint": payload.get("checkpoint"),
            "checkpoint_sha256": payload.get("checkpoint_sha256"),
            "mesh": payload.get("mesh"),
            "supercell_matrix": record.get("supercell_matrix", payload.get("supercell_matrix")),
            "volume_scale": float(volume),
            "minimum_frequency_THz": min_by_volume.get(volume),
            "negative_mode_fraction": frac_by_volume.get(volume),
            "negative_qpoint_fraction": qfrac_by_volume.get(volume),
            "negative_threshold_THz": threshold,
            "qh_reliable": frac_by_volume.get(volume) is not None and float(frac_by_volume[volume]) <= 0.01,
            "source_qh_summary": str(summary_path.relative_to(repo_root)) if summary_path.is_relative_to(repo_root) else str(summary_path),
        }
        path = output / f"{phase}_imaginary_diagnosis.json"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(item, indent=2) + "\n", encoding="utf-8")
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
