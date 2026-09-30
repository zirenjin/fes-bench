"""Run the frozen [2,2,2] metal imaginary-mode reference diagnostics — plan Table 2b."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from fes_bench.qh.diagnose import _diagnose, _frequencies


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--checkpoint", required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    root = args.repo_root.resolve()
    output = args.output_root if args.output_root.is_absolute() else root / args.output_root
    for system in ("hf", "ti", "zr"):
        for phase in ("hcp", "bcc"):
            config = {
                "data_root": str(root / "data/processed"),
                "system": system,
                "phase": phase,
                "checkpoint": args.checkpoint,
                "volume_scale": 1.0,
                "supercell_matrix": [2, 2, 2],
                "mesh": [12, 12, 12],
                "displacement_distance_A": 0.01,
                "negative_threshold_THz": -0.05,
                "gamma_radius_fraction": 0.05,
                "asr": False,
                "fc_symmetry": False,
            }
            qpoints, frequencies, cell = _frequencies(config, root / "data/processed" / system / phase, "Domains_Alloy")
            item = _diagnose(qpoints, frequencies, cell, config, "Domains_Alloy")
            payload = {"system": system, "phase": phase, "head": "Domains_Alloy", "supercell_matrix": [2, 2, 2], "mesh": [12, 12, 12], "diagnostic": item}
            path = output / system / f"{phase}_imaginary_diagnosis.json"
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
            print(f"diagnosed {system}/{phase}: min={item['minimum_frequency_THz']:.6g} THz fraction={item['negative_mode_fraction']:.6g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
