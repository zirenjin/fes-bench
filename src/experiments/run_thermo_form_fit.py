"""Generate unified-evaluator metrics for the training-only thermo-form floor."""

from __future__ import annotations

import json
import argparse
from datetime import datetime, timezone
from pathlib import Path

from _run_trivial_floor import _evaluate, _json_safe


ROOT = Path(__file__).resolve().parents[2]
DATA_ROOT = ROOT / "data/processed"
SPLITS = {
    "temp_extrap": ROOT / "data/processed/splits_v2/temp_extrap.json",
    "phase_lopo": ROOT / "data/processed/splits_v2/phase_lopo.json",
    "system_loso": ROOT / "data/processed/splits_v2/system_loso.json",
}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--split", choices=[*SPLITS, "all"], default="all")
    args = parser.parse_args()
    selected = SPLITS if args.split == "all" else {args.split: SPLITS[args.split]}
    for split_name, split_path in selected.items():
        split = json.loads(split_path.read_text(encoding="utf-8"))
        metrics = _evaluate(split, DATA_ROOT, "thermo_form_fit")
        output = ROOT / "result/experiments/external_baselines" / f"raw_{split_name}" / "thermo_form_fit" / "seed_none" / "metrics.json"
        output.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "provenance": {
                "split": split_name,
                "split_path": str(split_path.relative_to(ROOT)),
                "predictor": "thermo-form fit",
                "basis": "a + bT + cT log(T)",
                "fit_scope": "shared frozen training rows only",
                "timestamp_utc": datetime.now(timezone.utc).isoformat(),
            },
            "metrics": _json_safe(metrics),
        }
        output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
