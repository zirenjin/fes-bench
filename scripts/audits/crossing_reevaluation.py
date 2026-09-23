"""Re-evaluate frozen crossing records — plan experiment E1."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from common import write_audit


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument("--repo-root", type=Path, default=Path(".")); args = parser.parse_args(); root = args.repo_root.resolve()
    rows, inputs = [], []
    for path in sorted((root / "results/raw_runs/e1_full_grid").glob("*/seed_*/metrics.json")):
        inputs.append(path); record = json.loads(path.read_text(encoding="utf-8"))["metrics"]; basis = path.parts[-3]; seed = path.parts[-2].removeprefix("seed_")
        metrics = record.get("metrics", {}).get("formal_train_window", {}).get("full_range", {})
        for pair, item in metrics.get("pairs", {}).items():
            for index, error in enumerate(item.get("Tc_error_K", [])):
                rows.append({"basis": basis, "seed": seed, "pair": pair, "crossing_index": index + 1, "reference_Tc_K": item.get("reference_Tc_K", [])[index], "Tc_error_K": error, "Tc_err_from_dG_K": item.get("Tc_err_from_dG_K", [])[index]})
    write_audit(root, "crossing_reevaluation", ["basis", "seed", "pair", "crossing_index", "reference_Tc_K", "Tc_error_K", "Tc_err_from_dG_K"], rows, inputs, [])
    return 0


if __name__ == "__main__": raise SystemExit(main())
