#!/usr/bin/env bash
# Rebuild migration-derived experiment findings and all CSV tables; never trains a model.
set -euo pipefail

repo_root="$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)"
cd "$repo_root"
export PYTHONPATH="src/lib:src/experiments:src/tables${PYTHONPATH:+:$PYTHONPATH}"

python3 src/experiments/data_prep/normalize_raw_runs.py --manifest configs/raw_runs.json
for experiment in crossing_reevaluation calibration_window synthetic_recovery skill_floor mae_convention split_sign_structure imaginary_modes representative_structures reference_statistics leakage; do
  python3 "src/experiments/${experiment}.py"
done
for table in system_inventory phase_inventory split_definitions metric_definitions; do
  python3 "src/tables/${table}.py"
done
for split in temp_extrap phase_lopo system_loso; do
  python3 src/tables/predictor_comparison.py --split "$split"
  python3 src/tables/crossing_errors.py --split "$split"
done
python3 src/experiments/reproduction_compare.py
python3 src/experiments/data_prep/index_experiment_metrics.py
