#!/usr/bin/env bash
# Regenerate every reported metric, table and figure from the stored curves (no model inference, no DFT).
set -euo pipefail
cd "$(dirname "$0")"
python database_pairs/make_main.py
python database_pairs/audit_numerics.py
python finite_temperature/scripts/a_layer_score.py finite_temperature/data data/processed finite_temperature/results/A_LAYER.json
python analyses/curvature_correction/scripts/make_fig5.py analyses/curvature_correction/figures
echo "done: database_pairs/{tables,figures,audit}, finite_temperature/results, analyses/curvature_correction/figures"
