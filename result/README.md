# Result index

`result/` is the sole generated-evidence layer. Every table has a CSV for
review, a sibling `.meta.json` with commit and input hashes, and a producer in
`src/tables/`. Every experiment directory has `findings.csv`,
`findings.meta.json`, and a top-level `metrics.json` provenance index; raw
metrics remain below its `raw_runs/` directory when an experiment has them.

| Conclusion or output | Command | Source evidence |
| --- | --- | --- |
| system / phase / split / metric inventories | `python3 src/tables/{system,phase,split,metric}_inventory.py` | `experiments/data_prep/raw_inventory/` |
| predictor comparison and Tc errors | `python3 src/tables/{predictor_comparison,crossing_errors}.py --split <split>` | `experiments/external_baselines/raw_<split>/` |
| E1 / E2 / E3 | `python3 src/experiments/{crossing_reevaluation,calibration_window,synthetic_recovery}.py` | matching experiment directory |
| floor and skill audit | `python3 src/experiments/skill_floor.py` | `experiments/skill_floor/` |
| QH and imaginary-mode records | `python3 src/experiments/{quasi_harmonic,imaginary_modes}.py` | `experiments/quasi_harmonic/` |
| MAE, split-sign, representative, reference, leakage audits | `python3 src/experiments/<name>.py` | matching experiment directory |

All paths are relative to the repository root. Rebuild the complete derived
layer with `bash src/reproduce_all.sh`.

The Tc-error value for `interp_const` in
`tables/crossing_errors_temp_extrap.csv` comes from
`experiments/external_baselines/raw_temp_extrap/interp_const/seed_none/metrics.json`.
That file's `provenance` block records the imported source hash and frozen split
hash; `src/tables/crossing_errors.py --split temp_extrap` performs the read.

## N/A vocabulary

- `n/a:no_reference_crossing`: the frozen reference grid has no sign change.
- `n/a:pair_only_predictor`: the predictor does not produce phase-level `G`.
- `n/a:no_training_phase`: a held-out fold lacks a required training phase or system.

The matching `.meta.json` files retain field-specific explanations. The current
tables have 33, 15, and 32 cells of these respective types.

## Known issue recorded, not resolved

`constant_delta_g` is structurally unavailable on LOPO/LOSO. Therefore
temp-extrap skill uses `constant_delta_g`, whereas LOPO/LOSO use
`global_mean_delta_g`; no denominator was silently unified in this migration.

`_legacy/` contains six cleaned historical delivery records and `CLEANING.md`.
They are tracked for provenance but not read by table scripts.
