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
- `n/a:degenerate_prediction`: predicted ΔG is within `1e-12 eV/atom` of zero on the entire evaluation grid, so crossings and their derived errors are undefined.
- `n/a:input_unavailable`: a predictor's required immutable input is absent.

The matching `.meta.json` files retain field-specific explanations. The current
tables have 33, 15, and 32 cells of these respective types.

## E1 reproduction tolerance

The fixed policy is in `configs/reproduction_tolerance.yaml`; it is derived
from resolution, not fitted to a reproduction result. Crossing counts, false/
missed crossings, list cardinalities, `sign_accuracy`, and other ratios must
match exactly. `T_c` and `T_c`-error fields use an absolute tolerance of
0.05 K (one twentieth of the 1 K grid); G/ΔG and calibration-energy fields use
1e-6 eV/atom (one thousandth of the reported meV/atom precision). GPU
floating-point evaluation is not required to be bitwise identical across
runtime environments.

The canonical E1 checkpoint output is
`experiments/crossing_reevaluation/sio2_reanalysis.json`, produced by
`src/experiments/crossing_reevaluation_run.py`. Its field-level audit is
`experiments/crossing_reevaluation/reproduction_check.csv`; the canonical and
historical SHA-256 values, and the pass result, are recorded in its sibling
`sio2_reanalysis.meta.json`. The original input remains unchanged at
`experiments/legacy_support/phase5/sio2_reanalysis.json`.

E1 status: complete; acceptance criterion three did not pass because the
cristobalite–tridymite pair has one false crossing for each of 10/10
checkpoints. This is a model conclusion, not a reproduction failure.

## Known issue recorded, not resolved

`constant_delta_g` is structurally unavailable on LOPO/LOSO. Therefore
temp-extrap skill uses `constant_delta_g`, whereas LOPO/LOSO use
`global_mean_delta_g`; no denominator was silently unified in this migration.

无训练 predictor 在 phase_lopo 与 system_loso 上指标相同是预期的：两者的测试点并集都是全部相对的全温区，权重比例一致。

temp_extrap 的无训练 predictor 只在共享的 `T>T*` 测试点评测；zero 排除 crist–trid 后为 3.4692 meV/atom（全网格为 5.4372），按测试点评测后与参考口径一致。

`_legacy/` contains six cleaned historical delivery records and `CLEANING.md`.
They are tracked for provenance but not read by table scripts.
