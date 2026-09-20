# Phase 3 delivery — frozen splits and evaluator

Status: frozen reference and noise-regression evaluator path is validated on
V100 for all three split kinds. Predictor implementations beyond the two
acceptance predictors remain a later Phase 4–5 task.

## Files changed or added

- `fes_bench/splits/make.py`: deterministic temperature extrapolation, phase LOPO, and system LOSO JSON splits.
- `fes_bench/eval/run.py`: reference and smooth Gaussian-noise predictors,
  Delta-G/Tc metrics, crossing-local ΔG MAE for slope conversion,
  quantized-zero handling, and fold aggregation.
- `fes_bench/eval/plot.py`: Matplotlib Delta-G comparison rendering.
- `configs/splits.yaml`, `splits/{temp_extrap,phase_lopo,system_loso}.json`, and `tests/test_eval_roundtrip.py`.

## Actual V100 commands and logs

The r5 validation ran in `/share/jzr/fes-bench/runs/phase3_tests_20260919_r5/fes-bench`.

`PYTHONPATH=$PWD /root/fes-bench-env/bin/python -m pytest -q`

It returned `5 passed in 2.08s`. A reference-passthrough evaluator run then generated metrics for all three split kinds. The command log is `/share/jzr/fes-bench/runs/phase3_tests_20260919_r5/phase3_tests.log`.

The r4 V100 run also evaluates `reference_noise:0.005`, under `/share/jzr/fes-bench/runs/phase3_eval_20260919_r4/fes-bench/results/`.

After correcting LOPO's missing pair support, r12 was run in
`/share/jzr/fes-bench/runs/phase3_eval_20260919_r12/fes-bench` with the full
test suite plus `reference` evaluation and curve rendering for each of
`temp_extrap`, `phase_lopo`, and `system_loso`. It returned `12 passed in
4.26s`. The reference predictor has G MAE `0.0` for every split; its
materialized pair-curve counts are 4 (temperature extrapolation), 8 (phase
LOPO), and 4 (system LOSO).

The table-format enhancement was then isolated in r13 at
`/share/jzr/fes-bench/runs/phase3_eval_20260919_r13/fes-bench`. Its full suite
returned `13 passed in 4.23s`; all three reference evaluations completed, and
the LOPO summary was checked to contain both `Tc from ΔG` and
`test_and_train_partner` fields.

The frozen-split verifier was run in r14 at
`/share/jzr/fes-bench/runs/phase3_verify_20260919_r14/fes-bench`. It returned
`14 passed in 5.00s` and wrote a passing integrity report. The verifier
recomputes every embedded split SHA-256 and checks unique, disjoint
`(system, phase, T_index)` train/test triples in each materialized fold.

## Artifacts and acceptance evidence

- `splits/temp_extrap.json`: 3264 train and 4187 test frames; available crossings are at least 100 K into test.
- `splits/phase_lopo.json`: five materialized folds; `splits/system_loso.json`: two folds.
- `results/reference/temp_extrap/metrics.json`: G MAE, all Delta-G errors, all Tc errors, and false/missed crossing counts are zero.
- `results/reference/{phase_lopo,system_loso}/metrics.json`: reference
  passthrough metrics for frozen folds. LOPO now includes held-phase versus
  train-phase partner Delta-G/Tc curves, explicitly marked
  `pair_support: test_and_train_partner`; scalar G error remains test-only.
- `results/reference/{temp_extrap,phase_lopo,system_loso}/delta_g_curves.{png,pdf}`:
  labelled Hf and SiO2 Delta-G curves in meV/atom, including every LOPO fold.
- `results/reference/*/summary.md`: Table-5/6-style scalar and pairwise
  Markdown summaries with Delta-G errors, signs, direct and slope-converted
  Tc errors, and false/missed crossings.
- `results/phase3/split_integrity.json`: passing frozen-split hash and
  train/test-disjointness audit.
- `results/phase3_thu_hf_revalidate_20260920.log`: additional Hf-only
  revalidation on the dedicated thu-GenSi host; all three reference
  passthrough split paths returned exact zero error. The full two-system
  evidence remains the V100 r12--r14 run set.
- `results/phase3/v100_r12_r14.log`: local copy of the decisive V100 command
  excerpts and test outputs.
- `results/phase3/v100_r18_noise_20260921.log`: prior smooth-gauge 5 meV
  evaluation and all three rendered curve outputs. The fixture-level 20%
  slope check passes; active quantized tables show a worst-case 1.824 ratio,
  which is retained as a documented domain limitation. Per-pair crossing
  slopes, direct errors, and converted errors are tabulated in that log. The
  subsequent crossing-local-MAE change passed 14 full tests in the same
  isolated run.
- `results/phase3/v100_r19_constant_root_20260921.log`: constant-offset
  slope-control rerun and the renamed iid root-stability test. The latter
  reports Tc scatter without a pass threshold.
- `tests/test_eval_roundtrip.py`: exact passthrough and 20% Tc slope-conversion checks on a linear two-phase fixture.

## Deviations from design

- Splits cover Hf and three supported SiO2 phases only. Ti/Zr have no representative structures and CaSiO3 has no absolute per-phase G(T,P); neither is silently included.
- LOPO and LOSO are named folds in one frozen JSON, expanded and aggregated by the evaluator rather than emitted as one file per fold.
- The independent repository is now initialized.  The frozen splits and
  integrity report record the real provenance commit
  `b5b0b2e679b2da66a593d713afc677f1adeab61e`; embedded SHA-256 remains the
  content-level freeze evidence.
- LOPO pair curves necessarily use the in-domain training phase as the
  comparison partner. That partner is documented in each pair metric and is
  excluded from test-only scalar G metrics; without it a one-phase holdout
  cannot define Delta-G or Tc.
- The lightweight V100 image has Matplotlib but no seaborn. The plotter uses
  the equivalent built-in `seaborn-v0_8-whitegrid` style and explicit DejaVu
  Sans/panel styling when seaborn is unavailable, rather than adding a runtime
  dependency to the benchmark image.
- Hf's source has a quantized zero plateau at Tc. The evaluator merges only adjacent roots caused by that plateau. The 20% noise check uses a smooth linear fixture; Hf quantization remains visible in provenance/metrics.
- The fresh full-domain noise run confirms the limitation quantitatively:
  no false/missed crossings, but the Hf plateau can make direct Tc error differ
  from the local slope estimate by more than 20%. This is reported rather than
  altering the source table or relaxing the acceptance fixture.
