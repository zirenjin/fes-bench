# Overnight benchmark status — 2026-10-02

## Numeric summary

- Residual oracle covered 5 crossing pairs; metal static-energy errors are 125.030–250.988 meV/atom, while both SiO2 pairs are within 2.235 meV/atom of the effective reference.
- T3 completed 24 runs (4 systems × 2 bases × seeds 11/23/37); polynomial pair-train MAE was 22.026–22.741 meV/atom, and tlog was 5.901–5.995 meV/atom.
- The 20 meV/atom pair gate passes for all three tlog runs and fails for all three polynomial runs; the phase-energy gate fails for all six runs (0.121–0.439 eV/atom train MAE).
- The thermo-form floor reached 0.424817 meV/atom pair MAE, 0.994771 sign accuracy, and 0.993272 ranking accuracy on temp_extrap; T3 was worse on all three comparisons.
- T2 QH-residual training completed for seeds 11/23/37; train pair MAE was 7.289 meV/atom and test pair MAE was 35.748 meV/atom.

## Line status

### A — residual oracle and diagnostics

Complete and pushed in `c83b62c` (`Add thermo-form floor and correct energy-direction diagnostics`). The metal direction now follows the reference crossing slope: hcp is low and bcc is high. E_DPA provenance, checkpoint SHA-256, head, fitted crossings, and “no meaningful difference” language are recorded in [residual_oracle/summary.md](experiments/residual_oracle/summary.md).

The thermo-form floor uses `a + bT + cTlog(T)` and is trained only on eligible train rows. Phase-LOPO and system-LOSO correctly report `n/a:no_training_phase`.

### B — T3 representation regressors

All 24 planned runs are present, and the remote script hash matches the local script:

```text
75e558a5dcb418b2b318f11dcdb539b9befa78bb2797e5bd5c6472384fa84fb3
```

The unified tables and T3 acceptance evidence are in [training_acceptance.csv](experiments/t3_temp_extrap/training_acceptance.csv), [vs_thermo_form.csv](experiments/t3_temp_extrap/vs_thermo_form.csv), and `result/tables/predictor_comparison_*`. T3 evidence was committed and pushed in `2ed8470` (`Record T3 acceptance and thermo-form comparison`). The run is recorded honestly as not accepted because the phase-energy train gate fails for every run; no additional T3 training process remains active.

### C — T2 QH residual

The three-seed controlled ablation and `t2_temp_extrap.yaml` are present with provenance and unified tables. The result is committed in `64cddb0` (`Complete T2 QH residual controlled ablation`). The train pair error is below 20 meV/atom, but the held-out test pair error remains 35.748 meV/atom.

### D — DFT static-energy snapshots

Blocked/incomplete. Local VASP outputs exist under `/tmp/fes-bench-dft-static/results`, but the assembler requests `conv-k777` OUTCARs that are absent. Local snapshot inference also cannot start because the local environment has no `deepmd` module. The authorized untracked DFT files remain untouched and unstaged:

```text
src/experiments/dft_static/
result/experiments/dft_static/
```

No remote checkout was modified and no remote Git command was run.

### E — DeltaAI

Not run: no accessible DeltaAI endpoint or credentials were available in this workspace.

### F — integrity

`python3 src/notebooks/run_checks.py` completes with zero warnings. The checks cover both README files, canonical t1 pair names, predictor/t1 pair-set agreement, crossing slopes, and Figure 1 slope consistency. The local branch is synchronized with `origin/master`; only the explicitly authorized DFT directories remain untracked.

## Reproduction commands

```bash
python3 src/notebooks/run_checks.py
PYTHONPATH=src/lib python3 src/tables/predictor_comparison.py --split temp_extrap
PYTHONPATH=src/lib python3 src/tables/predictor_comparison.py --split phase_lopo
PYTHONPATH=src/lib python3 src/tables/predictor_comparison.py --split system_loso
```
