# Overnight benchmark status — 2026-10-02

## Numeric summary

- Residual oracle covered 5 crossing pairs; the recorded static-energy and residual diagnostics are in `result/experiments/residual_oracle/`.
- T3 completed the planned temp-extrapolation runs; the train-region pair gate retains tlog and excludes polynomial from the main comparison tables.
- The thermo-form floor remains the reference floor for temp extrapolation; T3 held-out metrics and exclusion reasons are recorded in the unified tables.
- T2 QH-residual training has three-seed train/test evidence with the held-out degradation retained in its table and report.
- DFT static assembly and DPA snapshot inference are partial and explicitly provenance-labeled; no missing output is treated as a successful run.

## Line status

### A — residual oracle and diagnostics

Complete and pushed in `c83b62c` (`Add thermo-form floor and correct energy-direction diagnostics`). The metal direction now follows the reference crossing slope: hcp is low and bcc is high. E_DPA provenance, checkpoint SHA-256, head, fitted crossings, and “no meaningful difference” language are recorded in [residual_oracle/summary.md](experiments/residual_oracle/summary.md).

The thermo-form floor uses `a + bT + cTlog(T)` and is trained only on eligible train rows. Phase-LOPO and system-LOSO correctly report `n/a:no_training_phase`.

### B — T3 representation regressors

All planned raw runs are present, and the remote script hash matches the local script:

```text
75e558a5dcb418b2b318f11dcdb539b9befa78bb2797e5bd5c6472384fa84fb3
```

The unified tables and T3 acceptance evidence are in [training_acceptance.csv](experiments/t3_temp_extrap/training_acceptance.csv), [vs_thermo_form.csv](experiments/t3_temp_extrap_v2/vs_thermo_form.csv), and `result/tables/predictor_comparison_*`. The train-region pair gate excludes polynomial and retains tlog; the exclusion reasons are in the temp-extrap table metadata. The table-gating correction was pushed in `20c084b`.

### C — T2 QH residual

The three-seed controlled ablation and `t2_temp_extrap.yaml` are present with provenance and unified tables. The result is committed in `64cddb0` (`Complete T2 QH residual controlled ablation`). The train pair error is below 20 meV/atom, but the held-out test pair error remains 35.748 meV/atom.

### D — DFT static-energy snapshots

Partial. Local VASP outputs exist under `/tmp/fes-bench-dft-static/results` and are assembled into `result/experiments/dft_static/`. Existing Expanse outputs supplied the comparison and convergence evidence; the final-adopted-converged outputs are absent, so the fallback final energies are explicitly recorded. Hf/SiO2 snapshot inference completed in an independent thu-GenSi directory; Ti/Zr inference remains unavailable there. The DFT artifacts were pushed in `2afe4ae`:

```text
src/experiments/dft_static/
result/experiments/dft_static/
```

No remote checkout was modified and no remote Git command was run.

### E — DeltaAI

Not run: no accessible DeltaAI endpoint or credentials were available in this workspace.

### F — integrity

`python3 src/notebooks/run_checks.py` completes with zero warnings. The checks cover both README files, canonical t1 pair names, predictor/t1 pair-set agreement, crossing slopes, and Figure 1 slope consistency. `result/STATUS.md` lists E1–E14 and all out-of-list overnight lines. The DFT snapshot CSV records 721 valid Hf snapshots per phase, 239/240/239 valid SiO2 snapshots, and explicit unavailable rows for Ti/Zr.

## Reproduction commands

```bash
python3 src/notebooks/run_checks.py
PYTHONPATH=src/lib python3 src/tables/predictor_comparison.py --split temp_extrap
PYTHONPATH=src/lib python3 src/tables/predictor_comparison.py --split phase_lopo
PYTHONPATH=src/lib python3 src/tables/predictor_comparison.py --split system_loso
```
