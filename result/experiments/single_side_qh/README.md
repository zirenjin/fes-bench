# Single-sided QH diagnostic

This diagnostic compares anchored thermal shapes on the QH-reliable side only. It does not evaluate phase pairs, ΔG, ranking, or T_c.

## Results

| Phase | Absolute-G RMSE | One-sided QH RMSE | R = QH / absolute | Improvement | Absolute slope | QH slope |
|---|---:|---:|---:|---:|---:|---:|
| Hf hcp | 421.351 ± 0.671 | 28.069 ± 0.001 | 0.067 ± 0.000 | +93.3% | -1.53214 | 0.10651 |
| Ti hcp | 214.456 ± 0.131 | 52.499 ± 0.000 | 0.245 ± 0.000 | +75.5% | -0.69241 | 0.15826 |
| Zr hcp | 118.946 ± 0.144 | 58.413 ± 0.000 | 0.491 ± 0.001 | +50.9% | -0.41045 | 0.16830 |
| SiO2 beta-quartz | 1382.801 ± 0.192 | 59.610 ± 0.003 | 0.043 ± 0.000 | +95.7% | -2.48975 | 0.10201 |

Each model has one train row and one test row per phase and seed in `per_phase_seed_metrics.csv`; the summary uses the common test seeds only.

The anchor is the final temperature in the frozen temp_extrap v2 training window for each phase. All curves are converted to anchored shape relative to that temperature before the primary RMSE is computed.

QH is quantum phonopy while the canonical DaRUS/TI reference is classical MD. Therefore the anchored thermal-shape error, not raw absolute-G MAE, is the primary comparison.

## Interpretation

- **Hf hcp:** QH reduces the anchored thermal residual (R = 0.067); the three seeds are consistent (R std = 0.0001). The improvement is mainly slope only (mean absolute offset: 8526.3 → 14170.4 meV; slope magnitude: 1.53214 → 0.10651 meV/atom/K).
- **Ti hcp:** QH reduces the anchored thermal residual (R = 0.245); the three seeds are consistent (R std = 0.0001). The improvement is mainly slope only (mean absolute offset: 7524.1 → 13285.1 meV; slope magnitude: 0.69241 → 0.15826 meV/atom/K).
- **Zr hcp:** QH reduces the anchored thermal residual (R = 0.491); the three seeds are consistent (R std = 0.0006). The improvement is mainly slope only (mean absolute offset: 8094.4 → 14981.5 meV; slope magnitude: 0.41045 → 0.16830 meV/atom/K).
- **SiO2 beta-quartz:** QH reduces the anchored thermal residual (R = 0.043); the three seeds are consistent (R std = 0.0000). The improvement is mainly offset and slope (mean absolute offset: 1830.3 → 1149.8 meV; slope magnitude: 2.48975 → 0.10201 meV/atom/K).

Across all four QH-reliable phases, QH gives R < 1, but this is a single-phase thermal-shape result only. It does not establish pairwise phase-stability improvement, does not validate half-QH on a dynamically unstable high-temperature partner, and does not identify a true entropy error.

## Reproduce

```bash
python src/experiments/single_side_qh.py --repo-root .
```

The script verifies QH reliability, the frozen split hash, exact temperature grids, common seed sets, eV/atom input units, and zero anchored residual at T_ref. Prediction CSVs are the canonical T3 tlog and T2 QH-residual outputs exported from the isolated runs; their source paths, checkpoint placeholders, heads, hashes, and commits are recorded in `per_phase_seed_metrics.csv` and `findings.meta.json`.
