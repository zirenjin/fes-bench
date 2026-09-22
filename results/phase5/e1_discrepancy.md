# E1 discrepancy and reproducibility record

## Hf: early result is not reproducible from a checkpoint

The only located Hf artifact is
`/GenSIvePFS/users/zirenj/fes_experiment_pbe_d3bj/reference/hf_hcp_bcc_pbe/hf_tlog_fewshot_results.json`.
It is a 20-entry metrics-only JSON for label counts 1, 2, 4, and 8 (five
seeds each).  Its leaf fields are only `status`, label indices/temperatures,
held-out scalar metrics, and reference/predicted crossing metrics.  It has no
model path, checkpoint, weights, state dict, or prediction curve.

The often quoted Hf ``≈47 K`` value is traceable to the legacy **8-label**
entry: its five stored crossing errors have mean 46.467 K and range
40.212–53.432 K (mean predicted Tc 1872.533 K against reference 1919 K).
It is an early, metrics-only result and cannot be rerun by `fes_bench.eval.run`.
It is therefore excluded from the reproducible E1 table and must not enter
the paper as a re-evaluated result.  No replacement training was performed.

## SiO2

The cited ``0.36 ± 0.34`` is also now localized: it is the
`ranking_mean=0.358035` and `ranking_std=0.338517` in the legacy
`continuous_sio2/piecewise/summary_corrected.json`.  The `piecewise/`
directory contains metrics JSON only—no `model.pth`—so it is not one of the
frozen checkpoints that can be re-evaluated.  It must not be compared as
though it were a result from the polynomial/tlog TorchScript checkpoints.

The current E1 rerun covers every loadable checkpoint: five `polynomial` and
five `tlog_polynomial` seeds.  On the three shared pairs and all ten seeds,
the new unified evaluator agrees with the adjacent legacy per-checkpoint
metrics: maximum absolute ΔG-MAE difference is `7.44e-08 eV/atom` (mean
`2.50e-08 eV/atom`), and maximum absolute predicted-Tc difference is
`0.0214 K` (mean `0.00531 K`).  The small Tc difference is numerical
root-interpolation precision; the energy curves and pair MAEs agree.  The
new report reverses the first pair's subtraction convention to its sorted
`cristobalite_minus_quartz` key, which changes ΔG sign but not Tc or MAE.
