# Residual Oracle Diagnostics

This report uses only canonical processed reference data, the frozen v2 temperature split, provenance tables, and canonical QH outputs. No model was trained or evaluated by this experiment.

## Input consistency

All 9 in-scope phases have matching checkpoint SHA-256 and head provenance between E_DPA, F_QH, and the head policy: **True**.

The QH CSV schema and `src/lib/fes_bench/qh/run.py` confirm that `F_QH` already includes `E_static + F_vib`; the baselines therefore use F_QH directly and do not add E_DPA again. The QH thermal properties are the quantum phonopy free energy and include zero-point energy. The QH runner records imaginary-mode diagnostics with a -0.05 THz cutoff; this script does not edit or recalculate those modes.

## Delta-E budget

| Pair (low/high) | Delta E DPA (meV/atom) | Delta E effective median [range] (meV/atom) | DPA minus effective | Reference Delta G SD | |error| / SD |
| --- | ---: | ---: | ---: | ---: | ---: |
| hf bcc/hcp | -152.150 | 98.839 [64.196, 132.365] | -250.988 | 11.987 | 20.938 |
| ti bcc/hcp | -93.749 | 64.237 [30.909, 96.317] | -157.986 | 6.111 | 25.851 |
| zr bcc/hcp | -65.487 | 59.543 [32.904, 82.505] | -125.030 | 8.581 | 14.571 |
| sio2 quartz_beta/cristobalite_beta | -9.996 | -9.955 [-11.590, -8.003] | -0.041 | 3.607 | 0.011 |
| sio2 quartz_beta/tridymite_p63mmc | -7.361 | -9.596 [-11.731, -7.905] | 2.235 | 3.576 | 0.625 |

The effective Delta E range is the six extrapolations from F1-F3 and the full/lowest-300-K windows. Values are an effective classical reference extrapolation, not a zero-point-inclusive ground-state energy.

## R3 crossing errors in the test region

| Pair | Reference Tc (K) | B0 error (K) | B1 error (K) | B2 error (K) | B0/B1/B2 test MAE (meV/atom) |
| --- | ---: | ---: | ---: | ---: | --- |
| hf bcc/hcp | 1919.000 | -6.108 | -6.132 | -6.071 | 0.526 / 0.528 / 0.523 |
| ti bcc/hcp | 1012.000 | -3.195 | -3.217 | -3.176 | 1.033 / 0.961 / 0.927 |
| zr bcc/hcp | 975.500 | -1.905 | -1.963 | -1.804 | 0.692 / 0.711 / 0.660 |
| sio2 quartz_beta/cristobalite_beta | 1542.709 | -1.800 | -1.842 | -14.487 | 0.149 / 0.151 / 1.347 |
| sio2 quartz_beta/tridymite_p63mmc | 1582.322 | -0.972 | -1.028 | -18.735 | 0.047 / 0.049 / 1.245 |

The reference-curve full-range R3 control RMSE is:
- hf bcc/hcp: 0.089 meV/atom
- ti bcc/hcp: 0.057 meV/atom
- zr bcc/hcp: 0.059 meV/atom
- sio2 quartz_beta/cristobalite_beta: 0.011 meV/atom
- sio2 quartz_beta/tridymite_p63mmc: 0.002 meV/atom

## Conclusions

- **hf bcc/hcp**: the absolute Delta-E discrepancy is larger than the reference Delta G standard deviation (20.938 times the standard deviation). Static energy difference is the main error source.
- **ti bcc/hcp**: the absolute Delta-E discrepancy is larger than the reference Delta G standard deviation (25.851 times the standard deviation). Static energy difference is the main error source.
- **zr bcc/hcp**: the absolute Delta-E discrepancy is larger than the reference Delta G standard deviation (14.571 times the standard deviation). Static energy difference is the main error source.
- **sio2 quartz_beta/cristobalite_beta**: the absolute Delta-E discrepancy is smaller or comparable than the reference Delta G standard deviation (0.011 times the standard deviation).
- **sio2 quartz_beta/tridymite_p63mmc**: the absolute Delta-E discrepancy is smaller or comparable than the reference Delta G standard deviation (0.625 times the standard deviation).

- **hf bcc/hcp**, equal three-parameter residual form: B0 |Tc error| = 6.108 K; B1 = 6.132 K; B2 = 6.071 K. Baseline improvement over B0: B2.
- **ti bcc/hcp**, equal three-parameter residual form: B0 |Tc error| = 3.195 K; B1 = 3.217 K; B2 = 3.176 K. Baseline improvement over B0: B2.
- **zr bcc/hcp**, equal three-parameter residual form: B0 |Tc error| = 1.905 K; B1 = 1.963 K; B2 = 1.804 K. Baseline improvement over B0: B2.
- **sio2 quartz_beta/cristobalite_beta**, equal three-parameter residual form: B0 |Tc error| = 1.800 K; B1 = 1.842 K; B2 = 14.487 K. Baseline improvement over B0: neither B1 nor B2.
- **sio2 quartz_beta/tridymite_p63mmc**, equal three-parameter residual form: B0 |Tc error| = 0.972 K; B1 = 1.028 K; B2 = 18.735 K. Baseline improvement over B0: neither B1 nor B2.

- **hf bcc/hcp**: compared with B0, B1 RMS ratio 24.104, B2 RMS ratio 55.632; bcc qh_reliable=False, hcp qh_reliable=True.
- **ti bcc/hcp**: compared with B0, B1 RMS ratio 20.684, B2 RMS ratio 84.528; bcc qh_reliable=False, hcp qh_reliable=True.
- **zr bcc/hcp**: compared with B0, B1 RMS ratio 14.406, B2 RMS ratio 54.670; bcc qh_reliable=False, hcp qh_reliable=True.
- **sio2 quartz_beta/cristobalite_beta**: compared with B0, B1 RMS ratio 6.558, B2 RMS ratio 302.957; quartz_beta qh_reliable=True, cristobalite_beta qh_reliable=False.
- **sio2 quartz_beta/tridymite_p63mmc**: compared with B0, B1 RMS ratio 6.661, B2 RMS ratio 310.773; quartz_beta qh_reliable=True, tridymite_p63mmc qh_reliable=False.

## Scope and deviations

- The cristobalite-beta/tridymite-beta pair has no reference crossing and was skipped.
- The optional classical-QH variants B1c/B2c were not run because the canonical QH artifacts preserve scalar F_QH outputs and minimum-frequency diagnostics, not the full positive phonon-frequency spectrum needed for a fresh classical vibrational sum.
- Canonical QH raw-run paths under `result/experiments/` were used because the current processed tree does not materialize every phase QH CSV; no legacy, archive, or invalid source was read.
