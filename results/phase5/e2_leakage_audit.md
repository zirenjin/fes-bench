# E2 — SiO2 calibration-window leakage audit

Calibration is one least-squares energy-gauge constant for the SiO2 predictor: `mean(G_reference − G_checkpoint)`. No model is retrained. Formal constants use 851–1442 K (canonical temp-extrap train indexes 0–591); audit constants improperly use 1443–2499 K (held-out indexes 592–1648).

The frozen historic checkpoint configurations have no `c_system` field. The two constants below are therefore post-hoc predictor-layer audit variants, not a claim about the original training protocol.

| predictor/seed | c_train (eV/atom) | c_heldout (eV/atom) | held-out G MAE train-c | held-out G MAE heldout-c | ΔMAE | relative ΔMAE | ΔG MAE change | Tc change |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| polynomial/seed11 | -3.1468 | 2.98839 | 6.13519 | 2.94756 | -3.18762 | -51.956% | 2.82e-17 | 2.36e-11 |
| polynomial/seed23 | -3.1157 | 3.0639 | 6.1796 | 2.96648 | -3.21312 | -51.996% | -2.95e-17 | 1.73e-10 |
| polynomial/seed37 | -3.12066 | 3.04307 | 6.16372 | 2.95949 | -3.20423 | -51.985% | 1.71e-17 | 1.46e-10 |
| polynomial/seed51 | -3.14559 | 2.99699 | 6.14258 | 2.95135 | -3.19123 | -51.953% | -4.45e-17 | 1.36e-10 |
| polynomial/seed67 | -3.14611 | 2.9864 | 6.13251 | 2.9464 | -3.1861 | -51.954% | 2.59e-17 | 1.91e-10 |
| tlog_polynomial/seed11 | -3.60731 | -1.62592 | 1.98139 | 0.80151 | -1.17988 | -59.548% | -8.67e-18 | 1.41e-11 |
| tlog_polynomial/seed23 | -3.60653 | -1.62234 | 1.98419 | 0.802554 | -1.18164 | -59.553% | -9.97e-18 | 2.27e-11 |
| tlog_polynomial/seed37 | -3.58687 | -1.57852 | 2.00836 | 0.811624 | -1.19673 | -59.588% | 1.99e-17 | 7.28e-12 |
| tlog_polynomial/seed51 | -3.59381 | -1.5933 | 2.00051 | 0.808476 | -1.19203 | -59.586% | -2.66e-17 | 3.41e-11 |
| tlog_polynomial/seed67 | -3.59366 | -1.59236 | 2.0013 | 0.808669 | -1.19263 | -59.593% | 5.78e-19 | 7.91e-11 |

A common system constant cancels exactly from every ΔG curve. The reported zero ΔG/Tc changes are an expected invariance that the two evaluator runs verify, while the scalar G MAE is sensitive to leakage.
