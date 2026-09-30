# Zero ΔG floor

Skill scores in the tables use the zero floor uniformly: `1 − MAE / MAE_zero`, with `MAE_zero` equal to the mean reference |ΔG| on evaluated points. `constant_delta_g` is fitted separately for each ordered ΔG = G(left) − G(right) from pair *training* points and is unavailable when a fold has no pairwise training labels.

Phase-level G MAE and coverage are not defined for these pair predictors.

## temp_extrap

### all

Fittable pairs: 0; aggregate ΔG MAE: — eV/atom; aggregate sign accuracy: —.

| Pair | status | train mean ΔG (eV/atom) | ΔG MAE | ΔG RMSE | sign accuracy | reference / predicted Tc (K) | false / missed |
|---|---|---:|---:|---:|---:|---|---:|
| hf:bcc_minus_hcp | degenerate_prediction | 0 | 0.00423642 | 0.00507363 | 0.00590551 | 1919 / n/a:degenerate_prediction | n/a:degenerate_prediction / n/a:degenerate_prediction |
| sio2:cristobalite_beta_minus_quartz_beta | degenerate_prediction | 0 | 0.00339517 | 0.00409094 | 0 | 1542.71 / n/a:degenerate_prediction | n/a:degenerate_prediction / n/a:degenerate_prediction |
| sio2:cristobalite_beta_minus_tridymite_p63mmc | degenerate_prediction | 0 | 0.000285824 | 0.000287805 | 0 | — / n/a:degenerate_prediction | n/a:degenerate_prediction / n/a:degenerate_prediction |
| sio2:quartz_beta_minus_tridymite_p63mmc | degenerate_prediction | 0 | 0.0031745 | 0.00387741 | 0 | 1582.32 / n/a:degenerate_prediction | n/a:degenerate_prediction / n/a:degenerate_prediction |
| ti:bcc_minus_hcp | degenerate_prediction | 0 | 0.00680569 | 0.00797444 | 0.00487805 | 1012 / n/a:degenerate_prediction | n/a:degenerate_prediction / n/a:degenerate_prediction |
| zr:bcc_minus_hcp | degenerate_prediction | 0 | 0.00661941 | 0.00773835 | 0.00310559 | 975.5 / n/a:degenerate_prediction | n/a:degenerate_prediction / n/a:degenerate_prediction |

## phase_lopo

### hf:bcc

Fittable pairs: 0; aggregate ΔG MAE: — eV/atom; aggregate sign accuracy: —.

| Pair | status | train mean ΔG (eV/atom) | ΔG MAE | ΔG RMSE | sign accuracy | reference / predicted Tc (K) | false / missed |
|---|---|---:|---:|---:|---:|---|---:|
| hf:bcc_minus_hcp | degenerate_prediction | 0 | 0.0113998 | 0.0145146 | 0.00239617 | 1919 / n/a:degenerate_prediction | n/a:degenerate_prediction / n/a:degenerate_prediction |

### hf:hcp

Fittable pairs: 0; aggregate ΔG MAE: — eV/atom; aggregate sign accuracy: —.

| Pair | status | train mean ΔG (eV/atom) | ΔG MAE | ΔG RMSE | sign accuracy | reference / predicted Tc (K) | false / missed |
|---|---|---:|---:|---:|---:|---|---:|
| hf:bcc_minus_hcp | degenerate_prediction | 0 | 0.0113998 | 0.0145146 | 0.00239617 | 1919 / n/a:degenerate_prediction | n/a:degenerate_prediction / n/a:degenerate_prediction |

### sio2:cristobalite_beta

Fittable pairs: 0; aggregate ΔG MAE: — eV/atom; aggregate sign accuracy: —.

| Pair | status | train mean ΔG (eV/atom) | ΔG MAE | ΔG RMSE | sign accuracy | reference / predicted Tc (K) | false / missed |
|---|---|---:|---:|---:|---:|---|---:|
| sio2:cristobalite_beta_minus_quartz_beta | degenerate_prediction | 0 | 0.00320909 | 0.0037702 | 0 | 1542.71 / n/a:degenerate_prediction | n/a:degenerate_prediction / n/a:degenerate_prediction |
| sio2:cristobalite_beta_minus_tridymite_p63mmc | degenerate_prediction | 0 | 0.00025389 | 0.000260579 | 0 | — / n/a:degenerate_prediction | n/a:degenerate_prediction / n/a:degenerate_prediction |

### sio2:quartz_beta

Fittable pairs: 0; aggregate ΔG MAE: — eV/atom; aggregate sign accuracy: —.

| Pair | status | train mean ΔG (eV/atom) | ΔG MAE | ΔG RMSE | sign accuracy | reference / predicted Tc (K) | false / missed |
|---|---|---:|---:|---:|---:|---|---:|
| sio2:cristobalite_beta_minus_quartz_beta | degenerate_prediction | 0 | 0.00320909 | 0.0037702 | 0 | 1542.71 / n/a:degenerate_prediction | n/a:degenerate_prediction / n/a:degenerate_prediction |
| sio2:quartz_beta_minus_tridymite_p63mmc | degenerate_prediction | 0 | 0.00313832 | 0.00367457 | 0 | 1582.32 / n/a:degenerate_prediction | n/a:degenerate_prediction / n/a:degenerate_prediction |

### sio2:tridymite_p63mmc

Fittable pairs: 0; aggregate ΔG MAE: — eV/atom; aggregate sign accuracy: —.

| Pair | status | train mean ΔG (eV/atom) | ΔG MAE | ΔG RMSE | sign accuracy | reference / predicted Tc (K) | false / missed |
|---|---|---:|---:|---:|---:|---|---:|
| sio2:cristobalite_beta_minus_tridymite_p63mmc | degenerate_prediction | 0 | 0.00025389 | 0.000260579 | 0 | — / n/a:degenerate_prediction | n/a:degenerate_prediction / n/a:degenerate_prediction |
| sio2:quartz_beta_minus_tridymite_p63mmc | degenerate_prediction | 0 | 0.00313832 | 0.00367457 | 0 | 1582.32 / n/a:degenerate_prediction | n/a:degenerate_prediction / n/a:degenerate_prediction |

### ti:bcc

Fittable pairs: 0; aggregate ΔG MAE: — eV/atom; aggregate sign accuracy: —.

| Pair | status | train mean ΔG (eV/atom) | ΔG MAE | ΔG RMSE | sign accuracy | reference / predicted Tc (K) | false / missed |
|---|---|---:|---:|---:|---:|---|---:|
| ti:bcc_minus_hcp | degenerate_prediction | 0 | 0.00669683 | 0.00777121 | 0.004329 | 1012 / n/a:degenerate_prediction | n/a:degenerate_prediction / n/a:degenerate_prediction |

### ti:hcp

Fittable pairs: 0; aggregate ΔG MAE: — eV/atom; aggregate sign accuracy: —.

| Pair | status | train mean ΔG (eV/atom) | ΔG MAE | ΔG RMSE | sign accuracy | reference / predicted Tc (K) | false / missed |
|---|---|---:|---:|---:|---:|---|---:|
| ti:bcc_minus_hcp | degenerate_prediction | 0 | 0.00669683 | 0.00777121 | 0.004329 | 1012 / n/a:degenerate_prediction | n/a:degenerate_prediction / n/a:degenerate_prediction |

### zr:bcc

Fittable pairs: 0; aggregate ΔG MAE: — eV/atom; aggregate sign accuracy: —.

| Pair | status | train mean ΔG (eV/atom) | ΔG MAE | ΔG RMSE | sign accuracy | reference / predicted Tc (K) | false / missed |
|---|---|---:|---:|---:|---:|---|---:|
| zr:bcc_minus_hcp | degenerate_prediction | 0 | 0.00759824 | 0.008687 | 0.0021978 | 975.5 / n/a:degenerate_prediction | n/a:degenerate_prediction / n/a:degenerate_prediction |

### zr:hcp

Fittable pairs: 0; aggregate ΔG MAE: — eV/atom; aggregate sign accuracy: —.

| Pair | status | train mean ΔG (eV/atom) | ΔG MAE | ΔG RMSE | sign accuracy | reference / predicted Tc (K) | false / missed |
|---|---|---:|---:|---:|---:|---|---:|
| zr:bcc_minus_hcp | degenerate_prediction | 0 | 0.00759824 | 0.008687 | 0.0021978 | 975.5 / n/a:degenerate_prediction | n/a:degenerate_prediction / n/a:degenerate_prediction |

## system_loso

### hf

Fittable pairs: 0; aggregate ΔG MAE: — eV/atom; aggregate sign accuracy: —.

| Pair | status | train mean ΔG (eV/atom) | ΔG MAE | ΔG RMSE | sign accuracy | reference / predicted Tc (K) | false / missed |
|---|---|---:|---:|---:|---:|---|---:|
| hf:bcc_minus_hcp | degenerate_prediction | 0 | 0.0113998 | 0.0145146 | 0.00239617 | 1919 / n/a:degenerate_prediction | n/a:degenerate_prediction / n/a:degenerate_prediction |

### ti

Fittable pairs: 0; aggregate ΔG MAE: — eV/atom; aggregate sign accuracy: —.

| Pair | status | train mean ΔG (eV/atom) | ΔG MAE | ΔG RMSE | sign accuracy | reference / predicted Tc (K) | false / missed |
|---|---|---:|---:|---:|---:|---|---:|
| ti:bcc_minus_hcp | degenerate_prediction | 0 | 0.00669683 | 0.00777121 | 0.004329 | 1012 / n/a:degenerate_prediction | n/a:degenerate_prediction / n/a:degenerate_prediction |

### zr

Fittable pairs: 0; aggregate ΔG MAE: — eV/atom; aggregate sign accuracy: —.

| Pair | status | train mean ΔG (eV/atom) | ΔG MAE | ΔG RMSE | sign accuracy | reference / predicted Tc (K) | false / missed |
|---|---|---:|---:|---:|---:|---|---:|
| zr:bcc_minus_hcp | degenerate_prediction | 0 | 0.00759824 | 0.008687 | 0.0021978 | 975.5 / n/a:degenerate_prediction | n/a:degenerate_prediction / n/a:degenerate_prediction |

## Deviations from design

No new model was trained. The design evaluator's phase-level G MAE/coverage cannot be applied to a pair-only ΔG predictor without inventing a phase-energy gauge, so they are explicitly N/A. In phase-LOPO and system-LOSO, a held phase/system has no pairwise training ΔG labels; `constant_delta_g` is consequently reported as unavailable rather than fitted with held-out reference values.


# Training-mean constant ΔG floor

Skill scores in the tables use the zero floor uniformly: `1 − MAE / MAE_zero`, with `MAE_zero` equal to the mean reference |ΔG| on evaluated points. `constant_delta_g` is fitted separately for each ordered ΔG = G(left) − G(right) from pair *training* points and is unavailable when a fold has no pairwise training labels.

Phase-level G MAE and coverage are not defined for these pair predictors.

## temp_extrap

### all

Fittable pairs: 6; aggregate ΔG MAE: 0.0100613 eV/atom; aggregate sign accuracy: 0.288841.

| Pair | status | train mean ΔG (eV/atom) | ΔG MAE | ΔG RMSE | sign accuracy | reference / predicted Tc (K) | false / missed |
|---|---|---:|---:|---:|---:|---|---:|
| hf:bcc_minus_hcp | ok | 0.016291 | 0.0199794 | 0.0202809 | 0.192913 | 1919 / — | 0 / 1 |
| sio2:cristobalite_beta_minus_quartz_beta | ok | 0.00287686 | 0.00620058 | 0.00664348 | 0.0946074 | 1542.71 / — | 0 / 1 |
| sio2:cristobalite_beta_minus_tridymite_p63mmc | ok | -0.000196871 | 8.90505e-05 | 9.51248e-05 | 1 | — / — | 0 / 0 |
| sio2:quartz_beta_minus_tridymite_p63mmc | ok | -0.00307373 | 0.00611163 | 0.00656944 | 0.13245 | 1582.32 / — | 0 / 1 |
| ti:bcc_minus_hcp | ok | 0.00583846 | 0.0119879 | 0.0130187 | 0.15935 | 1012 / — | 0 / 1 |
| zr:bcc_minus_hcp | ok | 0.00996805 | 0.0159996 | 0.0167179 | 0.153727 | 975.5 / — | 0 / 1 |

## phase_lopo

### hf:bcc

Fittable pairs: 0; aggregate ΔG MAE: — eV/atom; aggregate sign accuracy: —.

| Pair | status | train mean ΔG (eV/atom) | ΔG MAE | ΔG RMSE | sign accuracy | reference / predicted Tc (K) | false / missed |
|---|---|---:|---:|---:|---:|---|---:|
| hf:bcc_minus_hcp | not_fittable_without_pair_training_labels | — | n/a:no_training_phase | n/a:no_training_phase | n/a:no_training_phase | 1919 / n/a:no_training_phase | n/a:no_training_phase / n/a:no_training_phase |

### hf:hcp

Fittable pairs: 0; aggregate ΔG MAE: — eV/atom; aggregate sign accuracy: —.

| Pair | status | train mean ΔG (eV/atom) | ΔG MAE | ΔG RMSE | sign accuracy | reference / predicted Tc (K) | false / missed |
|---|---|---:|---:|---:|---:|---|---:|
| hf:bcc_minus_hcp | not_fittable_without_pair_training_labels | — | n/a:no_training_phase | n/a:no_training_phase | n/a:no_training_phase | 1919 / n/a:no_training_phase | n/a:no_training_phase / n/a:no_training_phase |

### sio2:cristobalite_beta

Fittable pairs: 0; aggregate ΔG MAE: — eV/atom; aggregate sign accuracy: —.

| Pair | status | train mean ΔG (eV/atom) | ΔG MAE | ΔG RMSE | sign accuracy | reference / predicted Tc (K) | false / missed |
|---|---|---:|---:|---:|---:|---|---:|
| sio2:cristobalite_beta_minus_quartz_beta | not_fittable_without_pair_training_labels | — | n/a:no_training_phase | n/a:no_training_phase | n/a:no_training_phase | 1542.71 / n/a:no_training_phase | n/a:no_training_phase / n/a:no_training_phase |
| sio2:cristobalite_beta_minus_tridymite_p63mmc | not_fittable_without_pair_training_labels | — | n/a:no_training_phase | n/a:no_training_phase | n/a:no_training_phase | — / n/a:no_training_phase | n/a:no_training_phase / n/a:no_training_phase |

### sio2:quartz_beta

Fittable pairs: 0; aggregate ΔG MAE: — eV/atom; aggregate sign accuracy: —.

| Pair | status | train mean ΔG (eV/atom) | ΔG MAE | ΔG RMSE | sign accuracy | reference / predicted Tc (K) | false / missed |
|---|---|---:|---:|---:|---:|---|---:|
| sio2:cristobalite_beta_minus_quartz_beta | not_fittable_without_pair_training_labels | — | n/a:no_training_phase | n/a:no_training_phase | n/a:no_training_phase | 1542.71 / n/a:no_training_phase | n/a:no_training_phase / n/a:no_training_phase |
| sio2:quartz_beta_minus_tridymite_p63mmc | not_fittable_without_pair_training_labels | — | n/a:no_training_phase | n/a:no_training_phase | n/a:no_training_phase | 1582.32 / n/a:no_training_phase | n/a:no_training_phase / n/a:no_training_phase |

### sio2:tridymite_p63mmc

Fittable pairs: 0; aggregate ΔG MAE: — eV/atom; aggregate sign accuracy: —.

| Pair | status | train mean ΔG (eV/atom) | ΔG MAE | ΔG RMSE | sign accuracy | reference / predicted Tc (K) | false / missed |
|---|---|---:|---:|---:|---:|---|---:|
| sio2:cristobalite_beta_minus_tridymite_p63mmc | not_fittable_without_pair_training_labels | — | n/a:no_training_phase | n/a:no_training_phase | n/a:no_training_phase | — / n/a:no_training_phase | n/a:no_training_phase / n/a:no_training_phase |
| sio2:quartz_beta_minus_tridymite_p63mmc | not_fittable_without_pair_training_labels | — | n/a:no_training_phase | n/a:no_training_phase | n/a:no_training_phase | 1582.32 / n/a:no_training_phase | n/a:no_training_phase / n/a:no_training_phase |

### ti:bcc

Fittable pairs: 0; aggregate ΔG MAE: — eV/atom; aggregate sign accuracy: —.

| Pair | status | train mean ΔG (eV/atom) | ΔG MAE | ΔG RMSE | sign accuracy | reference / predicted Tc (K) | false / missed |
|---|---|---:|---:|---:|---:|---|---:|
| ti:bcc_minus_hcp | not_fittable_without_pair_training_labels | — | n/a:no_training_phase | n/a:no_training_phase | n/a:no_training_phase | 1012 / n/a:no_training_phase | n/a:no_training_phase / n/a:no_training_phase |

### ti:hcp

Fittable pairs: 0; aggregate ΔG MAE: — eV/atom; aggregate sign accuracy: —.

| Pair | status | train mean ΔG (eV/atom) | ΔG MAE | ΔG RMSE | sign accuracy | reference / predicted Tc (K) | false / missed |
|---|---|---:|---:|---:|---:|---|---:|
| ti:bcc_minus_hcp | not_fittable_without_pair_training_labels | — | n/a:no_training_phase | n/a:no_training_phase | n/a:no_training_phase | 1012 / n/a:no_training_phase | n/a:no_training_phase / n/a:no_training_phase |

### zr:bcc

Fittable pairs: 0; aggregate ΔG MAE: — eV/atom; aggregate sign accuracy: —.

| Pair | status | train mean ΔG (eV/atom) | ΔG MAE | ΔG RMSE | sign accuracy | reference / predicted Tc (K) | false / missed |
|---|---|---:|---:|---:|---:|---|---:|
| zr:bcc_minus_hcp | not_fittable_without_pair_training_labels | — | n/a:no_training_phase | n/a:no_training_phase | n/a:no_training_phase | 975.5 / n/a:no_training_phase | n/a:no_training_phase / n/a:no_training_phase |

### zr:hcp

Fittable pairs: 0; aggregate ΔG MAE: — eV/atom; aggregate sign accuracy: —.

| Pair | status | train mean ΔG (eV/atom) | ΔG MAE | ΔG RMSE | sign accuracy | reference / predicted Tc (K) | false / missed |
|---|---|---:|---:|---:|---:|---|---:|
| zr:bcc_minus_hcp | not_fittable_without_pair_training_labels | — | n/a:no_training_phase | n/a:no_training_phase | n/a:no_training_phase | 975.5 / n/a:no_training_phase | n/a:no_training_phase / n/a:no_training_phase |

## system_loso

### hf

Fittable pairs: 0; aggregate ΔG MAE: — eV/atom; aggregate sign accuracy: —.

| Pair | status | train mean ΔG (eV/atom) | ΔG MAE | ΔG RMSE | sign accuracy | reference / predicted Tc (K) | false / missed |
|---|---|---:|---:|---:|---:|---|---:|
| hf:bcc_minus_hcp | not_fittable_without_pair_training_labels | — | n/a:no_training_phase | n/a:no_training_phase | n/a:no_training_phase | 1919 / n/a:no_training_phase | n/a:no_training_phase / n/a:no_training_phase |

### ti

Fittable pairs: 0; aggregate ΔG MAE: — eV/atom; aggregate sign accuracy: —.

| Pair | status | train mean ΔG (eV/atom) | ΔG MAE | ΔG RMSE | sign accuracy | reference / predicted Tc (K) | false / missed |
|---|---|---:|---:|---:|---:|---|---:|
| ti:bcc_minus_hcp | not_fittable_without_pair_training_labels | — | n/a:no_training_phase | n/a:no_training_phase | n/a:no_training_phase | 1012 / n/a:no_training_phase | n/a:no_training_phase / n/a:no_training_phase |

### zr

Fittable pairs: 0; aggregate ΔG MAE: — eV/atom; aggregate sign accuracy: —.

| Pair | status | train mean ΔG (eV/atom) | ΔG MAE | ΔG RMSE | sign accuracy | reference / predicted Tc (K) | false / missed |
|---|---|---:|---:|---:|---:|---|---:|
| zr:bcc_minus_hcp | not_fittable_without_pair_training_labels | — | n/a:no_training_phase | n/a:no_training_phase | n/a:no_training_phase | 975.5 / n/a:no_training_phase | n/a:no_training_phase / n/a:no_training_phase |

## Deviations from design

No new model was trained. The design evaluator's phase-level G MAE/coverage cannot be applied to a pair-only ΔG predictor without inventing a phase-energy gauge, so they are explicitly N/A. In phase-LOPO and system-LOSO, a held phase/system has no pairwise training ΔG labels; `constant_delta_g` is consequently reported as unavailable rather than fitted with held-out reference values.

