# Zero ΔG floor

`constant_delta_g` is fitted separately for each ordered ΔG = G(left) − G(right) from the intersection of that pair's *training* temperature points.  It is deliberately unavailable when a fold contains no pairwise training labels (rather than leaking held-out labels).

Phase-level G MAE and coverage are not defined for these pair predictors.

## temp_extrap

### all

Fittable pairs: 4; aggregate ΔG MAE: 0.00450029 eV/atom; aggregate sign accuracy: 0.000599042.

| Pair | status | train mean ΔG (eV/atom) | ΔG MAE | ΔG RMSE | sign accuracy | reference / predicted Tc (K) | false / missed |
|---|---|---:|---:|---:|---:|---|---:|
| hf:bcc_minus_hcp | ok | 0 | 0.0113998 | 0.0145146 | 0.00239617 | 1919 / 1701.5 | 0 / 0 |
| sio2:cristobalite_beta_minus_quartz_beta | ok | 0 | 0.00320909 | 0.0037702 | 0 | 1542.71 / 1675 | 0 / 0 |
| sio2:cristobalite_beta_minus_tridymite_p63mmc | ok | 0 | 0.00025389 | 0.000260579 | 0 | — / 1675 | 1 / 0 |
| sio2:quartz_beta_minus_tridymite_p63mmc | ok | 0 | 0.00313832 | 0.00367457 | 0 | 1582.32 / 1675 | 0 / 0 |

## phase_lopo

### hf:bcc

Fittable pairs: 1; aggregate ΔG MAE: 0.0113998 eV/atom; aggregate sign accuracy: 0.00239617.

| Pair | status | train mean ΔG (eV/atom) | ΔG MAE | ΔG RMSE | sign accuracy | reference / predicted Tc (K) | false / missed |
|---|---|---:|---:|---:|---:|---|---:|
| hf:bcc_minus_hcp | ok | 0 | 0.0113998 | 0.0145146 | 0.00239617 | 1919 / 1701.5 | 0 / 0 |

### hf:hcp

Fittable pairs: 1; aggregate ΔG MAE: 0.0113998 eV/atom; aggregate sign accuracy: 0.00239617.

| Pair | status | train mean ΔG (eV/atom) | ΔG MAE | ΔG RMSE | sign accuracy | reference / predicted Tc (K) | false / missed |
|---|---|---:|---:|---:|---:|---|---:|
| hf:bcc_minus_hcp | ok | 0 | 0.0113998 | 0.0145146 | 0.00239617 | 1919 / 1701.5 | 0 / 0 |

### sio2:cristobalite_beta

Fittable pairs: 2; aggregate ΔG MAE: 0.00173149 eV/atom; aggregate sign accuracy: 0.

| Pair | status | train mean ΔG (eV/atom) | ΔG MAE | ΔG RMSE | sign accuracy | reference / predicted Tc (K) | false / missed |
|---|---|---:|---:|---:|---:|---|---:|
| sio2:cristobalite_beta_minus_quartz_beta | ok | 0 | 0.00320909 | 0.0037702 | 0 | 1542.71 / 1675 | 0 / 0 |
| sio2:cristobalite_beta_minus_tridymite_p63mmc | ok | 0 | 0.00025389 | 0.000260579 | 0 | — / 1675 | 1 / 0 |

### sio2:quartz_beta

Fittable pairs: 2; aggregate ΔG MAE: 0.00317371 eV/atom; aggregate sign accuracy: 0.

| Pair | status | train mean ΔG (eV/atom) | ΔG MAE | ΔG RMSE | sign accuracy | reference / predicted Tc (K) | false / missed |
|---|---|---:|---:|---:|---:|---|---:|
| sio2:cristobalite_beta_minus_quartz_beta | ok | 0 | 0.00320909 | 0.0037702 | 0 | 1542.71 / 1675 | 0 / 0 |
| sio2:quartz_beta_minus_tridymite_p63mmc | ok | 0 | 0.00313832 | 0.00367457 | 0 | 1582.32 / 1675 | 0 / 0 |

### sio2:tridymite_p63mmc

Fittable pairs: 2; aggregate ΔG MAE: 0.00169611 eV/atom; aggregate sign accuracy: 0.

| Pair | status | train mean ΔG (eV/atom) | ΔG MAE | ΔG RMSE | sign accuracy | reference / predicted Tc (K) | false / missed |
|---|---|---:|---:|---:|---:|---|---:|
| sio2:cristobalite_beta_minus_tridymite_p63mmc | ok | 0 | 0.00025389 | 0.000260579 | 0 | — / 1675 | 1 / 0 |
| sio2:quartz_beta_minus_tridymite_p63mmc | ok | 0 | 0.00313832 | 0.00367457 | 0 | 1582.32 / 1675 | 0 / 0 |

## system_loso

### hf

Fittable pairs: 1; aggregate ΔG MAE: 0.0113998 eV/atom; aggregate sign accuracy: 0.00239617.

| Pair | status | train mean ΔG (eV/atom) | ΔG MAE | ΔG RMSE | sign accuracy | reference / predicted Tc (K) | false / missed |
|---|---|---:|---:|---:|---:|---|---:|
| hf:bcc_minus_hcp | ok | 0 | 0.0113998 | 0.0145146 | 0.00239617 | 1919 / 1701.5 | 0 / 0 |

### sio2

Fittable pairs: 3; aggregate ΔG MAE: 0.00220043 eV/atom; aggregate sign accuracy: 0.

| Pair | status | train mean ΔG (eV/atom) | ΔG MAE | ΔG RMSE | sign accuracy | reference / predicted Tc (K) | false / missed |
|---|---|---:|---:|---:|---:|---|---:|
| sio2:cristobalite_beta_minus_quartz_beta | ok | 0 | 0.00320909 | 0.0037702 | 0 | 1542.71 / 1675 | 0 / 0 |
| sio2:cristobalite_beta_minus_tridymite_p63mmc | ok | 0 | 0.00025389 | 0.000260579 | 0 | — / 1675 | 1 / 0 |
| sio2:quartz_beta_minus_tridymite_p63mmc | ok | 0 | 0.00313832 | 0.00367457 | 0 | 1582.32 / 1675 | 0 / 0 |

## Deviations from design

No new model was trained. The design evaluator's phase-level G MAE/coverage cannot be applied to a pair-only ΔG predictor without inventing a phase-energy gauge, so they are explicitly N/A. In phase-LOPO and system-LOSO, a held phase/system has no pairwise training ΔG labels; `constant_delta_g` is consequently reported as unavailable rather than fitted with held-out reference values.


# Training-mean constant ΔG floor

`constant_delta_g` is fitted separately for each ordered ΔG = G(left) − G(right) from the intersection of that pair's *training* temperature points.  It is deliberately unavailable when a fold contains no pairwise training labels (rather than leaking held-out labels).

Phase-level G MAE and coverage are not defined for these pair predictors.

## temp_extrap

### all

Fittable pairs: 4; aggregate ΔG MAE: 0.00528856 eV/atom; aggregate sign accuracy: 0.634019.

| Pair | status | train mean ΔG (eV/atom) | ΔG MAE | ΔG RMSE | sign accuracy | reference / predicted Tc (K) | false / missed |
|---|---|---:|---:|---:|---:|---|---:|
| hf:bcc_minus_hcp | ok | 0.016291 | 0.0124592 | 0.014471 | 0.672524 | 1919 / — | 0 / 1 |
| sio2:cristobalite_beta_minus_quartz_beta | ok | 0.00287686 | 0.00434745 | 0.00536716 | 0.419648 | 1542.71 / — | 0 / 1 |
| sio2:cristobalite_beta_minus_tridymite_p63mmc | ok | -0.000196871 | 7.26838e-05 | 8.18087e-05 | 1 | — / — | 0 / 0 |
| sio2:quartz_beta_minus_tridymite_p63mmc | ok | -0.00307373 | 0.00427486 | 0.00530447 | 0.443905 | 1582.32 / — | 0 / 1 |

## phase_lopo

### hf:bcc

Fittable pairs: 0; aggregate ΔG MAE: — eV/atom; aggregate sign accuracy: —.

| Pair | status | train mean ΔG (eV/atom) | ΔG MAE | ΔG RMSE | sign accuracy | reference / predicted Tc (K) | false / missed |
|---|---|---:|---:|---:|---:|---|---:|
| hf:bcc_minus_hcp | not_fittable_without_pair_training_labels | — | 0.0113998 | 0.0145146 | 0.00239617 | 1919 / 1701.5 | 0 / 0 |

### hf:hcp

Fittable pairs: 0; aggregate ΔG MAE: — eV/atom; aggregate sign accuracy: —.

| Pair | status | train mean ΔG (eV/atom) | ΔG MAE | ΔG RMSE | sign accuracy | reference / predicted Tc (K) | false / missed |
|---|---|---:|---:|---:|---:|---|---:|
| hf:bcc_minus_hcp | not_fittable_without_pair_training_labels | — | 0.0113998 | 0.0145146 | 0.00239617 | 1919 / 1701.5 | 0 / 0 |

### sio2:cristobalite_beta

Fittable pairs: 0; aggregate ΔG MAE: — eV/atom; aggregate sign accuracy: —.

| Pair | status | train mean ΔG (eV/atom) | ΔG MAE | ΔG RMSE | sign accuracy | reference / predicted Tc (K) | false / missed |
|---|---|---:|---:|---:|---:|---|---:|
| sio2:cristobalite_beta_minus_quartz_beta | not_fittable_without_pair_training_labels | — | 0.00320909 | 0.0037702 | 0 | 1542.71 / 1675 | 0 / 0 |
| sio2:cristobalite_beta_minus_tridymite_p63mmc | not_fittable_without_pair_training_labels | — | 0.00025389 | 0.000260579 | 0 | — / 1675 | 1 / 0 |

### sio2:quartz_beta

Fittable pairs: 0; aggregate ΔG MAE: — eV/atom; aggregate sign accuracy: —.

| Pair | status | train mean ΔG (eV/atom) | ΔG MAE | ΔG RMSE | sign accuracy | reference / predicted Tc (K) | false / missed |
|---|---|---:|---:|---:|---:|---|---:|
| sio2:cristobalite_beta_minus_quartz_beta | not_fittable_without_pair_training_labels | — | 0.00320909 | 0.0037702 | 0 | 1542.71 / 1675 | 0 / 0 |
| sio2:quartz_beta_minus_tridymite_p63mmc | not_fittable_without_pair_training_labels | — | 0.00313832 | 0.00367457 | 0 | 1582.32 / 1675 | 0 / 0 |

### sio2:tridymite_p63mmc

Fittable pairs: 0; aggregate ΔG MAE: — eV/atom; aggregate sign accuracy: —.

| Pair | status | train mean ΔG (eV/atom) | ΔG MAE | ΔG RMSE | sign accuracy | reference / predicted Tc (K) | false / missed |
|---|---|---:|---:|---:|---:|---|---:|
| sio2:cristobalite_beta_minus_tridymite_p63mmc | not_fittable_without_pair_training_labels | — | 0.00025389 | 0.000260579 | 0 | — / 1675 | 1 / 0 |
| sio2:quartz_beta_minus_tridymite_p63mmc | not_fittable_without_pair_training_labels | — | 0.00313832 | 0.00367457 | 0 | 1582.32 / 1675 | 0 / 0 |

## system_loso

### hf

Fittable pairs: 0; aggregate ΔG MAE: — eV/atom; aggregate sign accuracy: —.

| Pair | status | train mean ΔG (eV/atom) | ΔG MAE | ΔG RMSE | sign accuracy | reference / predicted Tc (K) | false / missed |
|---|---|---:|---:|---:|---:|---|---:|
| hf:bcc_minus_hcp | not_fittable_without_pair_training_labels | — | 0.0113998 | 0.0145146 | 0.00239617 | 1919 / 1701.5 | 0 / 0 |

### sio2

Fittable pairs: 0; aggregate ΔG MAE: — eV/atom; aggregate sign accuracy: —.

| Pair | status | train mean ΔG (eV/atom) | ΔG MAE | ΔG RMSE | sign accuracy | reference / predicted Tc (K) | false / missed |
|---|---|---:|---:|---:|---:|---|---:|
| sio2:cristobalite_beta_minus_quartz_beta | not_fittable_without_pair_training_labels | — | 0.00320909 | 0.0037702 | 0 | 1542.71 / 1675 | 0 / 0 |
| sio2:cristobalite_beta_minus_tridymite_p63mmc | not_fittable_without_pair_training_labels | — | 0.00025389 | 0.000260579 | 0 | — / 1675 | 1 / 0 |
| sio2:quartz_beta_minus_tridymite_p63mmc | not_fittable_without_pair_training_labels | — | 0.00313832 | 0.00367457 | 0 | 1582.32 / 1675 | 0 / 0 |

## Deviations from design

No new model was trained. The design evaluator's phase-level G MAE/coverage cannot be applied to a pair-only ΔG predictor without inventing a phase-energy gauge, so they are explicitly N/A. In phase-LOPO and system-LOSO, a held phase/system has no pairwise training ΔG labels; `constant_delta_g` is consequently reported as unavailable rather than fitted with held-out reference values.

