# External baseline report

Bartel coefficient audit: the original SI/Eq. 4 values and implementation values are identical: −2.48×10⁻⁴ for ln(V), −8.94×10⁻⁵ for m/V, +0.181 ln(T), and −0.882. For fixed-composition polymorphs the reduced mass is identical, so only relaxed per-atom volume distinguishes phases; this is the baseline's intrinsic limitation.

`global_mean_delta_g` is the mean of every available training pair label in a fold. It is pair-only and differs from `constant_delta_g`, which fits one constant separately for each pair.

## Skill relative to the zero floor

Skill is `1 − MAE_baseline / MAE_zero`, where `MAE_zero` is the mean reference |ΔG| on the evaluated points. Negative values are explicitly marked.

| split/fold | baseline | matched pair MAE (eV/atom) | zero floor MAE (eV/atom) | skill |
|---|---|---:|---:|---:|
| temp_extrap/all | bartel2018 | 0.0719694 | 0.00408617 | **-16.613 (negative)** |
| temp_extrap/all | interp_const | 0.00580824 | 0.00408617 | **-0.421 (negative)** |
| temp_extrap/all | phase_id_mlp | — | — | unavailable |
| phase_lopo/hf:bcc | bartel2018 | 0.128374 | 0.0113998 | **-10.261 (negative)** |
| phase_lopo/hf:hcp | bartel2018 | 0.128374 | 0.0113998 | **-10.261 (negative)** |
| phase_lopo/sio2:cristobalite_beta | bartel2018 | 0.0300918 | 0.00173149 | **-16.379 (negative)** |
| phase_lopo/sio2:quartz_beta | bartel2018 | 0.0564449 | 0.00317371 | **-16.785 (negative)** |
| phase_lopo/sio2:tridymite_p63mmc | bartel2018 | 0.0288456 | 0.00169611 | **-16.007 (negative)** |
| phase_lopo/ti:bcc | bartel2018 | 0.0985418 | 0.00669683 | **-13.715 (negative)** |
| phase_lopo/ti:hcp | bartel2018 | 0.0985418 | 0.00669683 | **-13.715 (negative)** |
| phase_lopo/zr:bcc | bartel2018 | 0.0627469 | 0.00759824 | **-7.258 (negative)** |
| phase_lopo/zr:hcp | bartel2018 | 0.0627469 | 0.00759824 | **-7.258 (negative)** |
| phase_lopo/hf:bcc | interp_const | — | — | unavailable |
| phase_lopo/hf:hcp | interp_const | — | — | unavailable |
| phase_lopo/sio2:cristobalite_beta | interp_const | — | — | unavailable |
| phase_lopo/sio2:quartz_beta | interp_const | — | — | unavailable |
| phase_lopo/sio2:tridymite_p63mmc | interp_const | — | — | unavailable |
| phase_lopo/ti:bcc | interp_const | — | — | unavailable |
| phase_lopo/ti:hcp | interp_const | — | — | unavailable |
| phase_lopo/zr:bcc | interp_const | — | — | unavailable |
| phase_lopo/zr:hcp | interp_const | — | — | unavailable |
| phase_lopo/hf:bcc | phase_id_mlp | — | — | unavailable |
| phase_lopo/hf:hcp | phase_id_mlp | — | — | unavailable |
| phase_lopo/sio2:cristobalite_beta | phase_id_mlp | — | — | unavailable |
| phase_lopo/sio2:quartz_beta | phase_id_mlp | — | — | unavailable |
| phase_lopo/sio2:tridymite_p63mmc | phase_id_mlp | — | — | unavailable |
| phase_lopo/ti:bcc | phase_id_mlp | — | — | unavailable |
| phase_lopo/ti:hcp | phase_id_mlp | — | — | unavailable |
| phase_lopo/zr:bcc | phase_id_mlp | — | — | unavailable |
| phase_lopo/zr:hcp | phase_id_mlp | — | — | unavailable |
| system_loso/hf | bartel2018 | 0.128374 | 0.0113998 | **-10.261 (negative)** |
| system_loso/ti | bartel2018 | 0.0985418 | 0.00669683 | **-13.715 (negative)** |
| system_loso/zr | bartel2018 | 0.0627469 | 0.00759824 | **-7.258 (negative)** |
| system_loso/hf | interp_const | — | — | unavailable |
| system_loso/ti | interp_const | — | — | unavailable |
| system_loso/zr | interp_const | — | — | unavailable |
| system_loso/hf | phase_id_mlp | — | — | unavailable |
| system_loso/ti | phase_id_mlp | — | — | unavailable |
| system_loso/zr | phase_id_mlp | — | — | unavailable |

## Environment and deviations

No new benchmark model was trained. Bartel now covers Hf using the Domains_Alloy-computed Hf E0; its Hf error is substantially larger than its SiO₂ error (the SiO₂-only value is 28.2 meV/atom). The independent thu-GenSi torch environment runs the two torch-dependent tests successfully.
