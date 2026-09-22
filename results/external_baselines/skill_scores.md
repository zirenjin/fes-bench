# External baseline report

Bartel coefficient audit: the original SI/Eq. 4 values and implementation values are identical: −2.48×10⁻⁴ for ln(V), −8.94×10⁻⁵ for m/V, +0.181 ln(T), and −0.882. For fixed-composition polymorphs the reduced mass is identical, so only relaxed per-atom volume distinguishes phases; this is the baseline's intrinsic limitation.

`global_mean_delta_g` is the mean of every available training pair label in a fold. It is pair-only and differs from `constant_delta_g`, which fits one constant separately for each pair.

## Skill relative to global_mean_delta_g

Skill is `1 − MAE_baseline / MAE_global_mean_delta_g`, using only pairs for which both values are available. Negative values are explicitly marked.

| split/fold | baseline | matched pair MAE (eV/atom) | global floor MAE (eV/atom) | skill |
|---|---|---:|---:|---:|
| temp_extrap/all | bartel2018 | 0.0194312 | 0.00505025 | **-2.848 (negative)** |
| temp_extrap/all | interp_const | 0.001969 | 0.00638272 | 0.692 |
| temp_extrap/all | phase_id_mlp | — | — | unavailable |
| phase_lopo/hf:bcc | bartel2018 | — | — | unavailable |
| phase_lopo/hf:hcp | bartel2018 | — | — | unavailable |
| phase_lopo/sio2:cristobalite_beta | bartel2018 | 0.01866 | 0.00472384 | **-2.950 (negative)** |
| phase_lopo/sio2:quartz_beta | bartel2018 | 0.0250602 | 0.00418241 | **-4.992 (negative)** |
| phase_lopo/sio2:tridymite_p63mmc | bartel2018 | 0.0145734 | 0.0033239 | **-3.384 (negative)** |
| phase_lopo/hf:bcc | interp_const | — | — | unavailable |
| phase_lopo/hf:hcp | interp_const | — | — | unavailable |
| phase_lopo/sio2:cristobalite_beta | interp_const | — | — | unavailable |
| phase_lopo/sio2:quartz_beta | interp_const | — | — | unavailable |
| phase_lopo/sio2:tridymite_p63mmc | interp_const | — | — | unavailable |
| phase_lopo/hf:bcc | phase_id_mlp | — | — | unavailable |
| phase_lopo/hf:hcp | phase_id_mlp | — | — | unavailable |
| phase_lopo/sio2:cristobalite_beta | phase_id_mlp | — | — | unavailable |
| phase_lopo/sio2:quartz_beta | phase_id_mlp | — | — | unavailable |
| phase_lopo/sio2:tridymite_p63mmc | phase_id_mlp | — | — | unavailable |
| system_loso/hf | bartel2018 | — | — | unavailable |
| system_loso/sio2 | bartel2018 | 0.0194312 | 0.0083536 | **-1.326 (negative)** |
| system_loso/hf | interp_const | — | — | unavailable |
| system_loso/sio2 | interp_const | — | — | unavailable |
| system_loso/hf | phase_id_mlp | — | — | unavailable |
| system_loso/sio2 | phase_id_mlp | — | — | unavailable |

## Environment and deviations

No new benchmark model was trained. `phase_id_mlp` is marked unavailable in this run because the remote torch installation fails while loading `libtorch_global_deps.so`; the runner continues and records the other baselines. Bartel Hf rows are unavailable because the canonical Hf meta files do not contain a DPA 0-K representative energy; no reference G value was substituted.
