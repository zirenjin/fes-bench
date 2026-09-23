# Phase-ID MLP mat-agent rerun

The phase-ID MLP was rerun in the working `/root/miniconda3/envs/mat-agent` environment on thu-GenSi, without training a benchmark model. The command was:

```text
PYTHONPATH=. /root/miniconda3/envs/mat-agent/bin/python scripts/run_external_baselines.py --data-root data --output-root results/external_baselines_matagent --splits splits/temp_extrap.json splits/phase_lopo.json splits/system_loso.json
```

Skill is `1 - MAE_phase_id_mlp / MAE_global_mean_delta_g`, using the matched pair set.

| split/fold | matched pairs | phase-ID pair MAE (eV/atom) | global floor MAE (eV/atom) | skill |
|---|---:|---:|---:|---:|
| temp_extrap/all | 4 | 0.0118580893 | 0.0063827153 | **−0.858 (negative)** |
| phase_lopo/hf:bcc | 0 | — | — | unavailable: held-out phase has no training rows |
| phase_lopo/hf:hcp | 0 | — | — | unavailable: held-out phase has no training rows |
| phase_lopo/sio2:cristobalite_beta | 0 | — | — | unavailable: held-out phase has no training rows |
| phase_lopo/sio2:quartz_beta | 0 | — | — | unavailable: held-out phase has no training rows |
| phase_lopo/sio2:tridymite_p63mmc | 0 | — | — | unavailable: held-out phase has no training rows |
| system_loso/hf | 0 | — | — | unavailable: held-out system has no training rows |
| system_loso/sio2 | 0 | — | — | unavailable: held-out system has no training rows |

The four temp-extrap pair MAEs were 0.0068977278 (Hf bcc−hcp), 0.0151739052 (cristobalite−quartz), 0.0050935862 (cristobalite−tridymite), and 0.0202671381 (quartz−tridymite) eV/atom. The remote raw files remain in the isolated workspace under `results/external_baselines_matagent/metrics.json`; the local baseline report records the portable summary and sign.

## Deviations from design

The local base environment could not load its broken Torch library, so this one rerun used the already-configured mat-agent environment. No data, split, head policy, or benchmark model was changed.
