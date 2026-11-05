# FES-Bench

```bash
pip install -r requirements.txt
```

Run every command from the repository root.

## Reproduce from stored curves (seconds, no GPU, no DFT)

| Test | Paper | Command | Output |
|---|---|---|---|
| Database polymorph pairs: metrics, tables, figures | main table, budget / profile / workflow / intervention figures, SI tables and route figure | `python database_pairs/make_main.py` | `database_pairs/tables/`, `database_pairs/figures/` |
| Crossing-budget numerical checks | SI, numerical robustness | `python database_pairs/audit_numerics.py` | `database_pairs/audit/NUMERICS.json` |
| Descriptor baseline (Bartel 2018) | SI descriptor table | `python analyses/descriptor_baseline/bartel_baseline.py` | `analyses/descriptor_baseline/results/BARTEL.json` |
| Descriptor vs DPA-3.1 on common pairs (run after the descriptor baseline) | SI descriptor paragraph | `python analyses/descriptor_baseline/common_subset.py` | `analyses/descriptor_baseline/results/BARTEL_COMMON_SUBSET.json` |
| Error scales and example curves (run after the descriptor baseline) | error-scale and example-curve figures | `python database_pairs/make_context_figs.py` | `database_pairs/figures/` |
| Finite-temperature reference systems | finite-temperature table, SI per-phase table | `python finite_temperature/scripts/a_layer_score.py finite_temperature/data data/processed finite_temperature/results/A_LAYER.json` | `finite_temperature/results/A_LAYER.json` |
| SiO2 reference-level sensitivity | introduction silica figure | `python finite_temperature/scripts/sio2_functionals.py` | `finite_temperature/results/`, `finite_temperature/figures/` |
| Curvature correction summary | curvature figure | `python analyses/curvature_correction/scripts/make_fig5.py analyses/curvature_correction/figures` | `analyses/curvature_correction/figures/` |

## Recompute from scratch (GPU and source data needed)

| Test | Scripts | Needs |
|---|---|---|
| Selection of the database pairs and model phonons | `analyses/selection_and_covariance/scripts/` (`select_entries.py`, `mlip_vs_dft_phonons.py`, `analyze_covariance.py`) | Alexandria phonon files, DPA-3.1 checkpoint |
| Second-model phonons and paired contrast | `analyses/second_model_replication/scripts/` | Alexandria phonon files, MatterSim checkpoint |
| Stored curves and manifest | `database_pairs/build.py` | outputs of the two rows above |
| Relaxed-geometry pipeline | `analyses/relaxed_geometry/scripts/` | DPA-3.1 checkpoint |
| Crossing sensitivity and tail audit | `analyses/crossing_sensitivity/scripts/` | stored curves, Alexandria phonon files |
| Learned correction routes A–C | `analyses/correction_routes/scripts/` | stored curves, spectra |
| Anchored temperature response | `analyses/anchored_response/scripts/` | stored curves |
| Curvature correction (gates 1, 1b, 1c), reference efficiency (RE1), identifiability (RE1-MECH) | `analyses/curvature_correction/scripts/` (`gate1*.py`, `re1.py`, `re1_mech.py`, then `*_score.py`) | Alexandria phonon files, both checkpoints, GPU |
| RE1 subset audit (paired displacements drawn together) | `analyses/curvature_correction/scripts/pm_pair_audit.py <Alexandria raw dir> database_pairs/manifest/ENTRIES.csv analyses/curvature_correction/results/PM_PAIR_AUDIT.json` | Alexandria phonon files (about 10 min) |
| Finite-temperature model curves | `finite_temperature/scripts/a_layer_phonons.py data/processed finite_temperature/data` | both checkpoints, GPU |

Each analysis folder contains its frozen `PROTOCOL.md` (with `.sha256`), `DEVIATIONS.md` where applicable, and `RESULT.md`.
