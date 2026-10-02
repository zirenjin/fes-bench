# T3 SiO2 polynomial seed 11 diagnosis

Status: the pre-fix run is invalid and was not admitted to a predictor table.
The corrected isolated rerun passes the training-region ΔG-MAE gate; no other
T3 seeds were started.

## Checks

| item | finding |
|---|---|
| structures | The generated training view uses the current canonical representatives: cristobalite 24 atoms, quartz 9 atoms, tridymite 12 atoms. The builder records their structure SHA-256 values in its provenance; no 192/243/216-atom MD snapshot is present. |
| source labels | `reference_G.csv:G_eV_per_atom` is read in eV/atom. The DeepMD FES loss divides labels by atom count internally, so the training/evaluation view now converts `N*G` to eV/cell and records `source_label_unit=eV/atom`. This conversion was missing in the failed run. |
| checkpoint/head | DPA-3.1-3M, `Domains_SSE_PBE`; checkpoint SHA-256 `86dd3a804d78ca5d203ebf98747e8f16dff9713ba8950097ceb760b161e19907`. |
| type map | The selected checkpoint branch uses the full periodic-table map (O index 7, Si index 13). The failed view wrote compact 0/1 IDs. A direct quartz single-point changed from `+1.5163857` eV/atom with compact IDs to `−4.0995074` eV/atom with checkpoint IDs, matching the canonical `E_DPA` record. |
| calibration | One constant was fit on training rows only. The corrected diagnostic gives `c_s = −9.3696092` eV/atom; this is a gauge offset, not a test-derived quantity. |
| optimizer | E1 settings were retained: polynomial basis, 3000 steps, exponential LR 1e−3 → 1e−5 at step 3000, seed 11. |

## Before/after gate

The invalid pre-fix result reported training ΔG MAE **79.35 meV/atom** and
test absolute-G MAE **11.638 eV/atom**. It is preserved under
`result/experiments/legacy_support/t3_invalid_type_map_20261001/`.

After fixing the checkpoint type IDs and the FES cell-label conversion, the
isolated polynomial/SiO2/seed-11 diagnostic reports training ΔG MAE
**1.850 meV/atom** (G MAE 98.41 meV/atom), passing the requested ≤20 meV/atom
training gate.  The canonical test view gives ΔG MAE **25.737 meV/atom**,
skill **−10.263**, sign accuracy **0.07569**, and G MAE **10.532 eV/atom**;
these are diagnostic only and are not promoted to the main table.  The
lightweight handoff is in `diagnostic_seed11/` (the checkpoint remains
external).  No batch T3/T2 run starts in this step.

The source split is unchanged and the train/test key intersection is empty.
