# Residual Oracle Diagnostics

This report uses only canonical processed reference data, the frozen v2 temperature split, provenance tables, and canonical QH outputs. No model was trained or evaluated by this experiment.

## Input consistency

All 9 in-scope phases have matching checkpoint SHA-256 and head provenance between E_DPA, F_QH, and the head policy: **True**.

The QH CSV schema and `src/lib/fes_bench/qh/run.py` confirm that `F_QH` already includes `E_static + F_vib`; the baselines therefore use F_QH directly and do not add E_DPA again. The QH thermal properties are the quantum phonopy free energy and include zero-point energy. The QH runner records imaginary-mode diagnostics with a -0.05 THz cutoff; this script does not edit or recalculate those modes.

## Metal static-energy audit

The reference crossing slope identifies hcp as the low-temperature phase for all three metals. The table prints the canonical E_DPA provenance used below.

| System | Phase | E_DPA (eV/atom) | Source file | Checkpoint SHA-256 | Head |
| --- | --- | ---: | --- | --- | --- |
| hf | hcp | -7.194502562284 | `external/data/raw/hf/unpacked/Hf_hcp_PBE/e.effective_QH_potentials/a.20K_4Shells_mediumCell_writeForceConstPot/POSCAR_supercell_1st_alat.vasp` | `86dd3a804d78ca5d203ebf98747e8f16dff9713ba8950097ceb760b161e19907` | Domains_Alloy |
| hf | bcc | -7.042352735996 | `external/data/raw/hf/unpacked/Hf_bcc_PBE/e.effective_QH_potentials/b.melting_T_5Shells_mediumCell_writeForceConstPot/POSCAR_supercell_1st_alat.vasp` | `86dd3a804d78ca5d203ebf98747e8f16dff9713ba8950097ceb760b161e19907` | Domains_Alloy |
| ti | hcp | -6.432060956955 | `structure.source.extxyz` | `86dd3a804d78ca5d203ebf98747e8f16dff9713ba8950097ceb760b161e19907` | Domains_Alloy |
| ti | bcc | -6.338312238455 | `structure.source.extxyz` | `86dd3a804d78ca5d203ebf98747e8f16dff9713ba8950097ceb760b161e19907` | Domains_Alloy |
| zr | hcp | -7.152492448688 | `structure.source.extxyz` | `86dd3a804d78ca5d203ebf98747e8f16dff9713ba8950097ceb760b161e19907` | Domains_Alloy |
| zr | bcc | -7.087005600333 | `structure.source.extxyz` | `86dd3a804d78ca5d203ebf98747e8f16dff9713ba8950097ceb760b161e19907` | Domains_Alloy |

For the metal rows, hcp is low and bcc is high; ΔE_DPA is therefore E_DPA(bcc) − E_DPA(hcp), with the same high-minus-low direction used for ΔE_eff.

## Delta-E budget

| Pair (low/high) | Delta E DPA (meV/atom) | Delta E effective median [range] (meV/atom) | DPA minus effective | Reference Delta G SD | |error| / SD |
| --- | ---: | ---: | ---: | ---: | ---: |
| hf hcp/bcc | 152.150 | -98.839 [-132.365, -64.196] | 250.988 | 11.987 | 20.938 |
| ti hcp/bcc | 93.749 | -64.237 [-96.317, -30.909] | 157.986 | 6.111 | 25.851 |
| zr hcp/bcc | 65.487 | -59.543 [-82.505, -32.904] | 125.030 | 8.581 | 14.571 |
| sio2 quartz_beta/cristobalite_beta | -9.996 | -9.955 [-11.590, -8.003] | -0.041 | 3.607 | 0.011 |
| sio2 quartz_beta/tridymite_p63mmc | -7.361 | -9.596 [-11.731, -7.905] | 2.235 | 3.576 | 0.625 |

The effective Delta E range is the six extrapolations from F1-F3 and the full/lowest-300-K windows. Values are an effective classical reference extrapolation, not a zero-point-inclusive ground-state energy.

## R3 crossing errors in the test region

| Pair | Reference Tc (K) | B0 error (K) | B1 error (K) | B2 error (K) | B0/B1/B2 test MAE (meV/atom) |
| --- | ---: | ---: | ---: | ---: | --- |
| hf hcp/bcc | 1919.000 | -6.108 | -6.132 | -6.169 | 0.526 / 0.528 / 0.531 |
| ti hcp/bcc | 1012.000 | -3.195 | -3.217 | -3.236 | 1.033 / 0.961 / 1.067 |
| zr hcp/bcc | 975.500 | -1.905 | -1.963 | -2.063 | 0.692 / 0.711 / 0.744 |
| sio2 quartz_beta/cristobalite_beta | 1542.709 | -1.800 | -1.842 | -14.487 | 0.149 / 0.151 / 1.347 |
| sio2 quartz_beta/tridymite_p63mmc | 1582.322 | -0.972 | -1.028 | -18.735 | 0.047 / 0.049 / 1.245 |

The reference-curve full-range R3 control RMSE is:
- hf hcp/bcc: 0.089 meV/atom
- ti hcp/bcc: 0.057 meV/atom
- zr hcp/bcc: 0.059 meV/atom
- sio2 quartz_beta/cristobalite_beta: 0.011 meV/atom
- sio2 quartz_beta/tridymite_p63mmc: 0.002 meV/atom

## Conclusions

- **hf hcp/bcc**: the absolute Delta-E discrepancy is larger than the reference Delta G standard deviation (20.938 times the standard deviation). Static energy difference is the main error source.
- **ti hcp/bcc**: the absolute Delta-E discrepancy is larger than the reference Delta G standard deviation (25.851 times the standard deviation). Static energy difference is the main error source.
- **zr hcp/bcc**: the absolute Delta-E discrepancy is larger than the reference Delta G standard deviation (14.571 times the standard deviation). Static energy difference is the main error source.
- **sio2 quartz_beta/cristobalite_beta**: the absolute Delta-E discrepancy is smaller or comparable than the reference Delta G standard deviation (0.011 times the standard deviation).
- **sio2 quartz_beta/tridymite_p63mmc**: the absolute Delta-E discrepancy is smaller or comparable than the reference Delta G standard deviation (0.625 times the standard deviation).

- **hf hcp/bcc**, equal three-parameter residual form: B0 |Tc error| = 6.108 K; B1 = 6.132 K; B2 = 6.169 K. Relative to B0: B1: no meaningful difference; B2: no meaningful difference.
- **ti hcp/bcc**, equal three-parameter residual form: B0 |Tc error| = 3.195 K; B1 = 3.217 K; B2 = 3.236 K. Relative to B0: B1: no meaningful difference; B2: no meaningful difference.
- **zr hcp/bcc**, equal three-parameter residual form: B0 |Tc error| = 1.905 K; B1 = 1.963 K; B2 = 2.063 K. Relative to B0: B1: no meaningful difference; B2: no meaningful difference.
- **sio2 quartz_beta/cristobalite_beta**, equal three-parameter residual form: B0 |Tc error| = 1.800 K; B1 = 1.842 K; B2 = 14.487 K. Relative to B0: B1: no meaningful difference; B2: worse.
- **sio2 quartz_beta/tridymite_p63mmc**, equal three-parameter residual form: B0 |Tc error| = 0.972 K; B1 = 1.028 K; B2 = 18.735 K. Relative to B0: B1: no meaningful difference; B2: worse.

- **hf hcp/bcc**: compared with B0, B1 RMS ratio 24.104, B2 RMS ratio 21.149; hcp qh_reliable=True, bcc qh_reliable=False.
- **ti hcp/bcc**: compared with B0, B1 RMS ratio 20.684, B2 RMS ratio 52.771; hcp qh_reliable=True, bcc qh_reliable=False.
- **zr hcp/bcc**: compared with B0, B1 RMS ratio 14.406, B2 RMS ratio 33.468; hcp qh_reliable=True, bcc qh_reliable=False.
- **sio2 quartz_beta/cristobalite_beta**: compared with B0, B1 RMS ratio 6.558, B2 RMS ratio 302.957; quartz_beta qh_reliable=True, cristobalite_beta qh_reliable=False.
- **sio2 quartz_beta/tridymite_p63mmc**: compared with B0, B1 RMS ratio 6.661, B2 RMS ratio 310.773; quartz_beta qh_reliable=True, tridymite_p63mmc qh_reliable=False.

## Scope and deviations

- The cristobalite-beta/tridymite-beta pair has no reference crossing and was skipped.
- The optional classical-QH variants B1c/B2c were not run because the canonical QH artifacts preserve scalar F_QH outputs and minimum-frequency diagnostics, not the full positive phonon-frequency spectrum needed for a fresh classical vibrational sum.
- Canonical QH raw-run paths under `result/experiments/` were used because the current processed tree does not materialize every phase QH CSV; no legacy, archive, or invalid source was read.
