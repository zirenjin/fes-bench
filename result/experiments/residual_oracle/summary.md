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

For every pair, ΔG and ΔE use the physical high-minus-low direction. For the metal rows, hcp is low and bcc is high; ΔE_DPA is therefore E_DPA(bcc) − E_DPA(hcp).

## Delta-E budget

DFT relaxed-static ΔE is the primary baseline below. ΔE_eff is only a reference extrapolation from the canonical reference curves.

| Pair (low/high) | Delta E DFT baseline | Delta E DPA | Delta E_eff reference [range] | DPA minus DFT | DPA/DFT sign | Reference Delta G SD |
| --- | ---: | ---: | ---: | ---: | --- | ---: |
| hf hcp/bcc | 182.793 | 152.150 | 98.839 [64.196, 132.365] | -30.643 | yes | 11.987 |
| ti hcp/bcc | 88.853 | 93.749 | 64.237 [30.909, 96.317] | 4.896 | yes | 6.111 |
| zr hcp/bcc | 70.067 | 65.487 | 59.543 [32.904, 82.505] | -4.580 | yes | 8.581 |
| sio2 quartz_beta/cristobalite_beta | 11.716 | -9.996 | 9.955 [8.003, 11.590] | -21.712 | no | 3.607 |
| sio2 quartz_beta/tridymite_p63mmc | 17.415 | -7.361 | 9.596 [7.905, 11.731] | -24.775 | no | 3.576 |

The effective Delta E range is the six positive extrapolations from F1-F3 and the full/lowest-300-K windows. Values are an effective classical reference extrapolation, not a zero-point-inclusive ground-state energy.

The DFT baseline shows that the two SiO2 DPA static-energy differences have the wrong sign: quartz_beta → cristobalite_beta is DPA -9.996 meV/atom versus DFT +11.716 meV/atom, and quartz_beta → tridymite_p63mmc is DPA -7.361 meV/atom versus DFT +17.415 meV/atom.

## R3 crossing errors in the test region

| Pair | Reference Tc (K) | B0 error (K) | B1 error (K) | B2 error (K) | B0/B1/B2 test MAE (meV/atom) |
| --- | ---: | ---: | ---: | ---: | --- |
| hf hcp/bcc | 1919.000 | -6.108 | -6.085 | -6.048 | 0.526 / 0.524 / 0.520 |
| ti hcp/bcc | 1012.000 | -3.195 | -3.173 | -3.154 | 1.033 / 1.104 / 0.999 |
| zr hcp/bcc | 975.500 | -1.905 | -1.848 | -1.747 | 0.692 / 0.674 / 0.641 |
| sio2 quartz_beta/cristobalite_beta | 1542.709 | -1.800 | -1.759 | 14.632 | 0.149 / 0.147 / 1.049 |
| sio2 quartz_beta/tridymite_p63mmc | 1582.322 | -0.972 | -0.915 | 23.341 | 0.047 / 0.045 / 1.151 |

The reference-curve full-range R3 control RMSE is:
- hf hcp/bcc: 0.089 meV/atom
- ti hcp/bcc: 0.057 meV/atom
- zr hcp/bcc: 0.059 meV/atom
- sio2 quartz_beta/cristobalite_beta: 0.011 meV/atom
- sio2 quartz_beta/tridymite_p63mmc: 0.002 meV/atom

## Conclusions

- **hf hcp/bcc**: DFT baseline ΔE = 182.793 meV/atom, DPA ΔE = 152.150 meV/atom, and |DPA − DFT| = 30.643 meV/atom; the signs agree. The positive ΔE_eff value is reference-only.
- **ti hcp/bcc**: DFT baseline ΔE = 88.853 meV/atom, DPA ΔE = 93.749 meV/atom, and |DPA − DFT| = 4.896 meV/atom; the signs agree. The positive ΔE_eff value is reference-only.
- **zr hcp/bcc**: DFT baseline ΔE = 70.067 meV/atom, DPA ΔE = 65.487 meV/atom, and |DPA − DFT| = 4.580 meV/atom; the signs agree. The positive ΔE_eff value is reference-only.
- **sio2 quartz_beta/cristobalite_beta**: DFT baseline ΔE = 11.716 meV/atom, DPA ΔE = -9.996 meV/atom, and |DPA − DFT| = 21.712 meV/atom; the signs are opposite. The positive ΔE_eff value is reference-only.
- **sio2 quartz_beta/tridymite_p63mmc**: DFT baseline ΔE = 17.415 meV/atom, DPA ΔE = -7.361 meV/atom, and |DPA − DFT| = 24.775 meV/atom; the signs are opposite. The positive ΔE_eff value is reference-only.

- **hf hcp/bcc**, equal three-parameter residual form: B0 |Tc error| = 6.108 K; B1 = 6.085 K; B2 = 6.048 K. Relative to B0: B1: no meaningful difference; B2: no meaningful difference.
- **ti hcp/bcc**, equal three-parameter residual form: B0 |Tc error| = 3.195 K; B1 = 3.173 K; B2 = 3.154 K. Relative to B0: B1: no meaningful difference; B2: no meaningful difference.
- **zr hcp/bcc**, equal three-parameter residual form: B0 |Tc error| = 1.905 K; B1 = 1.848 K; B2 = 1.747 K. Relative to B0: B1: no meaningful difference; B2: no meaningful difference.
- **sio2 quartz_beta/cristobalite_beta**, equal three-parameter residual form: B0 |Tc error| = 1.800 K; B1 = 1.759 K; B2 = 14.632 K. Relative to B0: B1: no meaningful difference; B2: worse.
- **sio2 quartz_beta/tridymite_p63mmc**, equal three-parameter residual form: B0 |Tc error| = 0.972 K; B1 = 0.915 K; B2 = 23.341 K. Relative to B0: B1: no meaningful difference; B2: worse.

- **hf hcp/bcc**: compared with B0, B1 RMS ratio 23.277, B2 RMS ratio 21.725; hcp qh_reliable=True, bcc qh_reliable=False.
- **ti hcp/bcc**: compared with B0, B1 RMS ratio 22.143, B2 RMS ratio 51.244; hcp qh_reliable=True, bcc qh_reliable=False.
- **zr hcp/bcc**: compared with B0, B1 RMS ratio 15.127, B2 RMS ratio 32.609; hcp qh_reliable=True, bcc qh_reliable=False.
- **sio2 quartz_beta/cristobalite_beta**: compared with B0, B1 RMS ratio 5.566, B2 RMS ratio 301.882; quartz_beta qh_reliable=True, cristobalite_beta qh_reliable=False.
- **sio2 quartz_beta/tridymite_p63mmc**: compared with B0, B1 RMS ratio 5.781, B2 RMS ratio 309.808; quartz_beta qh_reliable=True, tridymite_p63mmc qh_reliable=False.

## Scope and deviations

- The cristobalite-beta/tridymite-beta pair has no reference crossing and was skipped.
- The optional classical-QH variants B1c/B2c were not run because the canonical QH artifacts preserve scalar F_QH outputs and minimum-frequency diagnostics, not the full positive phonon-frequency spectrum needed for a fresh classical vibrational sum.
- Canonical QH raw-run paths under `result/experiments/` were used because the current processed tree does not materialize every phase QH CSV; no legacy, archive, or invalid source was read.
