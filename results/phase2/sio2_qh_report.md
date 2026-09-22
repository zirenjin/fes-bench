# SiO₂ QH baseline — `Domains_SSE_PBE`

The three-phase QH run completed in the isolated thu-GenSi workspace
`/share/jzr/codex_fes_overnight_20260922_1790052610710037915`. It used the frozen
DPA checkpoint `DPA-3.1-3M.pt`, head `Domains_SSE_PBE`, a `[1,1,1]` supercell,
`16×16×16` mesh, five volume scales `[0.96, 0.98, 1.00, 1.02, 1.04]`, and
`0.01 Å` finite displacements. Each phase has 1,649 temperatures (851–2499 K).

## QH outputs

| phase | output | selected volume scale over grid | `F_QH` at 851 K (eV/atom) | `F_QH` at 1442 K | `F_QH` at 2499 K |
|---|---|---:|---:|---:|---:|
| quartz_beta | `quartz_beta_fqh.csv` | 1.02 | −4.14302964 | −4.39080988 | −4.97597516 |
| cristobalite_beta | `cristobalite_beta_fqh.csv` | 1.02 | −4.17545559 | −4.42817293 | −5.02251496 |
| tridymite_p63mmc | `tridymite_p63mmc_fqh.csv` | 1.02 | −4.16820945 | −4.42062884 | −5.01436235 |

The raw files and `qh_summary.json` are in the remote path above. The existing
QH runner reports the minimum mesh frequency by volume but does not itself
remove imaginary modes; the explicit mode-fraction/BZ diagnosis is being run
serially after this production calculation.

## Calibrated comparison to reference `G`

`qh.compare` calibrated one common additive constant on the frozen training
range through 1442 K (`c_system = −3.95599693 eV/atom`) and evaluated the full
grid. The generated remote artifacts are
`results/phase2_sio2_sse_pbe/sio2/qh_comparison.json` and
`fqh_vs_reference.{png,pdf}`.

| phase | test MAE (meV/atom) | test max error (meV/atom) | positive-slope steps |
|---|---:|---:|---:|
| quartz_beta | 44.2479 | 59.1139 | 0 |
| cristobalite_beta | 8.12018 | 20.0685 | 0 |
| tridymite_p63mmc | 13.2290 | 28.0303 | 0 |
| all phases | **21.8657** | **59.1139** | — |

## Reliability / imaginary-mode status

The QH summary already shows large negative minimum frequencies at every
volume: quartz −3.274 to −2.174 THz, cristobalite −3.619 to −1.026 THz, and
tridymite −3.555 to −0.388 THz. Therefore a reliable-QH decision is not
expected for any phase. The exact negative-mode fraction and Brillouin-zone
distribution will be recorded from the serial `qh.diagnose` outputs before
this line is marked complete.

## Deviations from design

The first comparison attempt exposed a plotting bug in the generic 2×2 layout
for three phases; `fes_bench/qh/compare.py` was fixed to use a two-row,
`n_phases`-column layout, without changing any metric or split. The QH runner's
documented imaginary-mode cutoff is not applied inside `qh.run`; this is why an
explicit diagnostic pass is required rather than silently calling the raw
curves reliable.
