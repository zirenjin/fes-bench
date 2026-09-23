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

The raw files and `qh_summary.json` are in the remote path above. Phonopy's
`run_thermal_properties` default (`cutoff_frequency=None`, interpreted as
zero) excludes imaginary modes from the thermal free-energy sum, while the
runner separately records the raw minimum mesh frequency. The explicit
mode-fraction/BZ diagnosis below therefore reports reliability of the force
constants rather than silently treating the curves as dynamically stable.

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
tridymite −3.555 to −0.388 THz. The same-head diagnostic at volume scale 1.0
gives the following exact mode counts (threshold −0.05 THz):

| phase | min freq (THz) | negative modes / all modes | negative-mode fraction | negative q-points / q-points | negative q-point fraction | Γ-neighborhood fraction of negative modes | `qh_reliable` |
|---|---:|---:|---:|---:|---:|---:|---|
| quartz_beta | −2.63652 | 8,178 / 1,495,908 | 0.005467 | 2,051 / 2,052 | 0.999513 | 0 | **false** |
| cristobalite_beta | −1.70034 | 16,187 / 1,181,952 | 0.013695 | 2,052 / 2,052 | 1.000000 | 0.000494 | **false** |
| tridymite_p63mmc | −1.65362 | 19,359 / 1,329,696 | 0.014559 | 2,052 / 2,052 | 1.000000 | 0.000413 | **false** |

The total-mode denominator is `n_qpoints × n_modes`. Negative modes are distributed throughout the
Brillouin-zone mesh rather than being Γ-local, so no Γ exclusion or ASR
post-filter was applied. JSON evidence is in the isolated workspace at
`quartz_beta_imaginary_diagnosis.json`, `cristobalite_beta_imaginary_diagnosis.json`,
and `tridymite_p63mmc_imaginary_diagnosis.json`.

## Deviations from design

The first comparison attempt exposed a plotting bug in the generic 2×2 layout
for three phases; `fes_bench/qh/compare.py` was fixed to use a two-row,
`n_phases`-column layout, without changing any metric or split. The config's
`imaginary_frequency_cutoff_THz: -0.05` was not forwarded explicitly; the
phonopy default zero cutoff was used, which removes all negative modes and is
stricter than the configured threshold. This is recorded as a provenance
deviation; it does not change the raw reliability diagnosis.
