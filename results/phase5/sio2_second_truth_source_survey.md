# SiO2 second-tier truth-source survey (E9 preflight)

## Finding

The Forslund et al. 2026 cross-functional FEP work is publicly available through the DaRUS record [10.18419/DARUS-4999](https://darus.uni-stuttgart.de/dataset.xhtml?persistentId=doi%3A10.18419%2FDARUS-4999), associated with the 2026 *npj Computational Materials* article. The record is not a separate small CSV truth table: it is a 184.2 MB dataset containing two ZIP files.

The relevant formats are:

- `direct_upsampling_data_rungs_1-3.zip` (134.1 MB): effective harmonic potentials, MTPs, training OUTCARs, and thermodynamic-property tables for β-quartz, β-cristobalite, and the two investigated tridymites, across LDA, PBE, r2SCAN, PBE-D3(BJ), vdW-DF-cx, r2SCAN-D3(BJ), r2SCAN-D4, and r2SCAN-rVV10.
- `input_for_functional_evaluation.zip` (50.1 MB): inputs for the paper's cross-functional FEP evaluation procedure, rather than a ready-made canonical `reference_G.csv` table.

The dataset metadata identifies measured variables including Gibbs energy in meV/atom, entropy, equilibrium volume, and bulk modulus at 0.1 MPa and temperature in K. It is therefore a viable second truth/reference tier, but its tables must still be parsed and harmonized to the benchmark schema; no file was downloaded or processed in this survey.

## Scope decision

This is an availability/format survey only. It does not claim that the second tier is already aligned with the current PBE-D3(BJ) 851–2499 K grid, and it does not alter the frozen benchmark inputs.

## Deviations from design

No deviation: the requested preflight explicitly forbids downloading and processing this source.
