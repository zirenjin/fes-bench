# DFT static-energy error budget

The sign convention is delta E = E_high - E_low. Energies are in meV/atom.

## Numerical conclusions

- hf:hcp:bcc: |DPA - DFT| = 30.642889; reference Delta-G std = 11.991932. DPA static energy difference is the main error source for this pair.
  DPA = 152.149826; DFT@DPA = 181.811795; relaxed DFT = 182.792715; model error = -30.642889; structural error = -0.980920. Delta-bias is same-sign and has a magnitude ratio of 0.0685 relative to the DPA static error.
- ti:hcp:bcc: |DPA - DFT| = 4.895724; reference Delta-G std = 6.115905. The DPA static energy difference does not exceed the reference Delta-G standard deviation for this pair.
  DPA = 93.748719; DFT@DPA = 110.311845; relaxed DFT = 88.852995; model error = 4.895724; structural error = 21.458850. No usable snapshot bias was available for this pair.
- zr:hcp:bcc: |DPA - DFT| = 4.579802; reference Delta-G std = 8.585435. The DPA static energy difference does not exceed the reference Delta-G standard deviation for this pair.
  DPA = 65.486848; DFT@DPA = 84.987880; relaxed DFT = 70.066650; model error = -4.579802; structural error = 14.921230. No usable snapshot bias was available for this pair.
- sio2:quartz_beta:cristobalite_beta: |DPA - DFT| = 21.711989; reference Delta-G std = 3.607960. DPA static energy difference is the main error source for this pair.
  DPA = -9.995852; DFT@DPA = 28.648550; relaxed DFT = 11.716136; model error = -21.711989; structural error = 16.932413. Delta-bias is same-sign and has a magnitude ratio of 1.18 relative to the DPA static error.
- sio2:quartz_beta:tridymite_p63mmc: |DPA - DFT| = 24.775472; reference Delta-G std = 3.577461. DPA static energy difference is the main error source for this pair.
  DPA = -7.360924; DFT@DPA = 27.667406; relaxed DFT = 17.414548; model error = -24.775472; structural error = 10.252858. Delta-bias is same-sign and has a magnitude ratio of 1.04 relative to the DPA static error.

Across the five pairs, the mean absolute model error is 17.321175 meV/atom and the mean absolute structural error is 12.909254 meV/atom.

## Convergence and deviations

The convergence target meshes are 24x24x24 for Hf, Ti, and Zr, and 5x5x5 for SiO2. The available production final energies use the existing k666 metal and k444 SiO2 outputs because final-adopted-converged outputs were not present; this fallback is recorded in findings.meta.json. The adopted metal ENCUT target is 2.0 times the reference value because the required 1.3x check remained above 0.2 meV/atom; higher-cutoff checks are reported in convergence.csv.
Missing convergence levels are retained as n/a:missing_output in convergence.csv rather than being inferred from another k-point or ENCUT.
All nine relaxed structures retain their target space group according to spglib.
The actual VASP executable reported revision 5.4.4.18Apr17-6-g9f103f2a35, build Feb 06 2024. Reference archives reported the same executable revision family but different build metadata for Hf and SiO2; this is recorded in reference_settings.csv and findings.meta.json.
