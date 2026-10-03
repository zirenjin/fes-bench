# DFT static-energy error budget

The sign convention is delta E = E_high - E_low. Energies are in meV/atom.

## Numerical conclusions

- hf:hcp:bcc: |DPA - DFT| = 30.160429; reference Delta-G std = 11.991932. DPA static energy difference is the main error source for this pair.
  DPA = 152.149826; DFT@DPA = 181.811795; relaxed DFT = 182.310255; model error = -30.160429; structural error = -0.498460. Delta-bias is same-sign and has a magnitude ratio of 0.0696 relative to the DPA static error.
- ti:hcp:bcc: |DPA - DFT| = 16.325771; reference Delta-G std = 6.115905. DPA static energy difference is the main error source for this pair.
  DPA = 93.748719; DFT@DPA = 110.311845; relaxed DFT = 110.074490; model error = -16.325771; structural error = 0.237355. No usable snapshot bias was available for this pair.
- zr:hcp:bcc: |DPA - DFT| = 19.166422; reference Delta-G std = 8.585435. DPA static energy difference is the main error source for this pair.
  DPA = 65.486848; DFT@DPA = 84.987880; relaxed DFT = 84.653270; model error = -19.166422; structural error = 0.334610. No usable snapshot bias was available for this pair.
- sio2:quartz_beta:cristobalite_beta: |DPA - DFT| = 38.282017; reference Delta-G std = 3.607960. DPA static energy difference is the main error source for this pair.
  DPA = -9.995852; DFT@DPA = 28.648550; relaxed DFT = 28.286164; model error = -38.282017; structural error = 0.362386. Delta-bias is same-sign and has a magnitude ratio of 0.67 relative to the DPA static error.
- sio2:quartz_beta:tridymite_p63mmc: |DPA - DFT| = 34.654722; reference Delta-G std = 3.577461. DPA static energy difference is the main error source for this pair.
  DPA = -7.360924; DFT@DPA = 27.667406; relaxed DFT = 27.293798; model error = -34.654722; structural error = 0.373608. Delta-bias is same-sign and has a magnitude ratio of 0.743 relative to the DPA static error.

Across the five pairs, the mean absolute model error is 27.717872 meV/atom and the mean absolute structural error is 0.361284 meV/atom.

## Convergence and deviations

The adopted production meshes are 24x24x24 for Hf, Ti, and Zr, and 5x5x5 for SiO2. The production final energies use these adopted meshes. The adopted metal ENCUT is 2.0 times the reference value because the required 1.3x check remained above 0.2 meV/atom; higher-cutoff checks are reported in convergence.csv.
The neighboring k20-to-k24 pair changes were -0.564425 meV/atom (Hf), +0.061875 meV/atom (Ti), and +0.077300 meV/atom (Zr). Therefore the strict <0.2 meV/atom neighboring-mesh criterion was met for Ti and Zr but not Hf within the computed mesh range; Hf k24 is retained as the highest tested mesh and this deviation is explicit.
Missing convergence levels are retained as n/a:missing_output in convergence.csv rather than being inferred from another k-point or ENCUT.
All nine relaxed structures retain their target space group according to spglib.
The actual VASP executable reported revision 5.4.4.18Apr17-6-g9f103f2a35, build Feb 06 2024. Reference archives reported the same executable revision family but different build metadata for Hf and SiO2; this is recorded in reference_settings.csv and findings.meta.json.
