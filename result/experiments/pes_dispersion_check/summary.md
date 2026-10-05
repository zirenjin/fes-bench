# SiO2 analytic D3(BJ) dispersion check

D3(BJ) parameters: PBE s6=1, s8=0.7875, s9=0.0, a1=0.4289, a2=4.4407; these match the pairwise VASP IVDW=12 PBE-D3(BJ) convention.
The current-DPA estimate is E_DPA + analytic D3(BJ). DFT-relaxed and archived snapshot energies are already PBE-D3(BJ), so their D3 values are reported as corrections only and are not added again.

| Pair | ΔE(DPA+D3) (eV/atom) | ΔE(DFT final) (eV/atom) | model−DFT (meV/atom) | ΔE(snapshot) (eV/atom) | model−snapshot (meV/atom) |
|---|---:|---:|---:|---:|---:|
| quartz_minus_cristobalite | 0.00348405 | -0.02828616 | 31.770 | 0.03824561 | -34.762 |
| quartz_minus_tridymite | 0.00086930 | -0.02729380 | 28.163 | 0.03305306 | -32.184 |

Expanse job 54600945 (PBE, IVDW=0) was not accessible from this environment; no pure-PBE VASP output was found locally, so the direct PBE versus PBE-D3(BJ) comparison remains pending.
