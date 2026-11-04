# Ti hcp harmonic-QHA comparison layer v1

This versioned data package contains a seven-volume, single-phase Ti hcp harmonic-QHA comparison surface. It exports 217 15x15x15 F(V,T) rows and 155 15x15x15, degree-3 G(T,P) rows for 0–4 GPa and 0–1500 K. All 155 fitted minima lie inside the sampled seven-volume interval.

The phase-resolved static plus harmonic free energies use VASP PBE PAW Ti_pv, quantum harmonic nuclei, fixed c/a under isotropic volume scaling, and Te=0 K. This is a scoped QHA comparison layer, not a hcp–bcc ΔG, not an anharmonic TI truth layer, and not an uncertainty-qualified phase stability claim.

The per-state sensitivity table compares 9/12/15 q meshes and degree-2/3/4 fits. The qmesh12-to-15 difference reaches 0.2176 meV/atom, fit-degree sensitivity reaches 0.5263 meV/atom, and their summed sensitivity envelope reaches 0.7138 meV/atom. These are deterministic sensitivity estimates, not confidence intervals. Supercell/displacement/k convergence and systematic error remain unknown.

Raw inputs, mesh files, and their SHA-256 hashes are listed in `INPUT_PROVENANCE.json`; the full audit is `outputs/fes-bench-reference-v2/w9-ti-hcp-qha-volume-pilot-20261021/TI_HCP_QHA_INDEPENDENT_AUDIT_20261025.md`.

Each F(V,T) row includes its 9³→15³ and 12³→15³ vibrational q-mesh sensitivity. The maximum differences are 0.7993 and 0.2190 meV/atom, respectively; these remain deterministic sensitivities, not confidence intervals.
