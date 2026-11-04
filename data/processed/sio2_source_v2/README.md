# SiO2 four-phase source table layer v2

This is a separate versioned view of the checksum-verified DaRUS-4999 PBE-D3(BJ) direct-upsampled tables. It leaves the historical three-phase system and frozen split unchanged.

## Thermodynamic source admission

The four phase-resolved G tables and six same-state pair differences are **admitted for source-level thermodynamic evaluation** at the published source tier. They use the same composition (SiO2), PBE-D3(BJ) Hamiltonian, PAW/VASP source framework, direct-upsampling protocol, and controlled pressure. Same-source phase differences therefore do not require an independently assigned elemental energy zero. The data contain 1,649 exact common integer-temperature rows from 851 to 2499 K at 0.0001 GPa (0.1 MPa); no G interpolation was used. The published procedure uses phase-specific MTPs and finite-temperature Langevin NVT sampling/thermodynamic integration.

This is source-level admission, not an uncertainty-qualified precision claim. Per-row statistical/numerical uncertainties and pair covariance are not exposed; they are recorded as unknown, never zero. The source table does not separately expose electronic temperature, and quantum nuclear corrections are not separately reported. Do not infer a transition uncertainty or assert a QHA reconstruction from these tables.

## Structure status

`STRUCTURE_IDENTITY_STATUS.csv` tracks geometry separately. No exact rowwise thermodynamic geometry map is available for the four tables. This limits structure-based migration and DFT reproduction; it does not invalidate the published phase-resolved G values or same-source DeltaG values for thermodynamic baseline evaluation. The C2221 source FCP exists at matching PBE-D3(BJ) level, but its rowwise thermodynamic geometry association is still unknown; the other released effective harmonic FCPs are LDA.

## Data and splits

`pair_deltaG_candidates.csv` is retained as a compatibility filename; its row-level `admission_state` now records source-thermodynamic admission and unknown uncertainty. The four tables preserve exact raw source rows and checksums. `splits/phase_lopo_v3.json` is a standalone four-fold leave-one-phase-out split using the exact common T/P nodes. No historical split or canonical master manifest was changed.
