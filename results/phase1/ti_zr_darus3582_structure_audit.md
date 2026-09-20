# DaRUS 10.18419/darus-3582 Ti/Zr structure audit

Audit date: 2026-09-21.  The audit was performed read-only on the dedicated
`thu-GenSi` checkout under
`/GenSIvePFS/users/zirenj/fes_experiment_pbe_d3bj/external/darus_ti_zr_hf/`.

The directory contains exactly these relevant source artifacts:

| Material | Files found | Structure-bearing content |
|---|---|---|
| Ti | `Ti_bcc.tab`, `Ti_hcp.tab` | none |
| Zr | `Zr_bcc.tab`, `Zr_hcp.tab` | none |
| Hf (comparison) | `Hf_bcc_PBE.tar.gz`, `Hf_hcp_PBE.tar.gz` | DFT `OUTCAR` files inside the archives |

The Ti/Zr `.tab` files are quoted thermodynamic tables.  Their headers are
`T (K)`, `Gibbs_energy`, `Enthalpy`, `Entropy`, `Cv`, and `V^eq` (with the
corresponding phase labels in the filenames).  They contain temperature,
free-energy, and equilibrium-volume columns only.  No `POSCAR`, `CONTCAR`,
`OUTCAR`, coordinate block, lattice matrix, or structure archive is present in
the Ti/Zr source or in the copied `data/raw/ti` and `data/raw/zr` directories.

Therefore the missing input is the archive/structure itself, not a parser
format mismatch.  The existing tables cannot be converted into representative
atomic structures without inventing lattice and coordinate data.  No Ti/Zr
structure was fabricated or substituted, and the corresponding structure-based
QH/phonon task remains blocked pending the original structure files or an
explicitly authorized external source.
