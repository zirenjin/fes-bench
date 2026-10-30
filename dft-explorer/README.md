# DFT → Free Energy Explorer

This browser-only page visualizes existing FES-Bench computations. It submits no DFT, MD, model inference or training jobs. The source is the project compute ledger, independently hashed local VASP/phonon archives, and a timestamped read-only Expanse refresh when included.

## Evidence layers

1. **Workflow register.** Preserve scientific values, failures, retries, job arrays, model sampling and postprocessing; machine-specific path prefixes are redacted. A register row is not a unique completed DFT calculation. The original archived ledger and later scheduler refresh have separate provenance.
2. **Raw ionic outputs.** Every distinct parseable VASP XML in the two audited local archive roots is exported. Additional completed MgSiO3 pressure-matched calculations are extracted read-only from Expanse. A normal footer is reported independently of label admission. Truncated XML is listed as an extraction issue rather than repaired or fabricated.
3. **Trajectory playback.** CaSiO3 V319 includes the archived merged 1,600-frame trajectory plus its 600-frame continuation as a linked, overlapping segment. V323 is a separate bracket archive. The renderer uses recorded fractional coordinates and cells. Forces/energies bind to a trajectory only if the count matches every ionic frame one-to-one. Force arrows share a run-level scale. Missing SCF/energy series remain missing.
4. **Force constants and Γ modes.** Primitive/supercell geometry comes from archived Phonopy metadata; Ti uses its matched unitcell and perfect-supercell POSCAR. Both compact and full FORCE_CONSTANTS are supported. Compact row labels are supercell representative indices, not assumed consecutive primitive indices. Supercell coordinates map to primitive atoms modulo lattice translations, with a 1e-4 Å tolerance. At Γ, force constants are summed over matching supercell images, mass weighted, symmetrized for diagonalization, and converted with the archived Phonopy factor. Displacements use mass-unweighted eigenvectors and a visual amplitude. Negative frequencies remain negative. No arbitrary motion is substituted for an eigenvector.
5. **Phonon spectrum.** Histograms preserve q multiplicity. The vertical axis is weighted mode count per frequency bin, not a normalized DOS per atom. Raw Γ acoustic residuals remain in this view; this is not a replacement for the reference's non-acoustic stability gate. A Γ mode plot is not a full dispersion path.
6. **Thermodynamics.** Carbon v1 and the extended-volume v2 table remain separate. Ti uses mesh-15 comparison rows; kJ/mol of Ti2 is converted to eV/atom with division by two. Estatic + Fvib is stored data. Adding PV uses 1 GPa·Å³ = 0.006241509074 eV. The highlighted minimum is over sampled volumes only, not a new EOS fit or a new admitted Gibbs label.
7. **Phase comparison.** The carbon sign is Ggraphite − Gdiamond. Only recorded P–T points are evaluated. Error bars show documented numerical sensitivities, not total physical uncertainty. Sparse cells are not a continuous phase boundary.

## Scope boundaries

- The page covers all rows in the available project ledger and all parseable raw outputs in the declared archive roots, with the explicit recent Expanse export added. It is not a claim to scan every unrelated remote directory.
- New MgSiO3 force-grid output is not promoted to an admitted free-energy grid while its force-constant/QHA postprocessing remains pending or failed.
- The Ca V05 `phonon_assembly_v05_retry3` matrix is a historical superseded reconstruction. It can be inspected as a diagnostic, but is not the corrected thermodynamic reference. The later accepted force-response audit and thermal candidate are separate records.
- Independent source literature tables are not counted as DFT jobs run by this project. Fe model sampling and TI preparation are not counted as newly accepted DFT free energies.
- Harmonic/QHA labels do not certify general anharmonic, electronic or magnetic entropy. Fixed-volume Svib is not the pressure-minimized total entropy.
- No licensed pseudopotential bodies, wavefunctions, charge-density grids or credentials are published. A dedicated band path or density isosurface is not invented when the data is absent.

## Interactions

- Select a system and a raw calculation. Drag the structure to rotate it and scroll to zoom. Use the ionic slider or click an energy-trace point. Static runs may have only one ionic frame.
- Select a force-constant dataset. Click a Γ-frequency point or choose a mode. Inspect the force-constant block norms and the predicted force response to a 0.01 Å displacement of supercell atom 1.
- Select a thermodynamic layer, phase and stored temperature. Click a sampled volume to open the nearest available matching phonon dataset. It is explicitly a nearest-volume match, not a claim that two different geometries are identical.
- Change pressure to view the PV term. Use the separately admitted carbon grid for reference phase ordering.
- Filter the workflow and label tables. All matching rows are retained; there is no hidden display sampling.

## Rebuild and validation

`FES_ARCHIVE_PROJECT_ROOT=/path/to/project python3 build_data.py` reads the existing project archives and writes `data/`. The exporter requires NumPy and PyYAML. Optional `data/remote_refresh.json` and `data/remote_raw_export.json` preserve the timestamped remote export without requiring SSH during a rebuild. JSON and matching JavaScript data files support both static hosting and local playback. The source-hash inventory is `data/provenance.json`.

`VALIDATION.json` records numerical and browser checks. The browser test exercises desktop and mobile layouts, raw-frame selection/playback, mode selection, reference-grid pressure changes, table filters and missing-output handling. Screenshots are review artifacts. The plotting runtime is a vendored Plotly bundle with its original license header; the atomic renderer is local canvas code. Interactive plots are report views, not submission-ready scientific figure exports.

Publication adds this subtree and one navigation link to the existing GitHub Pages site. The latest user-edited main website is retained. Integration happens in an isolated publication checkout; the original working directories are preserved.

Public mirrors redact machine-specific path prefixes. Source hashes refer to the original archived files; the public mirror hashes are recorded separately in SHA256SUMS. Numerical values, units and job identities are unchanged.
