# Phase 1 delivery — data preparation and representative audit

Status: partial; no unsupported phase has been promoted into the frozen
benchmark.

## Files changed or added

- `fes_bench/data/import_thermo.py`, `align_grid.py`, and `crossings.py`:
  normalized PBE/PBE-D3(BJ) reference tables and generated common-grid crossing
  files.
- `fes_bench/data/standardize_structure.py` and
  `fes_bench/data/relax_representatives.py`: source-standardization and
  explicitly configured representative-relaxation paths.
- `configs/import/{hf,ti,zr}_*.yaml` and
  `configs/representatives/sio2_pes_*fixed_cell.yaml`.
- `data/{hf,sio2,ti,zr}/`: normalized reference tables, metadata, source
  manifests, and available representative structures.
- `data/sio2/exclusions.json`: C2221 exclusion record.
- `results/phase1/casio3_source_hold.md` and
  `results/phase1/audit_20260919.log`.
- `results/phase1_remote_audit_20260920.log`: a fresh dedicated remote audit
  confirming that the C2221 missing-source condition remains explicit.

## Actual commands and key logs

```bash
PYTHONPATH=. python3 -m fes_bench.data.crossings \
  --config configs/import/sio2_crossings.yaml
PYTHONPATH=. python3 -m fes_bench.data.audit --config configs/audit.yaml \
  > results/phase1/audit_20260919.log 2>&1
```

The crossing command reports quartz–cristobalite at 1542.7087 K.  The audit
log reports Hf hcp–bcc on the normalized common grid and enumerates the
representative provenance.  It intentionally exits nonzero because Ti and Zr
still lack source representative structures; the error is recorded rather
than hidden.

On the dedicated remote host, in the isolated directory
`<remote-repo>/fes-bench`, the SiO2 fixed-cell, SiO2-PES-adapted
relaxations converged at 0.001 eV/Angstrom:

| phase | final max force criterion | volume A3/atom | max displacement A |
| --- | ---: | ---: | ---: |
| quartz_beta | < 0.001 | 12.0000 | 1.4952 |
| cristobalite_beta | < 0.001 | 14.2857 | 1.7802 |
| tridymite_p63mmc | < 0.001 | 14.0000 | 1.4358 |

The command logs are preserved as
`results/phase1/sio2_pes_quartz_sanity.log` and
`results/phase1/sio2_pes_remaining_fixed_cell.log` on the remote host; the
per-phase results are in each `meta.json`.

## Artifacts

- `data/sio2/reference_crossings.json`: three-phase result; C2221 pairs are
  absent.
- `data/hf/reference_crossings.json`: hcp–bcc crossing at 1919 K on its
  quantized table grid.
- `data/sio2/{quartz_beta,cristobalite_beta,tridymite_p63mmc}/structure.extxyz`
  and corresponding `meta.json`: fixed representatives for all later runs.
- `data/hf/{hcp,bcc}/structure.extxyz`: standardized source structures.
- `results/phase1/casio3_source_hold.md`: source limitation with explicit
  evidence.

## Acceptance check

- SiO2 quartz–cristobalite crossing: **pass**, 1542.7087 K (target about
  1542.7 K).
- Hf hcp–bcc crossing: **pass**, 1919 K.
- C2221 uses no COD fallback: **pass**; it is excluded rather than replaced.
- All minimum systems have a validated common-grid reference plus a fixed
  representative: **not yet met**.  Ti/Zr lack source structures, and the
  available CaSiO3 archive lacks absolute per-phase G(T,P) tables.

## Deviations from design

- SiO2 C2221 is excluded from the active system phase list because no DaRUS
  source structure is available.  No COD replacement was used.
- SiO2 representatives use the SiO2-PES-adapted single-head checkpoint with
  a fixed cell, rather than a generic multi-task `Domains_SSE_PBE` full-cell
  relaxation.  The latter changed the structure excessively and is recorded
  only as a failed sanity attempt.  The formal FES representation remains
  `DPA-3.1-3M.pt` with `Domains_Alloy`; it is not used to generate 0-K SiO2
  structures.
- The available DaRUS snapshots are high-temperature/P1-like inputs, not
  low-temperature relaxed structures.  The large fixed-cell displacements and
  symmetry changes are preserved in metadata and prevent calling these
  structural representatives equivalently validated to a DFT 0-K source.
- Ti/Zr have normalized reference tables but no locally available source
  structures for representative construction.
- A 2026-09-19 read-only recheck found only large `T=1000 K` NEP titanium
  snapshots in `FES-data/clean_revision2157`, and no corresponding Zr
  representative. They are neither 0-K/lowest-temperature DaRUS structures
  nor appropriate DPA-relaxation inputs, so they remain excluded rather than
  being silently repurposed.
- CaSiO3 is held: its available Zenodo copy provides relative Delta-G data at
  50 GPa and phase-boundary points, not two absolute per-phase G(T,P) tables.
