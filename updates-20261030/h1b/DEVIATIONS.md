# Deviations from H1B_PROTOCOL.md

All recorded **before** any DPA output was scored and before any target free-energy label was used. None changes the readout, isobars, temperature window, criteria, structures or checkpoint.

## D1. Gate E2 supercell size (label-free)

- **Frozen text.** Target mean descriptors are compared between the unit cell and a 3×3×3 (coesite) or 4×4×6 (stishovite) supercell, tolerance 1×10⁻⁶.
- **Why it cannot run as written.**
  - Measured peak resident memory of `eval_descriptor` on CPU is about 25–27 MB per atom (cristobalite 192 atoms: 5.22 GB; tridymite 216 atoms: 5.78 GB).
  - The frozen supercells have 1,296 and 576 atoms (estimated ~35 GB and ~15 GB). WSL has 15 GB.
  - The first Stage 0 attempt exhausted memory, and this is the probable cause of the earlier WSL service failure.
- **Change.** Use the smallest repetition giving every perpendicular cell width ≥ 12 Å (twice the 6 Å cutoff): coesite 2×1×2 (192 atoms), stishovite 3×3×5 (90 atoms). Tolerance unchanged. Implemented in `scripts/embed_task_checkpoint_v2.py`, which differs from the frozen script only in the supercell lines (diff ignoring line endings).
- **Effect.** E2 still tests the same property: invariance of the mean descriptor to periodic replication beyond the cutoff.

## D2. Gate R1 comparison precision (reference only)

- **Frozen text.** Agreement ≤ 1×10⁻⁶ meV/atom with `round2/SIO2_ASSESSED_TIERS_DS62_SLB.json`.
- **Why.** That JSON stored ΔG rounded to 0.001 meV/atom, so 1×10⁻⁶ agreement is impossible by construction.
- **Change.** Agreement at the stored precision: |difference| ≤ 5×10⁻⁴ meV/atom. Implemented in `scripts/reference_ds62.py`.

## Runtime note

Stage 0 runs under an address-space cap (`ulimit -v`, 20 GB) so that a memory overrun fails the process rather than WSL.

## D3. Stage 0 run as one process per structure (implementation only)

- **Why.** In a single process, memory accumulated across structures: resident memory reached 15.0 GB, about all of WSL's memory, and the run failed. Quartz alone peaks at 7.8 GB.
- **Change.** `scripts/embed_one_save.py` embeds one structure per process. Each target's supercell, built by the D1 rule, runs in its own process too. `scripts/assemble_stage0.py` then performs the same E1, E2 and OOD computations as `embed_task_checkpoint_v2.py` and writes the same files.
- **Effect.** None on numerical definitions. Recorded before scoring; no target label had been used.

## D4. E2 evaluated per structure set; sensitivity set blocked

- **Why.** The D1 supercell of the 9.26 GPa stishovite (3×3×5, 270 atoms) exceeded the 18 GB address cap. Dense stishovite needs about 50 MB per atom; the ambient stishovite supercell already peaked at 13.1 GB resident. Raising the cap risks crashing WSL.
- **Change.** E2 is evaluated separately for the primary set (coesite_ambient, stishovite_ambient) and the sensitivity set (coesite_2.84GPa, stishovite_9.26GPa). A structure whose supercell could not run fails E2 for its set. Scoring requires E2 for the set being scored, so the sensitivity set is blocked and will not be scored. The verdict-governing primary set is unaffected.
- Recorded before scoring; no target label had been used.
