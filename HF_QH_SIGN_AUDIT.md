# Hf historical QH sign and convention audit

## Result

No bcc/hcp subtraction reversal was found from the phase-resolved source exports through the plotted pairs. All 1252 temperatures were checked. Both reference and QH use bcc minus hcp. The target is ΔR = ΔG_ref − ΔF_QH. Residual cancellation is real in these data. It is not evidence that the QH curve is physically valid.

At 1400 K (meV/atom):

| Quantity | Value |
| --- | ---: |
| Reference ΔG | 17.600000 |
| ΔF_QH | 297.428859 |
| Target ΔR | -279.828859 |
| ΔE_static | 152.032509 |
| ΔF_vib | 145.396350 |

The static bcc−hcp gap is about 152.03 meV/atom; the thermal part adds another positive gap. A sign flip is not a legal repair. Across this full support, baseline gap MAE is 342.364 meV/atom; reversing it gives 358.733 meV/atom. These are descriptive sign checks, not held-out model scores.

## Four checks

1. Independently subtract the two original QH exports. They match canonical QH and the plotted pair to numerical precision. No phase-specific sign transform occurs in materialization.
2. Independently import the original phase-specific reference tables (meV/atom × 0.001). They match canonical reference values at the common temperatures. The import applies no phase-specific additive alignment. The reference sources and structure metadata identify bcc and hcp separately. Structure metadata reports Im-3m (one-atom primitive cell) and P6_3/mmc (two atoms), respectively. This audit does not rerun symmetry finding or the original phonon force calculations.
3. Verify ΔR = ΔG_ref − ΔF_QH and ΔG_pred = ΔF_QH + Δcorrection. Maximum absolute check errors are in `data/hf_sign_audit/AUDIT.json`; all are below 1e-12 eV/atom.
4. Inspect the baseline construction. The actual canonical arrays match the early 2×2×2, three-volume prototype under legacy support, not the later enlarged-supercell run. F_QH = E_static + F_vib, per atom. Static energies are −7.0414142013 (bcc) and −7.1934467107 (hcp) eV/atom, constant along this selected-volume curve. Both select the upper sampled volume scale 1.02 throughout the temperature support. This is a discrete-volume diagnostic, not a converged EOS minimization.

## Energy content and limits

The archived QH runner adds the potential static energy and Phonopy thermal free energy. The harmonic thermal convention includes zero-point energy; this audit checks the implementation path but does not recompute its phonon integral or isolate a numerical ZPE from the high-temperature table. There is no explicit PV term or pressure-dependent minimization in this runner; materialization inherits the reference P = 0.0001 GPa column. At the stored volumes PV is approximately 0.014 meV/atom, far below the observed hundreds of meV/atom gap. Electronic and magnetic thermal terms are not explicitly included in this baseline.

The reference and model have different absolute energy zeros. The historical reconstruction applies one train-only scalar offset, common to both phases; it cancels in ΔG. No additional phase-specific offset is introduced here. Common-offset invariance cannot establish that the model static phase gap is physically accurate, or independently certify the original reference publisher's cross-phase energy convention.

The bcc legacy phonon report has negative frequencies, including about −1.963 THz at the selected volume. It explicitly marks this baseline diagnostic_only=true and qh_reliable=false. Ordinary harmonic thermodynamics is not valid for an imaginary oscillator. The exact historical mode filtering/weighting needs its original runtime or a mode-resolved cache to be reconstructed; it cannot be inferred from this aggregate CSV alone. Therefore residuals here must not be called pure anharmonic, electronic or magnetic corrections.

## Interpretation

The evidence supports a large mismatch between a diagnostic numerical baseline and the reference, plus failure of this historical learned correction to reproduce the needed cancellation. It does not identify which physical contribution causes the mismatch. The baseline generator uses a DPA alloy head; the separate historical residual-model initializer had a Si/O mapping problem. These are distinct checkpoints and should not be conflated.

[Four-panel audit figure](figures/hf_qh_sign_audit.pdf) · [All paired values](figures/source_data/hf_qh_sign_audit.csv) · [Checks and input hashes](data/hf_sign_audit/AUDIT.json).
