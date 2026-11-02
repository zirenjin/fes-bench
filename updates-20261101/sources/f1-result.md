# F1 result: forces-only fine-tuning fails the primary endpoint (protocol stop)

## Run

- **Protocol.** `PROTOCOL.md`, SHA-256 ec536cdb…, frozen 2026-11-01T20:12Z. `PROTOCOL_ADDENDUM_A1.md`, SHA-256 5a7062a1…, frozen 22:17Z, before any F1 evaluation.
- **Split.** 120 test formulas: 268 entries, 185 pairs. Training pool: 1,313 systems, 50,703 frames. Validation: 131 systems.
- **Training.** 20,000 steps. Validation force RMSE fell from 0.0385 at step 1 to 0.0188 at the end (`results/lcurve.out`), so no fallback was triggered.
- **Operational changes (none analytic).**
  - The input `type_map` was set to the 118-element map of the pretrained model, needed for data loading.
  - `dp` is run through `scripts/dp_safe.py`, which carries the same `torch.load` allowlist the evaluation script uses.
  - The evaluation script loads a single-task model (`head = "-"`).

## Primary endpoint (test formulas)

RMS of pair δΔF_th(1000 K): **8.41 → 15.06 meV/atom. Ratio 1.79, 95% CI 1.23–2.27. FAIL.**

| T | Pair δΔF_th | Pair total | Pair ΔE₀ | Single-phase δF_th | Wrong ordering | Coverage at α = 0.05 |
|---|---|---|---|---|---|---|
| 300 K | 2.15 → 4.16 | 8.78 → 15.24 | 8.62 → 12.76 | 3.24 → 20.28 | 9.2 → 17.8% | 0.45 → 0.03 |
| 1000 K | 8.41 → 15.06 | 11.94 → 23.66 | 8.62 → 12.76 | 14.45 → 92.5 | 11.1 → 23.5% | 0.29 → 0.08 |
| 1500 K | 13.06 → 23.0 | 15.59 → 30.72 | 8.62 → 12.76 | 22.81 → 145.9 | 11.4 → 26.6% | 0.24 → 0.11 |

Values are F0 → F1, in meV/atom except where marked.

- **In-sample is also worse.** On the training-pool formulas (166 pairs) the result is similar: pair δΔF_th at 1000 K goes from 5.40 to 8.93.
- **Addendum A1 at 1000 K.** The r(T₀ = 300 K) ratio is 1.74 (1.17–2.22) and the δΔS ratio 1.72 (1.14–2.20). Temperature dependence is not improved.

## Diagnosis: what the force loss actually fitted (`scripts/diag_forces.py`, four training-pool entries)

| Entry | DFT force RMS | F0 equilibrium residual | F0 response error | F1 equilibrium residual | F1 response error |
|---|---|---|---|---|---|
| mp-10080 | 0.0224 | 0.0317 | 0.0011 | 0.0061 | 0.0125 |
| mp-10219 | 0.0074 | 0.0489 | 0.0003 | 0.0065 | 0.0044 |
| mp-10251 | 0.0063 | 0.0387 | 0.0011 | 0.0100 | 0.0037 |
| mp-10281 | 0.0047 | 0.0085 | 0.0006 | 0.0011 | 0.0025 |

All values are eV/Å.

- **Equilibrium residual.** This is the MLIP force at the undisplaced DFT supercell. The DFT residual is about 0.
- **Response.** This is F(u) − F_MLIP(0), the part that determines the force constants.

What the diagnosis shows:

- **The pretrained model's raw-force error is almost all equilibrium residual.** Its geometry mismatch is 0.009–0.049 eV/Å, while its displacement response is already accurate to 0.3–1.1 meV/Å (about 4–17% of the DFT response).
- **The forces-only loss was therefore minimised mainly by removing the residual.** The residual fell about 5-fold, and the response error grew 4- to 15-fold. Phonons got worse even in sample.
- **Raw-force fine-tuning on small-displacement phonon data does not target the Hessian.** F1 tested raw-force fitting, not curvature fitting.

## Status

Under PROTOCOL.md the primary endpoint failed, so this route stops as specified, and F2 and F3 are not run.

A response-targeted variant would be a new experiment under a new frozen protocol. In it, labels are F_DFT(u) + R₀, where R₀ is the pretrained equilibrium residual, so the loss sees only response error. Running it would be a second attempt on the same test split, and must be reported as such.

**Side observation (descriptive).** The large equilibrium residuals mean the MLIP's own equilibrium geometry differs from DFT's. That geometry difference is the volume/relaxation channel that the fixed-geometry Alexandria tests do not see.
