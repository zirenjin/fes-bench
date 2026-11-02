# F1b: response-targeted fine-tuning (second attempt on the same test split)

Frozen before any F1b model exists. This is a new experiment motivated by the F1 diagnosis (`RESULT.md`). It is **the second fine-tuning attempt evaluated on the PROTOCOL.md test split**, and must be reported as such together with F1.

## Change from F1 (the only analytic change)

Force labels are F_DFT(u) + R₀, where R₀ is the pretrained `MP_traj_v024_alldata_mixu` force on the undisplaced DFT supercell of the same system, computed once before training.

- At the starting model the loss therefore sees only the displacement-response error, F(u) − F(0) against the DFT response, and not the equilibrium residual.
- The DFT residual is taken as about 0, the same assumption as in `scripts/diag_forces.py`.

Everything else is identical to F1:
- split (`data/SPLIT.json`);
- starting checkpoint and branch;
- forces-only loss;
- 20,000 steps, `auto:512`, learning rate 1e-4 → 3e-6, seed 20261101;
- the same single fallback rule (validation force RMSE ending above its step-1 value → one rerun at start learning rate 3e-5).

## Evaluation

Identical to PROTOCOL.md and Addendum A1, with F1b in place of F1, against the same F0.

- **Primary endpoint.** RMS of pair δΔF_th(1000 K) on test formulas, ratio F1b/F0.
- **Success requires both:** a ratio of 0.8 or less, and an upper 95% bound below 1.
- **Also reported:** F1b against F1.

**Diagnostic, not evidence.** `diag_forces.py` on the same four training-pool entries, reporting equilibrium-residual and response errors.

## Stop

If the primary endpoint fails, phonon-force fine-tuning is dropped for this project. No third attempt is made on this split.
