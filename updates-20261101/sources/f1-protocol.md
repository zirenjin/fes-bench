# Hessian fine-tuning test: does fitting displaced-supercell DFT forces improve held-out phase-relative free energies?

Frozen before any fine-tuned model is trained or evaluated. Arms F0 (no fine-tuning) and F1 (forces-only fine-tuning) only. F2 and F3 are deferred and depend on the F1 outcome.

## Data and split

- **Source.** Alexandria PBE phonon files already on thu-GenSi: 1,850 entries in `fes-alexandria-polymorph-audit-20261030/data/raw`, CC BY 4.0. Each displaced supercell carries DFT forces, and the training frames are taken from those. No new DFT.
- **Test formulas.** A random half of the 239 primary-selection formulas (`metadata/selection/FINAL_GROUPS.csv`): 120 formulas, `numpy.random.default_rng(20261101).permutation` of the sorted formula list, first 120.
- **Training pool.** Every downloaded entry whose formula is not a test formula. This includes entries the DFT-side selection removed, and formulas outside the primary selection.
  - No entry of a test formula is used in training, validation or model selection.
- **Validation.** A random 10% of the training-pool formulas (same generator, after the test draw), grouped by formula. It is used only to monitor training.

## Model and training (F1)

- **Starting point.** DPA-3.1-3M, branch `MP_traj_v024_alldata_mixu`, via `dp --pt train --finetune <ckpt> --model-branch MP_traj_v024_alldata_mixu --use-pretrain-script`.
- **Loss.** Forces only (`start_pref_f = limit_pref_f = 1`, energy and virial prefactors 0). No energy labels are supplied.
  - Same-composition ΔE₀ is not a training target. It is reported, to detect damage to the anchor.
- **Schedule.** Fixed a priori: 20,000 steps, `batch_size = "auto:512"`, exponential learning-rate decay from 1e-4 to 3e-6. Seed 20261101.
- **One fallback only.** If validation force RMSE ends above its step-0 value, there is one rerun at start lr 3e-5, decided on validation alone. No other tuning.

## Evaluation

Run `mlip_vs_dft_phonons.py`, unchanged except that it loads a single-task model, on every primary-selection entry of the test formulas. F0 values are the existing primary-head results for the same entries.

**Primary endpoint.**
- Quantity: RMS of the pair thermal error δΔF_th(1000 K) over the test-formula pairs, F1 against F0.
- Uncertainty: ratio with a 95% cluster bootstrap over test formulas, 2,000 draws.
- **Success requires both:** a point reduction of at least 20%, and an upper 95% bound on the ratio below 1.

**Secondary endpoints** (all on test formulas):
- pair δΔF_th RMS at 300, 600 and 1500 K;
- pair total-ΔF RMS;
- pair ΔE₀ RMS;
- single-phase δF_th RMS, which tests whether single-phase gains survive in pair differences;
- non-degenerate wrong-ordering rate (|ΔF_DFT| ≥ 1 meV/atom);
- coverage at α = 0.05 under the held-out-formula abstention rule (thresholds leave one formula out within the test set);
- dense-mesh minimum frequency, as an imaginary-mode count.

**Sanity checks** (not evidence):
- the same metrics on training-pool primary formulas (in-sample);
- validation force RMSE.

**Stop.** If the primary endpoint fails, the Hessian fine-tuning route is dropped, and F2 and F3 are not run.
