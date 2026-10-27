# Optimization Sanity Check and layer-probing protocol

## Current verdict

Source fit, transfer accuracy and phase selection answer different questions. C-DPA already fits the source relative-G curves closely in the two reported outer folds, yet its target-pair error rises. Historical Hf T2 also has a source-fit problem. These observations belong in different diagnostic categories.

The new audit uses frozen point predictions. It does not train models, refit offsets, tune checkpoints or run layer probes. Seven cases cover the archived matched temperature protocols and the two C outer folds. Each selected case shows complete source-fit and declared test error distributions, with separate seed curves and downloadable points. The full summary also retains validation and source-ID statistics.

| Fold | Readout | Source-fit ΔG MAE | Target-pair ΔG MAE | Source contrast vectors | Feature rank |
| --- | --- | ---: | ---: | ---: | ---: |
| SiO₂ cristobalite LOPO | C-DPA | 0.001155 | 1.181740 | 1 | 1 |
| SiO₂ cristobalite LOPO | C-simple | 0.001155 | 0.439594 | 1 | 1 |
| Ti-family LOSO | C-DPA | 0.034375 | 9.139169 | 2 | 2 |
| Ti-family LOSO | C-simple | 7.747816 | 7.026457 | 2 | 1 |

Errors are meV/atom. Within each fold, the source-fit statistic averages permitted source pair–temperature points; the outer score uses target-related pairs only. The source and target supports differ by the task definition. The table is not an average over independent materials.

C-DPA is a source-fit-success / target-degradation case in this audit. C-simple's metal source-fit error shows that its restricted readout has a different fit ceiling. Its ridge solver is a deterministic linear solve; a large residual alone does not diagnose a numerical optimization failure. Fit diagnostics must distinguish readout capacity, feature geometry and numerical optimization.

## Distribution definition and completeness

The plotted quantity is |ΔG_pred − ΔG_ref| at each available paired state. Tables give mean, median, P90, maximum and signed bias. The empirical CDF uses all paired points. A seed remains one curve; temperatures, pairs, domains and seeds are not additional materials. There is no invented convergence threshold.

Direct G and QH-residual models are evaluated after their recorded reconstruction and train-only calibration. A shared absolute-G gauge cancels in ΔG. For delta-G-trained models, phase-relative fit is the primary sanity check. This does not certify unconstrained absolute G. Deterministic QH, zero-gap and source-family controls have source-domain errors but no optimized training loss. Framework/task combinations without matching point predictions remain unassessed; no missing curve is filled with zero.

The loader selects the later explicitly matched SiO₂ predictions in place of older exports, and the multiseed Hf rows in place of overlapping historical rows. It retains the distinct Hf batch-2 and selected-8000 protocols separately. In the C outer score it excludes the diagnostic source-pair rows that are not target-related. Exact exclusion counts, input IDs, original hashes and per-model/domain summaries are in `data/optimization_check.json`. Duplicate canonical pair-state keys are rejected. Frozen original files remain unchanged.

## Update after the independent transfer audit

The Ti family-control outer pair direction is corrected in the active source-fit/test audit. Its revised test MAE is 15.8708 meV/atom. Original frozen predictions remain traceable.

For the current SiO₂ fold, one source contrast plus a linear ridge readout forces a rank-one coefficient map in exact arithmetic. A layer probe can change target scale factors but cannot generate an independent target temperature curve under this same design. Therefore prioritize independent admitted source contrasts before using layers to make a claim about independent thermodynamic coefficients. The saved DPA map is nearly rank one and C-simple is rank one; the observed shared roots are quantified in [the transfer audit](TRANSFER_P0_AUDIT.md).

## Layer probes: scientific question

Which frozen DPA stage supports source-fitted and transferable prediction of phase-relative temperature coefficients, gaps and slopes?

A failed readout is a result about that feature/readout/data combination. It is not a proof that the representation contains no thermodynamic information. This distinction is especially important here: the current source designs supply one and two contrast vectors for 128 DPA coordinates. Repeating a phase across temperature adds constraints on its four temperature coefficients; it does not create new structural feature vectors.

## Fixed protocol

1. **Feature provenance.** Reuse the archived outer manifests. Record checkpoint hash, branch/head, element mapping, exact stage names, pooling, geometry and any volume dependence. Identify extractable stages from the actual network before specifying layer counts. Extract embeddings without thermodynamic target labels. Do not treat the residual-model initializer and the alloy checkpoint used to generate a baseline as the same model.
2. **Coefficient supervision.** Use only source-fit temperatures to fit each phase with `G(T) = a + b x + c x ln(x) + d x²`, `x = T / 1000 K`. Supervise phase-level readouts through within-system coefficient contrasts, using one declared gauge. Reconstruct phase differences from phase outputs to preserve cycle consistency. Never fit target coefficients to train the probe.
3. **Slope target.** At a fixed pressure, the derivative implied by this basis is `dG/dT = [b + c(ln x + 1) + 2d x] / (1000 K)`. Evaluate the relative slope from phase contrasts and report its units separately. Keep this analytic-form derivative separate from the archived local-polynomial derivative diagnostic. Target reference curves are read only by the frozen scorer.
4. **Matched probes.** For each declared stage compare linear ridge readouts with C-simple, dimension-matched fixed random features, and the archived final descriptor. Use identical source scaling, coefficient metric, regularization grid and source-validation rule. Also report a matched low-dimensional projection with its transformation fitted on source structures only. Report feature width, rank and singular values; do not infer capacity from nominal width alone.
5. **Fit checks.** Report source-fit curve residual distributions, coefficient-fit error in the curve-induced metric, source-validation error and numerical solve residual. If a nonlinear probe is added, include a small source memorization check and learning curves before interpreting transfer. Use the same capacity budget at each stage.
6. **Negative control.** Permute coefficient targets across source structures, preserving complete coefficient vectors and the declared gauge. Do not shuffle temperatures independently. With only two source phases the permutation control has very few distinct assignments; report the actual count.
7. **Frozen evaluation.** Save all declared stage/readout configurations before scoring the outer target. Report ΔG, relative slope, stable-envelope accuracy and event metrics where reference events exist. Retain all stage outcomes rather than selecting a target-best layer.
8. **Scope of model selection.** Current source-temperature validation can select a ridge penalty; it cannot validate a choice of representation for unseen structures when there is only one source contrast. The existing outer folds have already been inspected, so new probes on them are exploratory. Reserve untouched phase/system folds for confirmatory selection and comparison when compatible references are admitted.

## Interpretation rules

| Observed pattern | Supported conclusion | Next discriminating test |
| --- | --- | --- |
| Source fit is poor at all stages | This probe setup does not reproduce source curves | Check coefficient scaling, solve residual, readout capacity and source feature rank |
| Source fit is accurate; held-out error rises | This fitted mapping loses accuracy in the declared transfer task | Compare stage, simple descriptors and dimension-matched random controls |
| One stage improves unseen gaps and slopes consistently | That stage supports a better readout under the fixed protocol | Confirm on untouched folds and matched feature budgets |
| DPA and random/simple controls behave alike | No demonstrated value from DPA features under this data budget | Add independent source structures and predeclared probes |

## Run status

- Frozen-prediction fit audit: completed; no new model training.
- Source design rank inspection: completed from the saved readout matrices.
- Intermediate-layer feature extraction and layer probes: **not run**.
- Additional untouched transfer folds: **not yet scored**.

Use `postprocess/build_optimization_check.py` for the current distribution audit. The layer experiment above is a protocol, not a reported result.
