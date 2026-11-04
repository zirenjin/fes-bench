# Frozen Reviewer 3 Report

## Overall assessment

The project addresses an important problem in computational materials science. It asks whether a physically motivated quasi-harmonic baseline can improve learned free-energy extrapolation across temperature, phases, and systems. The benchmark design contains several valuable elements, including frozen split definitions, explicit phase-transition metrics, reliability screening for the quasi-harmonic approximation, and negative controls.

The present evidence does not establish the title-level claim. On temperature extrapolation, the QH-residual model performs worse than the zero baseline in delta-G MAE, while the thermo-form comparator performs substantially better. Learned results are absent for phase and system holdouts. The packet also documents incomplete controlled ablations, uncertainty analysis, provenance, and external baselines. These limitations affect originality, importance, technical confidence, and the breadth of the readership case.

For nonspecialists, the current framing is difficult to follow because the physical hypothesis, benchmark task, modeling variants, reliability restrictions, and transition metrics are not yet connected through a complete narrative. Final manuscript readability cannot be assessed because no full manuscript was supplied.

## Who would be interested and why

Researchers in computational materials science would care about robust prediction of phase free energies beyond fitted temperature ranges.

Developers of machine-learned interatomic potentials would care about whether fixed physical components improve out-of-distribution generalization.

Researchers in lattice dynamics and thermodynamics would care about where quasi-harmonic information helps, fails, or becomes unusable because of unstable modes.

A broader machine-learning readership could care if the study establishes a general lesson about physics-based residual learning. That broader case currently remains weak because only one learned out-of-distribution axis has results, the strongest reported comparator is not yet adequately contextualized, and the cross-system evidence is incomplete.

## Major strengths

- The benchmark distinguishes temperature, phase, and system extrapolation rather than treating all held-out data as equivalent.

- Split hashes and explicit train and test counts improve auditability.

- The project measures transition behavior, missed crossings, ranking, and sign accuracy in addition to aggregate free-energy error.

- The QH reliability rule and single-side diagnostic expose an important applicability boundary.

- The project records negative findings, including the limited no-learning simplification of metal residuals and the poor aggregate performance of QH-only and learned residual variants.

- The systems include elemental metals and silica polymorphs, which gives the eventual study some potential interdisciplinary reach.

## Major Concerns

### R3-M1

- Severity | Major
- Blocking | Yes
- Claim pointer | Title claim and project narrative lines 24 to 54
- Evidence pointer | Temperature extrapolation results, LOPO and LOSO results, task status E5 and E10
- Concern | The central claim that quasi-harmonic baselines make learned free-energy surfaces extrapolate is not supported by the current benchmark. Temperature-extrapolation delta-G MAE is 31.239 meV per atom for QH-residual and 4.521 meV per atom for the zero baseline. Learned LOPO and LOSO rows are not trained.
- Why it matters | The claim covers temperature, phases, and systems, but the learned evidence currently covers only temperature and does not outperform the simplest error baseline on the principal MAE.
- Resolution test | Complete the frozen three-way controlled ablation for all three out-of-distribution axes, report uncertainty and pair coverage, and either demonstrate a reproducible advantage attributable to QH or narrow the title and central claim to the supported result.

### R3-M2

- Severity | Major
- Blocking | Yes
- Claim pointer | Proposed use of QH as a physical baseline and physics-regime stress tests in project narrative lines 38 to 52 and 64 to 122
- Evidence pointer | Phase imaginary-mode fractions, single-side QH diagnostic, and no-learning thermo-form ratios
- Concern | The applicability boundary of QH is central rather than peripheral. Five of the nine phases exceed the stated one-percent imaginary-mode threshold. The no-learning diagnostic shows little simplification for Hf, Ti, and Zr, while mixed one-side QH produces 644.868 meV per atom pair-level error.
- Why it matters | A baseline that is unavailable or unreliable on one side of many target transitions cannot straightforwardly support a general extrapolation claim. The learned anchored-shape results on reliable phases are promising, but they do not resolve the two-phase prediction problem.
- Resolution test | Predefine how unreliable phases are handled, complete the planned regime stress tests, and show pairwise results separately for both-side reliable, one-side reliable, and neither-side reliable cases. The final claim must state the demonstrated applicability domain.

### R3-M3

- Severity | Major
- Blocking | Yes
- Claim pointer | Claimed value of the learned QH-residual formulation
- Evidence pointer | Temperature results for thermo-form, QH-residual, representation tlog, and zero predictors
- Concern | The thermo-form fit achieves 0.401 meV per atom delta-G MAE, 0.9951 sign accuracy, and 2.796 K transition error, substantially outperforming both learned variants. The packet does not provide enough information to determine whether this comparison uses an equally constrained training protocol or what scientific function remains uniquely served by the learned formulation.
- Why it matters | The strongest available result appears to favor a simpler thermodynamic form. Without a fair and clearly explained comparison, the proposed learned contribution and its scientific importance remain unclear.
- Resolution test | Specify the thermo-form inputs, fitting data, degrees of freedom, and leakage controls. Compare all models under matched split and training-information constraints. Explain what capability the learned model provides that the thermo-form model does not, or revise the contribution accordingly.

### R3-M4

- Severity | Major
- Blocking | Yes
- Claim pointer | Reproducible benchmark and repository availability
- Evidence pointer | Notebook, code, and provenance state in the packet
- Concern | The supplied export, tracked notebook, repository state, and generated artifacts do not form one reproducible release. The tracked notebook has no saved outputs and differs materially from the export. Sections A to C are not included in the exported tables or manifest. GitHub acquisition is not pinned before retrieval. The working tree is dirty. QH-residual records show an unknown commit, and referenced checkpoints are absent.
- Why it matters | Readers cannot reconstruct which code and data produced the central numbers. This prevents independent verification and weakens every quantitative claim.
- Resolution test | Produce a clean pinned release containing the exact notebook or script, immutable input hashes, all producers, checkpoints or durable substitutes, generated outputs, and a manifest covering every source used by tables, figures, and appended diagnostics. Execute it independently and archive the resulting logs.

### R3-M5

- Severity | Major
- Blocking | Yes
- Claim pointer | Generalization across systems and seed-based uncertainty reporting
- Evidence pointer | Four-system inventory, three-metal LOSO restriction, seed results, and missing uncertainty analysis
- Concern | The evidence base is too small and heterogeneous to support broad cross-system language. LOSO is restricted to Hf, Ti, and Zr, with the Hf fold additionally confounded by temperature range. The reported QH-residual seed standard deviation is about 5.4e-9 eV per atom, yet the planned seed uncertainty work is incomplete and checkpoint provenance is missing.
- Why it matters | Random-seed variation is not a substitute for uncertainty across material systems. With three closely related elemental-metal systems in the LOSO task, the effective inferential unit count is limited.
- Resolution test | Report per-system and per-pair effects, distinguish optimization variability from system-level variability, complete the seed audit, and qualify population-level claims. Resolve the Hf temperature confounding through the overlap-temperature analysis and clearly separate that result from unrestricted LOSO.

### R3-M6

- Severity | Major
- Blocking | Yes
- Claim pointer | Cross-system comparison and aggregate benchmark conclusions
- Evidence pointer | PBE thermodynamic-integration references for metals, PBE-D3(BJ) direct upsampling for SiO2, reference standard deviations, and unstarted functional or truth sensitivity
- Concern | The reference construction differs between metals and SiO2, and reference variability spans more than two orders of magnitude. Functional and truth sensitivity has not been tested.
- Why it matters | Aggregate performance may reflect differences in reference method, energy scale, or system family rather than transferable model behavior.
- Resolution test | Present stratified results by system family and reference-generation route, justify any aggregation, and complete a sensitivity analysis sufficient to show which conclusions are stable to the reference definition.

### R3-M7

- Severity | Major
- Blocking | Yes
- Claim pointer | Accurate phase-transition prediction
- Evidence pointer | Crossing metrics, transition slopes, missed-crossing counts, and SiO2 slope conversion
- Concern | Transition-temperature claims are especially sensitive for the two SiO2 pairs, where 1 meV per atom corresponds to about 132 to 134 K. The packet reports mean numeric transition errors together with aggregate missed crossings, but it does not establish how missing crossings enter the summary or provide complete per-pair and per-seed coverage.
- Why it matters | A low error among surviving numeric crossings can coexist with poor overall transition recovery. Aggregate means can therefore mislead nonspecialist readers.
- Resolution test | Report every reference crossing by pair and seed, including predicted crossing, error, false or missed status, and denominator. Define aggregation before analysis and present coverage beside every transition-error summary.

### R3-M8

- Severity | Major
- Blocking | No
- Claim pointer | Interdisciplinary readership and title-level framing
- Evidence pointer | Missing full introduction, methods narrative, figures, legends, prior-work comparison, and discussion
- Concern | The project is presently written for readers already familiar with QH, TI, TDEP, DPA, LOPO, LOSO, residual learning, and crossing metrics. The conceptual sequence from physical decomposition to benchmark prediction task is not yet available in a coherent manuscript.
- Why it matters | Readers outside computational free-energy modeling will struggle to understand what is predicted, what information is available at training time, why each split is scientifically meaningful, and why the reported error scale matters.
- Resolution test | Add a self-contained conceptual account with defined terms, a simple schematic of inputs and outputs, an explanation of each extrapolation axis, and physical interpretation of error magnitudes. Test the abstract, opening, primary figure, and conclusion with readers outside the immediate specialty.

## Minor Comments

### R3-m1

- Severity | Minor
- Blocking | No
- Claim pointer | Model and benchmark terminology
- Evidence pointer | Definitions of QH-only, QH-residual, representation tlog, and thermo-form
- Concern | The relation among these predictors is not immediately transparent.
- Why it matters | Readers may attribute performance differences to the wrong component.
- Resolution test | Provide a compact model table listing inputs, trainable terms, fitting data, system constants, and intended control function.

### R3-m2

- Severity | Minor
- Blocking | No
- Claim pointer | Skill and accuracy metrics
- Evidence pointer | Metric definitions and central results
- Concern | Negative skill values and the low QH-residual ranking accuracy require interpretation. Metric directionality and chance levels are not evident from the packet.
- Why it matters | Nonspecialists cannot judge whether a reported number represents minor degradation or severe failure.
- Resolution test | State the range, optimum, baseline, and physical meaning for each principal metric.

### R3-m3

- Severity | Minor
- Blocking | No
- Claim pointer | Reported energy and slope values
- Evidence pointer | Results in eV, meV, micro-eV, and kelvin
- Concern | Frequent unit changes make comparisons unnecessarily demanding.
- Why it matters | Conversion burden increases the risk of misreading the practical scale of errors.
- Resolution test | Use one primary energy unit throughout and place explicit conversions beside slope-based temperature interpretations.

### R3-m4

- Severity | Minor
- Blocking | No
- Claim pointer | Inventory and split descriptions
- Evidence pointer | Four systems, nine phases, six within-system pairs, and three extrapolation axes
- Concern | System, phase, pair, frame, fold, and crossing are distinct units but are easy to conflate.
- Why it matters | The sample size and level of generalization can appear larger than they are.
- Resolution test | Include a glossary and one diagram showing the nesting of frames, phases, pairs, systems, folds, and crossings.

### R3-m5

- Severity | Minor
- Blocking | No
- Claim pointer | One-percent imaginary-mode rule
- Evidence pointer | QH reliability classifications
- Concern | The rationale and sensitivity of the one-percent threshold are not supplied.
- Why it matters | The threshold determines which phases can participate in the proposed physical baseline.
- Resolution test | Justify the threshold and show whether reasonable alternative thresholds change the classification or conclusions.

### R3-m6

- Severity | Minor
- Blocking | No
- Claim pointer | Availability and completion status
- Evidence pointer | Task list and rows marked not trained, absent, ongoing, or unstarted
- Concern | Several distinct states of missingness are mixed across the project record.
- Why it matters | Readers need to distinguish an inapplicable analysis from an unfinished run, unavailable artifact, or negative result.
- Resolution test | Use a controlled status vocabulary and preserve those distinctions in tables, figures, and repository manifests.

## Technical failings that need to be addressed before the case is established

- Complete learned evaluation on temperature, phase, and system holdouts under a frozen controlled ablation. R3-M1

- Establish and test the valid QH applicability regime for both sides of each transition. R3-M2

- Resolve the scientific and procedural comparison with the much stronger thermo-form result. R3-M3

- Release one clean, pinned, executable evidence package with complete provenance. R3-M4

- Separate seed variability from uncertainty across systems and resolve LOSO confounding. R3-M5

- Quantify sensitivity to heterogeneous reference construction and truth definition. R3-M6

- Report transition recovery with full per-pair coverage and explicit treatment of missed crossings. R3-M7

## Assessment against Nature-style criteria

| Criterion | Assessment | Basis |
|---|---|---|
| Originality | Not assessable | No complete prior-work comparison or evidence supporting the novelty wording was supplied. The physics-residual framing may be interesting, but its distinct contribution cannot be judged from the packet alone. |
| Scientific importance | Weak | Reliable phase-transition extrapolation would be important, but the current learned method underperforms the zero baseline on temperature MAE and lacks learned phase and system holdout results. |
| Interdisciplinary readership | Weak | The question could connect thermodynamics, lattice dynamics, and machine learning. Current evidence and framing remain too specialized and incomplete to establish broad interest. |
| Technical soundness | Weak | Frozen splits, multiple metrics, and negative controls are strengths. Missing OOD evaluations, uncertainty, sensitivity analysis, and reproducible provenance prevent the central case from being established. |
| Readability for nonspecialists | Not assessable | No full manuscript, primary figures with legends, or complete discussion was supplied. The available material uses dense specialist terminology without a complete explanatory arc. |

## Recommendation posture

The package is promising as a benchmark under development, particularly because it exposes negative results and physically meaningful failure regimes. The current evidence supports continued investigation, not the broad title-level claim. The next review should occur only after the controlled three-axis evaluation, comparator clarification, transition-coverage accounting, uncertainty analysis, and pinned reproducible release are available.

This is not an editorial decision. The review covered the full supplied packet and all five requested criteria. It did not assess unavailable manuscript sections, independently rerun computations, inspect external literature, or verify artifacts outside the packet. No missing evidence was treated as a negative result, and no unavailable claim was inferred to be resolved.
