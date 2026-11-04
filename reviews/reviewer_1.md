# Reviewer 1

## Overall assessment

This is a potentially valuable benchmark project with unusually candid negative controls, explicit split definitions, crossing-aware metrics, and a documented QH reliability rule. The current evidence does not establish the title claim.

For temperature extrapolation, QH-residual has a delta-G MAE of 31.239 meV/atom, compared with 31.557 for representation-only and 4.521 for the zero baseline. Its sign accuracy is higher than representation-only, but its ranking accuracy, numeric transition-temperature error, and missed-crossing count are worse. QH-only performs substantially worse still. Learned results are absent for phase and system extrapolation. The evidence therefore supports the usefulness of the benchmark design and the strength of the thermo-form comparator, but not the claimed general benefit of a QH baseline across temperature, phases, and systems.

This is a bounded assessment of a proposal and partial evidence package. Several conclusions remain not assessable because the packet lacks complete methods, executable frozen outputs, uncertainty analysis, and prior-work comparison.

## Who would be interested and why

Researchers in machine-learned interatomic potentials, computational thermodynamics, phase-transition prediction, and out-of-distribution generalization would care about a benchmark that tests free-energy models through crossings rather than interpolation error alone. The distinction between a physically motivated baseline that helps, one that fails because its regime is invalid, and a compact thermo-form that dominates learned models could be broadly informative. That interest currently attaches more strongly to the benchmark and its negative findings than to the proposed QH-residual method.

## Major strengths

- The project separates temperature, phase, and system extrapolation and records split hashes.
- The metric panel includes free-energy error, sign and ranking measures, crossing errors, missed crossings, and coverage.
- The QH reliability rule is explicit and exposes the severe imaginary-mode problem in several phases.
- The single-side diagnostic and no-learning control directly test whether QH simplifies the target before attributing gains to learning.
- The project records incomplete tasks and negative results rather than concealing them.
- The thermo-form result provides an important challenge to claims that a learned residual is necessary.

## Major Concerns

### R1-M1

- **Severity** Major
- **Blocking** Yes
- **Axis** Technical soundness, scientific importance
- **Claim pointer** The title claim and the stated model progression from a QH-residual construction to the broader cheap-physics decomposition.
- **Evidence pointer** Central current results for temperature extrapolation and the single-side QH diagnostic.
- **Concern** The available results do not show that adding QH makes the learned surface extrapolate better. QH-residual improves delta-G MAE over representation-only by only 0.318 meV/atom, while both are much worse than the zero baseline. QH-residual also has poorer ranking accuracy, larger numeric transition-temperature error, and more aggregate missed crossings than representation-only. QH-only is markedly worse. The evidence is mixed even on the one completed OOD axis.
- **Why it matters** The central claim requires a consistent advantage attributable to the QH baseline on scientifically relevant outcomes. A small MAE difference accompanied by worse crossing behavior does not establish that advantage.
- **Resolution test** Complete a controlled comparison in which QH-residual and representation-only use matched data, architecture, capacity, optimization, selection rules, and seed handling. Report per-pair and aggregate results for every prespecified metric. The claim should be narrowed or reformulated unless QH provides a robust advantage on the metrics that support phase stability and transitions.

### R1-M2

- **Severity** Major
- **Blocking** Yes
- **Axis** Technical soundness, interdisciplinary readership
- **Claim pointer** The title language covering temperature, phases, and systems.
- **Evidence pointer** LOPO and LOSO tables, project task status, and missing-material list.
- **Concern** Learned QH-residual and representation-only results are absent for LOPO and LOSO. The three-way ablation and cross-system learned evaluation are incomplete. The family LOSO setting contains only Hf, Ti, and Zr, and the Hf fold is confounded by temperature range.
- **Why it matters** A claim spanning three OOD axes cannot be supported by learned-model results from temperature extrapolation alone. The small and chemically narrow system set also limits the meaning of across-system generalization.
- **Resolution test** Run the frozen learned-model protocol across all nine LOPO folds and all defined LOSO folds, including the Hf overlap-temperature analysis. Report failures and coverage as outcomes. Treat systems as inferential units and state the generalization domain supported by four systems, rather than implying unrestricted cross-system transfer.

### R1-M3

- **Severity** Major
- **Blocking** Yes
- **Axis** Technical soundness
- **Claim pointer** The proposed physical interpretation that QH supplies a useful cheap baseline whose residual is easier to learn.
- **Evidence pointer** QH reliability inventory, project narrative discussion of missing entropy mechanisms, and single-side QH results.
- **Concern** Five of nine phases fail the stated QH reliability rule. These include all three metal bcc phases and two SiO2 phases. The no-learning diagnostic shows clear simplification for beta-quartz but little simplification for the reliable metal hcp phases. Mixed one-side QH produces 644.868 meV/atom pair-level error, while the project itself acknowledges missing entropy mechanisms.
- **Why it matters** The proposed mechanism is not portable across the benchmark domain if the baseline is physically unreliable on one side of many target transitions. A learner may then be correcting baseline failure rather than learning a simpler anharmonic residual.
- **Resolution test** Predefine how unreliable phases are gated, replaced, or modeled. Demonstrate residual simplification separately for every phase and pair under that rule. Compare matched models with and without each cheap-physics component, and show that any gain persists when one or both sides fail the QH reliability criterion.

### R1-M4

- **Severity** Major
- **Blocking** Yes
- **Axis** Technical soundness
- **Claim pointer** The use of frozen DPA energies and heterogeneous reference free energies as inputs and targets.
- **Evidence pointer** Data and physical domain, frozen DPA pairwise static-energy errors, and the missing fine-tuned PES and functional sensitivity studies.
- **Concern** Frozen DPA static-energy errors of 16.3 to 38.3 meV/atom are large relative to several reported reference delta-G variations. The packet does not show how those errors propagate into pairwise free energies or crossings. The fine-tuned PES acceptance file is absent. Metals and SiO2 also use different reference constructions, while functional and truth sensitivity remain untested.
- **Why it matters** Representation error or reference-construction differences may dominate the residual being modeled. Without an error budget, apparent extrapolation behavior cannot be attributed cleanly to the QH baseline or residual learner.
- **Resolution test** Provide a pair-resolved error budget that propagates static-energy and reference uncertainties into delta-G and transition predictions. Apply predefined PES acceptance criteria. Repeat key conclusions under the planned truth and functional sensitivity analyses, or restrict claims to the verified reference construction.

### R1-M5

- **Severity** Major
- **Blocking** Yes
- **Axis** Technical soundness, readability
- **Claim pointer** The numerical tables and diagnostics presented as repository-derived results.
- **Evidence pointer** Notebook, code, and provenance state.
- **Concern** The numerical package is not presently reproducible as a frozen artifact. The tracked notebook has no execution outputs and differs materially from the supplied export. Sections A to C are appended after the export step. Section C lacks its intended producer. GitHub mode is not commit-pinned. Some processed inputs are unhashed. The working tree is modified, QH-residual provenance records an unknown commit, and referenced checkpoints are absent.
- **Why it matters** Independent verification cannot determine which source, data state, and checkpoints produced the reported values. Syntax checks and internal check counts do not substitute for reproducing the complete analysis from a clean state.
- **Resolution test** Release one clean, commit-pinned package containing synchronized source and notebook files, all producer scripts, exact inputs, complete hashes, checkpoints, environment specification, executed outputs, and a regenerated manifest. A clean rerun should reproduce every reported table and figure.

### R1-M6

- **Severity** Major
- **Blocking** Yes
- **Axis** Technical soundness
- **Claim pointer** Comparisons among zero, interpolation constant, global mean, Bartel, QH-only, learned models, and thermo-form.
- **Evidence pointer** Benchmark design, task status showing the leakage audit is ongoing, and central current results.
- **Concern** Comparative validity is not yet established. The leakage audit remains incomplete, and the packet does not provide enough protocol detail to determine whether the dominant thermo-form fit and learned models use equivalent train-only fitting, model-selection, and tuning constraints.
- **Why it matters** The thermo-form fit outperforms all learned models by a large margin. Its interpretation depends on whether every fitted quantity, hyperparameter choice, and stopping rule was restricted to training information under the same split.
- **Resolution test** Complete the leakage audit and publish model-by-model provenance for feature construction, fitted constants, hyperparameter selection, early stopping, and checkpoint choice. Recompute the comparison after confirming that test temperatures, phases, systems, crossings, and labels played no role in selection.

### R1-M7

- **Severity** Major
- **Blocking** No
- **Axis** Originality, scientific importance
- **Claim pointer** The novelty and scope implied by the title.
- **Evidence pointer** Missing-material list stating that no complete prior-work comparison or support for the novelty wording was supplied.
- **Concern** Originality relative to existing free-energy models, physical residual learning, and phase-transition benchmarks is not assessable from the packet.
- **Why it matters** A technically sound benchmark can be valuable, but Nature-style originality and broad scientific importance require a clear distinction from prior approaches.
- **Resolution test** Supply a claim-specific prior-work comparison showing what capability, benchmark axis, or mechanistic result is genuinely new. Align the title and significance claims with that demonstrated distinction.

## Minor Comments

### R1-m1

- **Severity** Minor
- **Axis** Technical soundness, readability
- **Affected element** Transition-temperature reporting
- **Evidence pointer** Central current results and benchmark metrics.
- **Issue** Mean numeric transition-temperature error is reported alongside aggregate missed crossings, but the packet does not give the denominator and pair-by-seed coverage underlying each mean.
- **Required correction** Report transition errors with explicit denominators, per-pair and per-seed coverage, false crossings, missed crossings, and a stated rule for models that predict no numeric crossing.

### R1-m2

- **Severity** Minor
- **Axis** Technical soundness
- **Affected element** Seed uncertainty
- **Evidence pointer** QH-residual provenance and seed standard deviation.
- **Issue** The reported seed standard deviation of about 5.4e-9 eV/atom is unusually small. The packet does not establish whether the three seeds generated meaningfully independent fits.
- **Required correction** Verify initialization, data order, optimizer state, and checkpoint identity for each seed. Report both aggregate and pair-level dispersion without treating the present value as an error unless the provenance check proves one.

### R1-m3

- **Severity** Minor
- **Axis** Readability for nonspecialists
- **Affected element** Model hierarchy and metric interpretation
- **Evidence pointer** Stated project logic and benchmark design.
- **Issue** The distinction among QH-only, QH-residual, representation-only, and thermo-form is concise but insufficient for a nonspecialist to understand what information and flexibility each model receives.
- **Required correction** Add a compact model table specifying inputs, fitted terms, parameter counts, training scope, physical assumptions, and expected failure regime.

### R1-m4

- **Severity** Minor
- **Axis** Readability, technical soundness
- **Affected element** Reliability status
- **Evidence pointer** Phase inventory and QH reliability rule.
- **Issue** Readers could miss that most transition pairs contain at least one QH-unreliable phase.
- **Required correction** Display the reliability flag beside every phase-level and pair-level QH result and state whether each reported predictor used, held, or replaced the unreliable contribution.

### R1-m5

- **Severity** Minor
- **Axis** Readability
- **Affected element** Completion status
- **Evidence pointer** Project task list and missing-material list.
- **Issue** Completed results, diagnostic upper bounds, planned experiments, and unavailable evidence are interleaved.
- **Required correction** Separate frozen evidence from planned work and label every table or conclusion as completed, provisional, or not yet evaluated.

## Technical failings that need to be addressed before the case is established

The blocking technical failings are R1-M1 through R1-M6. The most consequential are the absence of a demonstrated QH advantage, missing learned LOPO and LOSO results, failure of the QH reliability assumption for five phases, an unresolved representation and reference error budget, incomplete leakage controls, and the lack of a reproducible frozen artifact.

## Assessment against Nature-style criteria

- **Originality** Not assessable. The benchmark construction may be distinctive, but no prior-work comparison is available.
- **Scientific importance** Potentially strong. Reliable OOD free-energy prediction and transition detection are important problems. The current positive method claim is weak, while the benchmark and negative findings may still be significant.
- **Interdisciplinary readership** Plausible but not yet demonstrated. The work could connect machine learning, statistical mechanics, and materials modeling. The present four-system evidence base and incomplete phase and system tests limit the reach of the conclusion.
- **Technical soundness** Not established. The central comparison is unfavorable or mixed, two claimed OOD axes lack learned results, QH is unreliable for many phases, error propagation is missing, leakage review is unfinished, and provenance is incomplete.
- **Readability for nonspecialists** Partly assessable. The packet states the project logic clearly, but a full manuscript, figures, legends, model comparison table, and complete methods narrative are absent.

## Recommendation posture

The project requires substantial technical completion and claim recalibration before a strong Nature-style case can be assessed. The benchmark, reliability diagnostics, and negative controls are promising. The current evidence supports neither the broad title claim nor a causal conclusion that QH baselines enable extrapolation across temperature, phases, and systems.
