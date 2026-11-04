# Frozen Reviewer 2 Report

| Field | Review setting |
|---|---|
| Review profile | Nature-style project review with emphasis on originality and scientific importance |
| Mode | Review only |
| Source basis | Immutable review packet only |
| Package stage | Project or proposal with partial evidence |
| Readiness | The benchmark framework is promising, but the title-level scientific case is not established |

## Overall assessment

The project addresses a consequential problem. Free-energy errors of only 1 meV per atom correspond to roughly 132 to 134 K for the two reported SiO2 transitions, so extrapolative accuracy has direct physical importance. The package also shows unusually useful candour about negative results and incomplete work.

The current evidence does not support the central claim that quasi-harmonic baselines make learned free-energy surfaces extrapolate across temperature, phases and systems. On the completed temperature split, QH-residual has a delta-G MAE of 31.239 meV per atom, compared with 4.521 for the zero baseline. Its ranking accuracy is 0.1045, and it has 12 aggregate missed crossings across three seeds. Representation-only performs similarly in MAE. In contrast, the thermo-form fit reaches 0.401 meV per atom and 2.796 K transition-temperature error. This result is scientifically interesting, but it redirects the contribution toward constrained thermodynamic forms rather than establishing the title claim.

Originality cannot yet be judged because the packet contains no complete prior-work comparison or evidence supporting the novelty wording. Importance of the scientific problem is supported. Importance of the proposed learned solution remains weak because the central method currently underperforms simple baselines and has not been evaluated on the phase and system OOD splits.

## Who would be interested and why

Computational materials scientists would care about reliable phase-stability prediction under temperature extrapolation. Researchers developing machine-learned interatomic potentials and thermodynamic surrogates would care about whether physical baselines improve out-of-distribution learning. Researchers studying Hf, Ti, Zr and SiO2 phase transitions would care about the reported crossing sensitivities and failure modes. Benchmark and reproducibility researchers would care about the frozen split definitions, leakage controls, provenance design and distinction between interpolation and several forms of extrapolation.

Broader interdisciplinary interest is plausible because the project tests a general proposition about combining physical structure with machine learning. That interest will depend on establishing the proposition across the declared OOD axes and explaining the physical logic for nonspecialists.

## Major strengths

- The project exposes negative evidence rather than hiding it. The narrative records that QH does not clearly simplify the Hf, Ti and Zr residuals in the no-learning diagnostic.

- The benchmark separates temperature, phase and system extrapolation and records hashes for the split definitions.

- Evaluation extends beyond aggregate MAE to sign, ranking, crossing, transition-temperature and pair-coverage measures.

- The crossing slopes make the physical consequences of apparently small energy errors legible.

- The single-side QH experiment and no-learning negative control begin to distinguish learned correction from genuine physical simplification.

- The thermo-form result identifies a potentially important alternative scientific direction.

## Major Concerns

### R2-M1 Major, Blocking Yes

- **Claim pointer** The title-level claim that quasi-harmonic baselines make learned free-energy surfaces extrapolate across temperature, phases and systems.

- **Evidence pointer** Temperature extrapolation results and single-side QH diagnostics.

- **Evidence status** Confirmed conflict between the broad claim and the reported completed results.

- **Concern** QH-residual performs substantially worse than the zero baseline on delta-G MAE. It also has poor ranking accuracy and multiple missed crossings. Its MAE differs little from representation-only, while the thermo-form fit is much better than both.

- **Why it matters** The central intervention is not shown to deliver the claimed benefit. The strongest result currently supports a different model class.

- **Resolution test** Demonstrate a prespecified and practically meaningful QH-residual advantage over representation-only and relevant simple baselines on the declared OOD tasks, including crossing outcomes. If that test is not met, narrow the title and contribution to the observed conditional findings and foreground the thermo-form result.

### R2-M2 Major, Blocking Yes

- **Claim pointer** Novelty language in the title and project framing.

- **Evidence pointer** Missing material affecting confidence, which states that no complete prior-work comparison or support for the novelty wording was supplied.

- **Evidence status** Not assessable.

- **Concern** The packet does not establish what is original relative to existing free-energy fitting, physics-informed residual learning, QH correction or extrapolation benchmarks.

- **Why it matters** A Nature-style originality judgment requires a defined advance over the strongest relevant prior approaches, not only a new internal benchmark.

- **Resolution test** Provide a claim-matched prior-work comparison that identifies the nearest methodological and benchmarking precedents, states the precise difference, and tests against the strongest relevant available baselines.

### R2-M3 Major, Blocking Yes

- **Claim pointer** Extrapolation across temperature, phases and systems.

- **Evidence pointer** Learned QH-residual and representation-only results are absent for LOPO and LOSO. The three-way controlled ablation is incomplete.

- **Evidence status** Confirmed missing evidence.

- **Concern** Only temperature extrapolation has learned-model results. Phase and system extrapolation currently contain only nonlearned baselines.

- **Why it matters** Evidence from one OOD axis cannot establish performance across all three axes. It also cannot isolate whether any benefit comes from QH, representation, model form or split-specific temperature structure.

- **Resolution test** Complete the same controlled QH-residual, representation-only and relevant constrained-model comparisons on every frozen OOD axis, with identical reporting of errors, rankings, crossings, coverage and seed uncertainty.

### R2-M4 Major, Blocking Yes

- **Claim pointer** Generalization across phases using QH baselines.

- **Evidence pointer** Only four phases satisfy the 1 percent imaginary-mode reliability rule. Metal bcc phases and two SiO2 phases exceed the threshold. Mixed one-side QH gives 644.868 plus or minus 465.591 meV per atom.

- **Evidence status** Confirmed.

- **Concern** QH is reliable on only part of the stated phase domain, and the mixed reliable versus held setting performs extremely poorly. The physical applicability boundary is therefore central rather than peripheral.

- **Why it matters** A baseline cannot support a broad across-phase claim if it is unavailable or unreliable for many target phases unless the method explicitly handles that regime and demonstrates success.

- **Resolution test** Justify and sensitivity-test the 1 percent rule, define the prediction protocol for unreliable phases before evaluation, and show performance separately for reliable, unreliable and mixed-regime pairs. Scope the claim to the supported regime if robust coverage is not achieved.

### R2-M5 Major, Blocking Yes

- **Claim pointer** Out-of-distribution extrapolation and model comparison.

- **Evidence pointer** The leakage audit is ongoing. Hf LOSO is confounded by temperature range, and the overlap-temperature subset is only partially reported.

- **Evidence status** Confirmed unresolved validity risk.

- **Concern** The package has not yet established that preprocessing, constants, hyperparameter choices and derived inputs are fit strictly within training boundaries. One system split also changes system identity and temperature range together.

- **Why it matters** Leakage or distributional confounding could make apparent extrapolation success uninterpretable.

- **Resolution test** Complete the leakage audit with explicit provenance for every fitted quantity. Report the Hf overlap-temperature evaluation for all learned and nonlearned models and distinguish system transfer from temperature-range transfer.

### R2-M6 Major, Blocking Yes

- **Claim pointer** Transition prediction and seed robustness.

- **Evidence pointer** QH-residual reports 102.379 K mean numeric transition-temperature error together with 12 aggregate missed crossings. Seed standard deviation is about 5.4e-9 eV per atom, while checkpoint files are absent and metric provenance records an unknown commit.

- **Evidence status** Confirmed reporting gap with unverified run independence.

- **Concern** The summary does not establish how missed crossings affect the reported transition-temperature error, and the almost identical seed aggregates cannot be independently verified.

- **Why it matters** A conditional error over successful crossings can appear favorable while omitting failures. Unverified seed independence weakens uncertainty claims.

- **Resolution test** Report per-pair and per-seed crossing status, the denominator for every transition metric, an explicit policy for missed and false crossings, and confidence summaries that include coverage. Restore checkpoint-level provenance and verify that each seed represents an independent run.

### R2-M7 Major, Blocking Yes

- **Claim pointer** Reproducible quantitative evidence.

- **Evidence pointer** The tracked notebook has no saved outputs and differs materially from the supplied export. The GitHub source is not commit-pinned. The local tree is dirty. Sections A to C are appended after export. Section C has no standalone producer. QH-residual metrics record an unknown commit, and referenced checkpoints are absent.

- **Evidence status** Confirmed.

- **Concern** No single immutable, executable release currently links source data, code, checkpoints, tables, figures and manuscript claims.

- **Why it matters** The main numerical evidence cannot be independently reconstructed or tied unambiguously to a versioned analysis.

- **Resolution test** Produce a clean commit-pinned release that executes end to end, saves outputs, exports every cited diagnostic, hashes every input used by figures and tables, includes required checkpoints or lawful immutable references, and reproduces the reported values.

### R2-M8 Major, Blocking Yes

- **Claim pointer** Extrapolation across systems and broad scientific importance.

- **Evidence pointer** The active domain contains four systems, while family LOSO is restricted to Hf, Ti and Zr. CaSiO3 and SiO2 C2221 are excluded. No uncertainty analysis treats systems as inferential units.

- **Evidence status** Confirmed scope limitation.

- **Concern** The evidence base is too narrow and structurally uneven to support a general across-system conclusion. The available family transfer test covers three related metal systems, while SiO2 belongs to a different reference-data route.

- **Why it matters** Broad scientific importance depends on distinguishing a general transferable principle from behavior specific to a small selected set.

- **Resolution test** Define the target population of systems, justify the sampling and exclusions, report system-level variability, and frame inference around the observed systems unless additional independent systems support broader generalization.

### R2-M9 Major, Blocking Yes

- **Claim pointer** Accuracy of the learned free-energy workflow.

- **Evidence pointer** Frozen DPA pairwise static-energy errors are 16.3 to 38.3 meV per atom. The fine-tuned PES acceptance file and TI upper bound are absent. Functional and truth sensitivity is unstarted.

- **Evidence status** Confirmed missing validation.

- **Concern** The package does not yet bound the contribution of upstream potential error or reference choice to the free-energy and transition errors.

- **Why it matters** The reported static-energy discrepancies are comparable to or larger than several target free-energy scales. Without sensitivity and acceptance evidence, failure cannot be assigned cleanly to the residual-learning hypothesis.

- **Resolution test** Supply the prespecified PES acceptance evidence, quantify propagation of upstream errors into delta-G and crossings, and complete functional or reference-truth sensitivity before attributing performance to the proposed architecture.

## Minor Comments

### R2-M10 Minor, Blocking No

- **Claim pointer** The project moves from a QH baseline to a sum of multiple cheap physical channels.

- **Evidence pointer** Stated model evolution in the project logic.

- **Evidence status** Confirmed.

- **Concern** The scientific object has shifted from a QH-specific hypothesis to a broader multi-physics residual framework, but the title remains QH-specific.

- **Why it matters** Readers need to know whether the paper tests one physical baseline or proposes a general architecture.

- **Resolution test** State one primary hypothesis and align the title, model definition, ablation plan and conclusions with it.

### R2-M11 Minor, Blocking No

- **Claim pointer** Accessibility to nonspecialists.

- **Evidence pointer** No full introduction, complete methods narrative, saved figures, legends or submission-ready discussion was supplied.

- **Evidence status** Not assessable.

- **Concern** Readability and conceptual accessibility cannot be judged from the partial package.

- **Why it matters** Interdisciplinary readers will need clear explanations of free-energy crossings, QH reliability, residual learning and why each OOD axis represents a distinct scientific challenge.

- **Resolution test** Provide a complete narrative and figures, then verify that each central physical quantity, comparator and failure mode is understandable without notebook context.

### R2-M12 Minor, Blocking No

- **Claim pointer** Interpretation of aggregate benchmark scores.

- **Evidence pointer** Reference delta-G variability spans very different scales across systems, and skill is defined relative to the zero baseline.

- **Evidence status** Confirmed.

- **Concern** Aggregate metrics may be dominated by systems or pairs with larger energy variation.

- **Why it matters** A single pooled MAE or skill value may obscure whether improvement is consistent across physically distinct transitions.

- **Resolution test** Report per-system and per-pair metrics alongside macro-averaged summaries and state the weighting used in every aggregate.

### R2-M13 Minor, Blocking No

- **Claim pointer** External benchmark strength.

- **Evidence pointer** External baselines are ongoing, and the complete panel is missing.

- **Evidence status** Confirmed incomplete comparison.

- **Concern** Bartel and simple constants do not by themselves establish competitiveness against the strongest relevant alternatives.

- **Why it matters** Scientific importance depends partly on whether the proposed insight changes practice relative to credible existing options.

- **Resolution test** Complete the prespecified external-baseline panel and explain which scientific capability each comparator tests.

## Technical failings that need to be addressed before the case is established

- Resolve the conflict between the title claim and the completed temperature results under R2-M1.

- Establish originality through a claim-matched prior-work comparison under R2-M2.

- Complete learned evaluations and controlled ablations across all OOD axes under R2-M3.

- Define and validate the QH applicability regime under R2-M4.

- Close the leakage and temperature-confounding audits under R2-M5.

- Make crossing metrics and seed uncertainty auditable under R2-M6.

- Release one clean, immutable and executable evidence package under R2-M7.

- Bound system-level generalization under R2-M8.

- Quantify upstream PES and reference sensitivity under R2-M9.

## Assessment against Nature-style criteria

| Criterion | Assessment | Basis |
|---|---|---|
| Originality | Not assessable | No complete prior-work comparison or evidence for the novelty wording is available |
| Scientific importance | Potentially high problem importance, weak contribution evidence | Transition temperatures are highly sensitive to small energy errors, but the proposed learned method currently underperforms simple baselines |
| Interdisciplinary readership | Promising but not established | The physics-informed learning question could interest materials, machine-learning and benchmark communities, but the evidence and narrative remain incomplete |
| Technical soundness | Weak at present | Frozen splits and broad metrics are strengths, while incomplete OOD tests, unresolved leakage, conditional QH applicability, missing uncertainty structure and broken provenance remain material |
| Readability for nonspecialists | Not assessable | A complete manuscript, saved figures, legends and methods narrative were not supplied |

## Recommendation posture

No editorial decision is offered. On the supplied evidence, the project should be regarded as an exploratory benchmark and hypothesis-testing programme rather than an established demonstration of QH-enabled extrapolation. The problem is important, the negative results are informative, and the thermo-form result may support a strong alternative contribution. A broad Nature-style case would require either decisive completion of the title-level tests or a narrower reframing around the conditional failure of QH residual learning and the apparent advantage of constrained thermodynamic forms.

## Functional-completeness retrospective

| Dimension | Status |
|---|---|
| Scope covered | Project narrative, benchmark design, current results, missing evidence and provenance as represented in the immutable packet |
| Authority and evidence | Limited strictly to the packet |
| Alignment | The title claim, methods, available results and missing tests were traced |
| Modules applied | General scientific integrity and computational or ML evaluation |
| Change propagation | Not applicable because no artifact was edited |
| Deterministic and visual checks | Not performed because source files and saved figures were outside the permitted evidence basis |
| Open issues | Nine Major blocking concerns and four Minor nonblocking concerns |
| Readiness | Not ready for a title-level Nature-style scientific case |
