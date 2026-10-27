# FES-Bench: reliability limits and controlled failure diagnosis

FES-Bench separates source-fit quality, transfer of phase-relative information and thermodynamic decisions. Two findings anchor the narrative: lower mean ΔG error can select the wrong phase; high envelope accuracy can coexist with incorrect pair-crossing topology. Source-fit distributions and the audited target-phase participation explain what these scores measure.

## Where -> What -> Why

### Where do methods fail?

Compare tasks under a stated input budget. Record measured phase errors, event misses, applicability limits and missing evaluations.

I · Failure landscape

### What relative information is lost?

Inspect the unanchored phase gap, relative entropy slope and the margin between the lowest competing phases.

II · Thermodynamic information loss

### Why, and can physical priors help?

Use matched controls to separate missing inputs, fit quality, reference mismatch and baseline validity.

III · Physical baseline contribution

## Four Results sections

1. Failure landscape: methods by system × task, with legal inputs and reported decision metrics.
2. Thermodynamic information loss: relative gap, entropy slope and the critical low-phase margin.
3. Physical baseline contribution: implementation and fit checks, then matched input replacements.
4. Applicability and transfer limits: compatible reference layers and state-specific physical assumptions.

## Landscape definitions

- OBS: a reported score; it is not an invented success threshold.
- NR: no matched reported evaluation in this view.
- NA: the named control needs target inputs forbidden by this task.
- RF: the compatible reference for the named new grid is missing.
- NE: no reference event occurs in the evaluation domain. False predicted events still matter.

Each row retains its own protocol and label budget. Seed extrema are descriptive ranges, not uncertainty intervals. Colour uses the worst reported seed. Source-specific SiO₂ residual and QH scores supplement the older unified table without changing it.

## One diagnostic framework

### Single-phase shape is accurate, but ΔG is wrong.

Candidate issue: Phase-dependent offsets or reference alignment.

Required test: Replace the static gap; evaluate unanchored ΔG.

Evidence limit: Diagnostic hypothesis; not a cause assignment.

### Train fit is good, but the temperature tail drifts.

Candidate issue: Relative entropy or higher-order temperature response.

Required test: Compare slope and a matched temperature-function control.

Evidence limit: Check train fit before calling this transfer failure.

### Mean ΔG MAE is low, but the stable phase is wrong.

Candidate issue: A small margin between the lowest competing phases.

Required test: Inspect the responsible lower-envelope pair and error direction.

Evidence limit: SiO₂ has a counterexample; its pair mechanism remains open.

### Performance drops when a phase is held out.

Candidate issue: Structure-to-thermodynamics mapping for an unseen phase.

Required test: Compare C-simple, C-DPA and a matched few-shot extension.

Evidence limit: Two current outer folds; no general cause established.

### QH residual does not improve a matched control.

Candidate issue: The baseline may add no usable information for that budget.

Required test: Compare the same temperature prior and verify the QHA definition.

Evidence limit: No consistent known-phase gain; causal transfer test open.

### Performance drops across systems.

Candidate issue: Energy references, representation or thermal-response differences.

Required test: Use family folds, static-gap replacement and matched support.

Evidence limit: Do not infer the cause from one Ti fold.

## Proposed matched replacement test

M₀: structure with common state inputs. M₁: the same pipeline plus a simple thermal prior. M₂: the same pipeline plus verified QHA. Keep labels, support, capacity and selection matched. Record the extra physical input and cost of QHA. This full causal comparison remains proposed.

## Evidence and supplementary history

The supplied first-principles outline sets the explanation order. Measurements retain their archived protocol and source. This guide does not promote historical scores into a matched leaderboard.

## Physical starting point

An MLIP predicts energy and force for configurations. Free energy requires their statistical weights under a defined Hamiltonian and thermal convention. QHA, sampling, integration and thermodynamic supervision provide different information. A derivative of a fitted G curve is its implied entropy; it does not establish where that information came from.

The current question is: can atomistic representations and physical priors predict phase-relative free-energy information for a new phase or system without its thermodynamic labels?

## Supplementary history of the question

### Start with residual learning

Question: Can a physical baseline simplify G?

Observation: Early SiO₂ tests motivated the idea. Their targets and protocols differ from the later matched tests.

Limit: Keep early metrics out of the current leaderboard.

Evidence: E019

### Check the formula before the physics

Question: Does inference implement the declared model?

Observation: Historical static-energy double counting required a correction. The selected Hf reconstruction later passed.

Limit: One passed checkpoint does not validate every historical run.

Evidence: E003, E004, E006

### Add the analytic control

Question: Is the temperature function itself difficult to fit?

Observation: Thermo4 fits the known-phase tails well. Adding the declared QH baseline gives no consistent pairwise gain.

Limit: The control uses target training labels. It is not zero-shot.

Evidence: E001, E005

### Separate error from phase decisions

Question: Does lower MAE imply the correct stable phase?

Observation: The two SiO₂ seeds provide a counterexample. Event scores also require an actual reference event.

Limit: Phase classification in a tail is not boundary recovery.

Evidence: E001, E007, E009

### Test unseen phases and systems

Question: Can structure predict coefficients without target labels?

Observation: The C readout now has one SiO₂ LOPO and one Ti LOSO result. C-simple has lower pair MAE in both.

Limit: Two folds do not establish broad transfer or its cause.

Evidence: E002, E008

### Make the reference target comparable

Question: Do the labels describe compatible physics?

Observation: Carbon has a scoped paired QHA grid. Ti has a single-phase QHA layer. Other source and production layers remain separate.

Limit: Eight indexed systems do not mean eight complete grids.

Evidence: E010, E011, E012, E013, E014, E015, E016, E017, E018

## Hypotheses and required controls

### H0 · QH + residual should give a large accuracy gain.

State: Not established

Question: Does the same function fit better after a physical baseline is added?

Observation: The matched analytic controls show no consistent pairwise gain. Early SiO₂ tests used historical protocols and cannot establish a common ranking.

Control: Compare direct Thermo4 with baseline + the same four-term residual. Keep labels, splits and selection fixed.

Next test: Test whether the baseline helps with less target information or with an unseen phase.

Evidence: E001, E003, E019

### H1 · QH supplies useful temperature dependence.

State: Clue, not a transfer result

Question: Can a thermal-shape gain preserve the phase gap?

Observation: Historical anchored single-phase curves motivated this idea. Phase-specific offsets can still change ΔG and reverse the phase order.

Control: Measure single-phase shape, pairwise gap and slope separately. State the actual QH definition.

Next test: Use a valid QHA layer and compare a matched smooth temperature prior.

Evidence: E003, E007, E012, E015

### H2 · Mean ΔG error does not determine phase accuracy.

State: Bounded counterexample

Question: Can a lower mean error still select the wrong stable phase?

Observation: Yes. In the two SiO₂ test seeds, the Phase-ID MLP has lower mean MAE than T3 but zero stable-envelope accuracy. T3 reaches 100%.

Control: Inspect the gap between the lowest two phases and the direction of its prediction error.

Next test: Locate the pair that reverses the stable envelope. The test has no reference transition.

Evidence: E001, E009

### H3 · Structure can predict transferable thermodynamic coefficients.

State: Two folds tested

Question: Can a model predict a new phase or system without its free-energy labels?

Observation: C-simple has lower target-pair MAE than C-DPA in the two completed folds. This does not prove a general descriptor ranking.

Control: Use source-only selection. Match target-label budgets, baselines and energy references.

Next test: Complete more LOPO/LOSO folds and seeds, plus matched T2/T3 transfer.

Evidence: E002, E008

### H4 · A physical baseline helps when it adds missing task information.

State: Open causal question

Question: Does QHA add transferable information beyond a smooth temperature function?

Observation: The present results do not isolate a causal vibration, electronic or magnetic contribution to model error.

Control: Match the Hamiltonian, nuclear statistics, electronic temperature, reference gauge and legal inputs.

Next test: Compare direct shape, QHA residual and simple structural controls on an admitted layer.

Evidence: E001, E007, E011, E012, E015

## Experiment register

### E001 · Matched temperature-tail tests

Ti / Zr / SiO₂ | Completed, bounded

Question: Can a learned model beat a simple function under the same label budget?

Observation: Direct Thermo4 has lower test ΔG MAE than the compared learned models. An analytic residual fit shows no consistent gain from the declared QH baseline.

Claim limit: These are fixed-pressure, known-phase tests. The test tails contain no reference transitions. Do not use them to rank zero-label transfer.

Direct sources:

- [S00001 · CROSS_SYSTEM_DIAGNOSTIC_REPORT.md](assets/sources/661483fc2ecfc9_CROSS_SYSTEM_DIAGNOSTIC_REPORT.md) · SHA-256 `661483fc2ecfc916eb3d1094b9627e63046a6d860c7dfdae56989ba6bf1e5310`
- [S00002 · MODEL_COMPARISON_BY_SEED.csv](assets/sources/eaa473677cf850_MODEL_COMPARISON_BY_SEED.csv) · SHA-256 `eaa473677cf850ad1b89f8af008187a8ff090a106cec4ceb3efa518f9b5a3aad`
- [S00003 · UNIFIED_TAIL_SCORES_20261020.csv](assets/sources/b8147d06672223_UNIFIED_TAIL_SCORES_20261020.csv) · SHA-256 `b8147d066722230507eea2e175d82d9eb8478c1751c5bd757c1f13a97400fd59`

### E002 · Structure-to-coefficient transfer

SiO₂ / Ti | Two folds completed

Question: Can structure predict thermodynamic coefficients without target free-energy labels?

Observation: C-simple has lower target-pair MAE than C-DPA in the SiO₂ cristobalite LOPO and Ti LOSO folds. Both match the SiO₂ envelope event. Neither matches the Ti event.

Claim limit: One held-out phase and one held-out system are not a broad material survey. Matched T2/T3 transfer remains open.

Direct sources:

- [S00004 · README_STATUS.md](assets/sources/24b3796f609442_README_STATUS.md) · SHA-256 `24b3796f6094423607606dfeecae4be3c96628bff4ce80b0936911cd60360ac1`
- [S00005 · pair_metrics.csv](assets/sources/8c633e77fdaf25_pair_metrics.csv) · SHA-256 `8c633e77fdaf25a20f7c8aecf2c85d2b50d86cd5714156d943fb21a367763159`
- [S00006 · envelope_events.csv](assets/sources/1ada7b3a46458b_envelope_events.csv) · SHA-256 `1ada7b3a46458b34cd3282592f3d8f85a6ddf0bb75db736121a1a2c2a1279bd6`
- [S00007 · fold_manifest.json](assets/sources/027d80739bfa20_fold_manifest.json) · SHA-256 `027d80739bfa20fe4e70ba4e225ba068cf09b2a76f194067d0328e4c07054df8`
- [S00008 · model_selection.json](assets/sources/e3c5d3973e6373_model_selection.json) · SHA-256 `e3c5d3973e63738363605fde838f09189c40fb4002c6e1d665f9bd58062e98bb`

### E003 · Rebuild the historical Hf prediction

Hf | Historical checkpoint

Question: Does inference apply the baseline, atom count and learned correction correctly?

Observation: The selected seed-11, step-1500 reconstruction agrees with archived inference. The first large error appears in the learned correction. The old initializer had only Si/O types.

Claim limit: This result applies to one old checkpoint. It does not describe later compatible initializers or prove that every historical run was correct.

Direct sources:

- [S00009 · DIAGNOSIS.md](assets/sources/d2eec2149f6e40_DIAGNOSIS.md) · SHA-256 `d2eec2149f6e400e6551cb361af9759df7e02b419dd1f5d0678fdc8ec373f29f`
- [S00010 · component_metrics.csv](assets/sources/56ec75d6f03b39_component_metrics.csv) · SHA-256 `56ec75d6f03b3989665467cc20d7e74762d8d1b55c8a0e59da1ba66f41b68c3c`
- [S00011 · component_pair_predictions.csv](assets/sources/a7fa26b7988978_component_pair_predictions.csv) · SHA-256 `a7fa26b7988978bb28c8b6d1ba92eb7185066dcd8ce80c59f91f35d3ea76c79c`
- [S00012 · qh_static_decomposition.csv](assets/sources/ee13d387ce436c_qh_static_decomposition.csv) · SHA-256 `ee13d387ce436c0d36aeccd7e66daf4a288ea050228e4e8b55cdf9ce5fb6091e`

### E004 · Hf-compatible initialization

Hf | Completed, one system

Question: Does a compatible element map remove the fit gap?

Observation: The compatible checkpoint loads and uses the correct Hf type index. Matched T2/T3 runs still trail analytic controls.

Claim limit: A valid type map does not prove useful thermodynamic pretraining. Keep each training recipe and seed separate.

Direct sources:

- [S00013 · DIAGNOSIS.md](assets/sources/8630884710c0a6_DIAGNOSIS.md) · SHA-256 `8630884710c0a60cb789a48259d14a3267c0758993ab38e7ced735ea73786312`
- [S00014 · CORRECTION_TO_PRIOR_HFINIT_REPORT.md](assets/sources/909cad92292e28_CORRECTION_TO_PRIOR_HFINIT_REPORT.md) · SHA-256 `909cad92292e2852c6b3531ba5dd445c254d081815cb6f1de77ba50e6f05e959`

### E005 · Batch size and training budget

Hf | Completed, bounded

Question: Does the training recipe explain part of the error?

Observation: Batch 8 improves T3 in the two tested seeds. T2 test error does not improve. A longer T2 run reduces training error without a matching test gain.

Claim limit: The longer run also changes learning-rate decay. This is not a pure step-count ablation. Seed spread measures training variability, not predictive uncertainty.

Direct sources:

- [S00015 · BATCH8_DIAGNOSTIC.md](assets/sources/a97742be653a6a_BATCH8_DIAGNOSTIC.md) · SHA-256 `a97742be653a6a91d03f512cc5c5d359e404315b9b711c85842cd6b1554b25af`
- [S00013 · DIAGNOSIS.md](assets/sources/8630884710c0a6_DIAGNOSIS.md) · SHA-256 `8630884710c0a60cb789a48259d14a3267c0758993ab38e7ced735ea73786312`
- [S00001 · CROSS_SYSTEM_DIAGNOSTIC_REPORT.md](assets/sources/661483fc2ecfc9_CROSS_SYSTEM_DIAGNOSTIC_REPORT.md) · SHA-256 `661483fc2ecfc916eb3d1094b9627e63046a6d860c7dfdae56989ba6bf1e5310`

### E006 · Task branches and alternate initialization

Hf | Partly complete

Question: Do task branches explain the pairwise prediction error?

Observation: Branches in the same frozen shared descriptor give nearly identical pair predictions.

Claim limit: This is not a contrast between independent descriptors. An HfO₂ branch is not proof of elemental-Hf pretraining.

Direct sources:

- [S00016 · DIAGNOSTIC_REPORT.md](assets/sources/a18692d118e5f3_DIAGNOSTIC_REPORT.md) · SHA-256 `a18692d118e5f31b63f99526c02c86620a000c31842a30e0a7a87e76c5766891`
- [S00001 · CROSS_SYSTEM_DIAGNOSTIC_REPORT.md](assets/sources/661483fc2ecfc9_CROSS_SYSTEM_DIAGNOSTIC_REPORT.md) · SHA-256 `661483fc2ecfc916eb3d1094b9627e63046a6d860c7dfdae56989ba6bf1e5310`

### E007 · Relative slope and curvature

Ti / Zr / Hf / SiO₂ | Descriptive evidence

Question: Does a model preserve the temperature response of the phase gap?

Observation: The archived Hf T2 curve has a large relative-slope error. Curvature error changes with the smoothing window.

Claim limit: These derivatives come from smooth reference tables. Window sensitivity prevents a stable curvature ranking or a causal attribution.

Direct sources:

- [S00017 · EIGHT_SYSTEM_EVIDENCE_REPORT_20261020.md](assets/sources/838b8d0b8c87b1_EIGHT_SYSTEM_EVIDENCE_REPORT_20261020.md) · SHA-256 `838b8d0b8c87b14e6d385c296946213c6385230eb76c3e59a85c29cdb87d8417`
- [S00018 · RELATIVE_SHAPE_WINDOW_SENSITIVITY_20261020.csv](assets/sources/fc484515a87b1e_RELATIVE_SHAPE_WINDOW_SENSITIVITY_20261020.csv) · SHA-256 `fc484515a87b1ef4f6fbb9eac4426031ee8cc4a6c3fe3a1eec4f1c9fe57cddc7`
- [S00019 · RELATIVE_SLOPE_CURVATURE_SUMMARY_20261020.csv](assets/sources/8eb92285f73f7a_RELATIVE_SLOPE_CURVATURE_SUMMARY_20261020.csv) · SHA-256 `8eb92285f73f7a0ecab4082e84ba8046888bb177bb41266b1456652f483c99b9`

### E008 · Frameworks and input budgets

Eight-system scope | Historical matrix

Question: Which methods can enter each task with legal inputs?

Observation: The A–G register separates method applicability, available results and information budgets.

Claim limit: The old matrix predates the two C transfer folds. Its entry that C had not run is now historical.

Direct sources:

- [S00020 · METHOD_FRAMEWORK_A_G_APPLICABILITY_FULL_20261020.csv](assets/sources/bedccd7dc57aee_METHOD_FRAMEWORK_A_G_APPLICABILITY_FULL_20261020.csv) · SHA-256 `bedccd7dc57aee8f44096786d868bce8a9c3af363142fbdb547466c0b7d63086`
- [S00021 · METHOD_FRAMEWORK_A_G_MATRIX_20261020.png](assets/sources/3cd81db21496ab_METHOD_FRAMEWORK_A_G_MATRIX_20261020.png) · SHA-256 `3cd81db21496ab157fe75c0ce266d76cf863d4e59a297f01369f5ed6b29166a5`

### E009 · Frozen splits and event definitions

Core systems | Implemented; comparisons open

Question: Can the test measure a stable transition without selecting it from the reference curve?

Observation: The temperature-tail split uses T/P grids. Pair crossings and stable-envelope transitions have separate scores.

Claim limit: With no reference test event, recall and matched Tc error are not applicable. False predicted events can still be counted.

Direct sources:

- [S00022 · temp_extrap.json](assets/sources/30ede38197139d_temp_extrap.json) · SHA-256 `30ede38197139d5e2e7aebb33580f2e050c71d1062b68dbe807abdeba9c9d845`
- [S00023 · phase_lopo.json](assets/sources/ddce0b06a8653a_phase_lopo.json) · SHA-256 `ddce0b06a8653aca50dec71278671ab5215e7fe65b984fdef9a0923216e7ea78`
- [S00024 · system_loso.json](assets/sources/136a377ef10c77_system_loso.json) · SHA-256 `136a377ef10c77c67b8bc0657037eb8ed3146ae82ce212b6a304b569c90becc0`
- [S00025 · metric_definitions.csv](assets/sources/23bf25b69cb911_metric_definitions.csv) · SHA-256 `23bf25b69cb911880da6c62809de2579c904e2b459364843337f19e60e51107b`
- [S00001 · CROSS_SYSTEM_DIAGNOSTIC_REPORT.md](assets/sources/661483fc2ecfc9_CROSS_SYSTEM_DIAGNOSTIC_REPORT.md) · SHA-256 `661483fc2ecfc916eb3d1094b9627e63046a6d860c7dfdae56989ba6bf1e5310`

### E010 · Static DFT and reference checks

Ti / Zr / Hf / SiO₂ | Completed, static only

Question: Do structures and Hamiltonians support a meaningful energy comparison?

Observation: Ten new W1 static calculations support targeted Hf k-point, SiO₂ Hamiltonian and structure checks.

Claim limit: Static total energies do not provide finite-temperature G. Compatible zero points and settings remain essential.

Direct sources:

- [S00026 · EXECUTION_UPDATE_20261020.md](assets/sources/84c343b72e4771_EXECUTION_UPDATE_20261020.md) · SHA-256 `84c343b72e47718dbc75537e250cc4a4c8afd9f61c61c556464d1113fcded535`
- [S00027 · STRUCTURE_SELECTION.md](assets/sources/e8200b6ba23703_STRUCTURE_SELECTION.md) · SHA-256 `e8200b6ba23703f9dc5896ba41839f68b233fe0efda5629f1c70d851cd1fe2b9`
- [S00028 · COMPUTE_LEDGER.csv](assets/sources/31647756ebdfe6_COMPUTE_LEDGER.csv) · SHA-256 `31647756ebdfe662fb0303978617d9e241f506ad29b7574463780ccdab79176d`

### E011 · Carbon paired QHA grid

C | Scoped layer admitted

Question: Can both phases be compared on the same temperature-pressure grid?

Observation: The admitted grid contains 32 phase-G and 16 paired-ΔG states. The volume, phonon and independent reconstruction records are available.

Claim limit: This is quantum-harmonic QHA at Te = 0. Numerical bounds omit full systematics. The sign at 10 GPa and 1200 K is unresolved within the recorded bound.

Direct sources:

- [S00029 · DFT_QHA_MULTIPRESSURE_LAYER_REPORT.md](assets/sources/a4c44e2e3504c4_DFT_QHA_MULTIPRESSURE_LAYER_REPORT.md) · SHA-256 `a4c44e2e3504c455bdb71e0ce4eeb281a991b11f979cff018a486b61b005d8e6`
- [S00030 · EXTENDED_MULTIPRESSURE_G_DELTA_CANDIDATE.csv](assets/sources/78fddb748e3c3f_EXTENDED_MULTIPRESSURE_G_DELTA_CANDIDATE.csv) · SHA-256 `78fddb748e3c3faefdad120c360deab0fdc8e3429de24ebcc5160848027178ae`
- [S00031 · MASTER_MANIFEST_VALIDATION_REPORT.json](assets/sources/11e34404a8da60_MASTER_MANIFEST_VALIDATION_REPORT.json) · SHA-256 `11e34404a8da60b1eb275b263e1cc232c7b1c6f18107857020d6bdf16f5e45d3`

### E012 · Ti hcp QHA comparison surface

Ti | Comparison layer admitted

Question: What single-phase QHA response can be reproduced?

Observation: The package contains 217 F(V,T) and 155 G(T,P) rows, with q-mesh and fit sensitivities.

Claim limit: The calculation fixes c/a and Te = 0. It does not provide hcp–bcc ΔG or an anharmonic reference. Several convergence errors remain unknown.

Direct sources:

- [S00032 · README.md](assets/sources/5b520dff572eb1_README.md) · SHA-256 `5b520dff572eb1a74a4592eb64843f3457d920e5c489d94d95166d83657ba065`
- [S00033 · G.csv](assets/sources/ae71812def1dd1_G.csv) · SHA-256 `ae71812def1dd1b6a88d0c31b293c314a07ca727ed6fd1bfe2a8e7e267eb83b0`
- [S00034 · NUMERICAL_SENSITIVITY.csv](assets/sources/8103aff75a5594_NUMERICAL_SENSITIVITY.csv) · SHA-256 `8103aff75a559448281683a9ca56814d8bcd70363b2c3312898553b69a45f211`

### E013 · Four-phase SiO₂ source tables

SiO₂ | Source layer admitted

Question: Does the fourth phase lack G, or does it lack a matched input contract?

Observation: The source-v2 package contains 1,649 G rows per phase and 9,894 derived pair rows. Its exact pressure is 0.0001 GPa.

Claim limit: Source reuse is not new DFT. Uncertainty and covariance are unknown. Structural matching remains a separate requirement.

Direct sources:

- [S00035 · THERMODYNAMIC_SOURCE_ADMISSION_20261025.md](assets/sources/f7810e915b3d91_THERMODYNAMIC_SOURCE_ADMISSION_20261025.md) · SHA-256 `f7810e915b3d917a7574f7e8765474f967a0cc62b2da06eaaabdfe20dfefda14`
- [S00036 · README.md](assets/sources/dfc98dc22dd6c3_README.md) · SHA-256 `dfc98dc22dd6c335b97d2f7e691c4988f15ded88c59626f83ff613c97afa3ae8`
- [S00037 · MANIFEST_SIO2_SOURCE_V2_INTEGRATION_20261025.json](assets/sources/27f8e8d344b65f_MANIFEST_SIO2_SOURCE_V2_INTEGRATION_20261025.json) · SHA-256 `27f8e8d344b65f7b369cb2c6d9c1b50ad52684ba4ce4ba027a832c6d96b97985`

### E014 · Ti/Zr/Hf source archives

Ti / Zr / Hf | Provenance checked

Question: Can the original TI and effective-harmonic assets be traced?

Observation: Six phase archives were recovered and checked. The audit corrected the use of a representative structure hash as a source-G structure hash.

Claim limit: G is an ensemble quantity. It need not map to one unique microscopic structure per row. Some reference conventions and errors remain unknown.

Direct sources:

- [S00038 · DARUS3582_SOURCE_STRUCTURE_REVIEW_20261024.md](assets/sources/efed65deaa38fe_DARUS3582_SOURCE_STRUCTURE_REVIEW_20261024.md) · SHA-256 `efed65deaa38fe2ef7f6b73f9cad09a7bab91e02029921e6318a05af4e7f472c`
- [S00039 · DARUS3582_SOURCE_CONTRACT_20261024.md](assets/sources/f0535e5b5a50b4_DARUS3582_SOURCE_CONTRACT_20261024.md) · SHA-256 `f0535e5b5a50b42f21b2941dcb3c490174e045adcbc1586782589ab0b1dca30a`

### E015 · Six-phase phonon checks

MgSiO₃ / C / CaSiO₃ | Completed; interpretation updated

Question: Which phases need a different free-energy route?

Observation: Fifty displacement outputs passed the force-output audit. Soft modes were found. Later fine-k graphite checks removed the old Γ shear instability.

Claim limit: One stable central volume does not validate a whole QHA grid. The old graphite rejection was superseded by later numerical checks.

Direct sources:

- [S00040 · SIX_PHASE_PHONON_ACCEPTANCE_20261022.md](assets/sources/ab83fbe5f98ddf_SIX_PHASE_PHONON_ACCEPTANCE_20261022.md) · SHA-256 `ab83fbe5f98ddf20e2d795f284c99786d81339b15f256530b0b7df7cfbaeebf1`
- [S00041 · GRAPHITE_SOFT_MODE_DECISION_20261024.md](assets/sources/ac56bf43781ec0_GRAPHITE_SOFT_MODE_DECISION_20261024.md) · SHA-256 `ac56bf43781ec0b64c715af89c5afda613723f0ebd843ac0ef706f8951b5d448`

### E016 · MgSiO₃ PPv soft-mode route

MgSiO₃ | Not closed

Question: Can a compatible two-phase reference be completed?

Observation: Small-amplitude fixed-cell curvature and the supercell contrast support a soft mode under the tested settings.

Claim limit: Do not remove imaginary modes and call the result ordinary QHA. The new finite-temperature or integration route has not closed the paired grid.

Direct sources:

- [S00042 · STATE_ALIGNED_NEXT_STEP_REVIEW_20261023.md](assets/sources/e4c89d9ca79e87_STATE_ALIGNED_NEXT_STEP_REVIEW_20261023.md) · SHA-256 `e4c89d9ca79e87dbd781032e61ab17a3e97bc5b0a908cf468bb670e9e2c9caac`
- [S00043 · MASTER_STATE.md](assets/sources/cf910281244890_MASTER_STATE.md) · SHA-256 `cf91028124489097d6affb8ecd8144d110a0c7aa62e0bad399c15ba099659844`

### E017 · Ca volume brackets and phase identity

CaSiO₃ | Not closed

Question: Does sampling represent the target state and yield a two-phase reference?

Observation: The V319 mean T/P is near the target, but distortion coordinates have few effective samples. Two tetragonal volumes have F(V,T) candidates.

Claim limit: Short MD is not a free energy. Average P1 symmetry alone does not reject the phase. Keep new PBE-D3(BJ) and source DP-LDA layers separate.

Direct sources:

- [S00044 · EXECUTION_REFRESH_20261026_2100UTC.md](assets/sources/c503ad4ba19d34_EXECUTION_REFRESH_20261026_2100UTC.md) · SHA-256 `c503ad4ba19d342ae0a2c3fe2a48a6c03604b4ab96f19687357622d0aa0dab42`
- [S00045 · EXECUTION_REFRESH_20261026_1921UTC.md](assets/sources/dc2f5546c4b825_EXECUTION_REFRESH_20261026_1921UTC.md) · SHA-256 `dc2f5546c4b825b76fcbefb1534f7e33048b9215a9dffe4c12afce3745a4411a`
- [S00046 · README.md](assets/sources/2d49115b9ee302_README.md) · SHA-256 `2d49115b9ee302f24ee78647616116d13419f279f1a1bb016a574e94ca3e87b1`

### E018 · Fe sampling and TI paths

Fe | Path revised; no new integral

Question: Can the sampled reference path connect to the target phases?

Observation: Twelve source fcc/hcp pairs have scoped admission. After disorder in the SC→DP midpoint, a direct Einstein→DP route was started.

Claim limit: No complete new TI/FEP label has been delivered. The anchor does not fill missing bcc data. QE and VASP reference layers remain distinct.

Direct sources:

- [S00047 · README.md](assets/sources/c8415b406e6a23_README.md) · SHA-256 `c8415b406e6a232e004f26aa466e7a76c8ef1e4cfe416631c659deae2ab92d23`
- [S00048 · EXECUTION_REFRESH_20261026_1946UTC.md](assets/sources/8a6772f8b75f6b_EXECUTION_REFRESH_20261026_1946UTC.md) · SHA-256 `8a6772f8b75f6b90c6b9a7c652be3809fc8ebe18a717f3d74aba789bb7af045d`
- [S00049 · FE_TI_CALCULATOR_VALIDATION_20261024.json](assets/sources/a3ffbca45818e3_FE_TI_CALCULATOR_VALIDATION_20261024.json) · SHA-256 `a3ffbca45818e3fcc8867bf6a602909dd41e9c8a17c170ecee1f9741c0ce8c3d`

### E019 · Historical figures and manuscript gaps

Paper-wide | Historical / open

Question: Which claims have valid evidence, and which remain placeholders?

Observation: Historical seed-23 figures, old QH definitions, corrected runs and pending result slots remain explicitly marked.

Claim limit: A figure placeholder is not a result. A single-system diagnostic must not become a general model claim.

Direct sources:

- [S00050 · historical_figure_register.csv](assets/sources/ea71e6f7261735_historical_figure_register.csv) · SHA-256 `ea71e6f7261735cccbf1b218009b8bf110e59d55758b7393d37b87e9aadbf5c7`
- [S00051 · fill_in_register.csv](assets/sources/3224aef3004e2b_fill_in_register.csv) · SHA-256 `3224aef3004e2be1fe5f1042c830f41357904c7a68be656731a352a8536c9529`
- [S00052 · 05_results.tex](assets/sources/6ea088af5fdf2e_05_results.tex) · SHA-256 `6ea088af5fdf2ec309c2fd7bd23bd022eb6cb1fb723ce1129f607434ad243e92`
- [S00053 · A1_supplementary.tex](assets/sources/079b79ec716cb3_A1_supplementary.tex) · SHA-256 `079b79ec716cb351f4ea8b1b09766ca170c377ef16ea7b1a7e144299d1d6e2f3`

