# QH generator audit — 2026-10-27

## Verdict

The current legacy QH comparison cannot establish the physical efficacy of verified QHA. This audit confirms a missing imaginary-mode rejection gate and a contradictory quartz reliability flag. It also traces known unit and minimization defects to their specific historical branches. It does not establish that all QH scripts are wrong or that repairing them will improve learned transfer.

No original file, frozen prediction, score or production calculation was changed. No DFT, training, or full phonon recomputation was run.

## 1. Confirmed: the imaginary-frequency threshold does not reject a volume

`data/hf_sign_audit/qh_run.py` is the preserved source of the historical generator (SHA256 a1405fc03c19401479b067621d9c0a0f1c8c0f0e4b9373ae98f331e9a2987ce6).

- Lines 75–83 calculate thermal properties before computing negative-mode statistics.
- Lines 87–95 count negative frequencies but return a finite thermal table without rejecting unstable modes.
- Lines 169–172 add all volumes to the candidate list and choose the smallest E+F without using the negative-mode fields.

`check_archived_control_flow.py` executes the exact archived function/expression with synthetic fixtures. A spectrum containing −1.96 THz still returns a finite thermal value, and the exact minimum-selection expression chooses an unstable candidate if its reported energy is lower. This checks control flow, not Phonopy numerical thermodynamics.

Phonopy's documented default thermal cutoff excludes nonpositive frequencies from oscillator sums. Therefore obtaining a numerical thermal table is not proof that imaginary oscillators were physically treated. The exact archived library version and its integrated-mode count still need recovery; current documentation alone does not prove every detail of the historical execution.

Documentation: https://phonopy.github.io/phonopy/setting-tags.html#cutoff-frequency

## 2. Confirmed: quartz reliability metadata conflicts with its own raw output

Current source config `configs/qh/sio2_sse_pbe_materialize.yaml` explicitly supplies quartz_beta=true. `qh/materialize.py` lines 98–104 validates the boolean type; line 159 copies that boolean into the report. It does not derive the acceptance flag from the spectra. The same code copies Hf/hcp-specific explanation text to SiO2 reports.

All 1,649 raw quartz rows choose V/V0=1 and report minimum frequency −1.05403160233 THz, while `qh_reliable=true`. Raw F_QH values exactly match the canonical quartz fqh.csv on their common temperatures. This flag is not a verified harmonic-stability certificate. Whether the instability is structural, a potential error, force construction, acoustic error or insufficient numerical convergence needs mode-resolved investigation.

The other active SiO2 raw rows also select unstable settings: cristobalite −0.193488962524 THz at V/V0=1.02; P63/mmc tridymite −0.354003196531 THz at V/V0=1.02. Their flags already say false. Do not silently replace these selected points with a different volume merely to obtain a favorable result.

`BASELINE_LINEAGE_AUDIT.csv/json` records the source reports, selected settings, hashes and exact raw-to-canonical scalar check. The previous statement that quartz passed a phonon validity gate must be withdrawn for this baseline.

## 3. Confirmed historical unit bug: C2221 fixed-cell generator

`generate_sio2_qh_batched.py` line 53 converts kJ/mol primitive cell to eV/cell, but line 58 labels it eV/atom without dividing by primitive-cell atom count. Its identity primitive choice preserves the 24-atom input in the recorded C2221 branch. This caused the previously identified factor-24 error.

A corrected script and correction CSV already exist. This branch was excluded from the three active-phase scores. Therefore this bug does NOT automatically explain the Ti/Zr/Hf/three-phase SiO2 comparison table. It also produces fixed-cell F_vib only, not static+vibrational volume-minimized G.

Phonopy's thermal energy unit is per mole of the defined primitive cell, not automatically per mole of atoms. The main historical qh_run.py does contain an atom-count division; its primitive/unitcell equivalence still needs explicit assertion for arbitrary future cells.

Documentation: https://phonopy.github.io/phonopy/setting-tags.html#tprop

## 4. Confirmed historical minimizer defect; actual-data impact must stay separate

The pre-fix Ti pilot minimized among interior local minima and considered endpoints only if no local minimum existed. This can miss a lower boundary value. The exact archived block returns x=0 for f=x²−x³ on [−1,2], whereas the bounded global minimum is x=2 with f=−4.

A bounded-minimum fix comparing all stationary points and endpoints already exists. The earlier 155-row table reconstruction found no value change under that correction; this counterexample does not establish a new change to those actual rows.

A separate already-completed Ti raw-spectrum audit found that excluding residual Γ acoustic modes changes the early pilot G by up to 2.416 meV/atom at the reported 1200 K/0 kbar example. Later scoped Ti reference exports use the corrected handling. Do not confuse this independent DFT-QHA layer with the historical MLIP-QH input in the method-comparison table.

## 5. Deficient QHA approximation, rather than a demonstrated subtraction bug

The preserved MLIP generator uses a finite set of scaled cells, no internal-coordinate relaxation step, and selects the discrete minimum of E+F. It has no explicit PV term, no EOS fit, and no boundary-expansion acceptance gate. Hf selects the upper 1.02 endpoint through the recorded support. It is a discrete diagnostic, not a converged pressure-minimized QHA.

At the reference pressure 0.0001 GPa, Hf's missing PV is only about 0.014 meV/atom in the prior audit; it cannot explain hundreds of meV/atom alone. The sign and static+vibrational reconstruction checks previously passed. None of those checks certifies the phonon integration or physical validity.

Documentation: https://phonopy.github.io/phonopy/qha.html

## 6. What remains trustworthy, and what must be done next

- Existing numerical model scores remain records for their exact declared inputs. Do not delete unfavorable results or relabel them as verified QHA performance.
- The current four-system QH-vs-direct Thermo4 table measures a historical diagnostic-input effect. It cannot establish that valid vibrational physics fails to help.
- Carbon's independently reconstructed DFT harmonic-QHA layer and the later Ti scoped layer use separate scripts, spectra and acceptance. These legacy defects do not automatically invalidate them. Carbon's spectral/entropy reconstruction is a useful independent check, within its harmonic and convergence limits.

Minimal repair/validation order:
1. Freeze hashes and map every scored fqh.csv to the exact generator/config/cache, including archived Phonopy version and primitive atom count.
2. Recompute F_vib from existing full spectra with independent q-weighted oscillator sums; compare with historical thermal output and integrated mode counts. No new DFT is needed if the spectra are retained.
3. Reject genuine non-acoustic unstable modes; handle only verified translational zero modes specially. Do not clip, take absolute values or silently remove soft modes to manufacture QHA.
4. On stable reference states, check per-primitive/per-atom units, ZPE, static term counted once, E(V) geometry relaxation, EOS/PV, and interior minima/convergence.
5. If the same reference branch is physically unstable after numerical validation, declare ordinary QHA inapplicable; do not attempt to repair physics by changing a sign.
6. Only then rerun same-function controls and decide whether learned residual retraining is needed. Preserve original and corrected runs separately.

This first-pass audit is sufficient to suspend a physical interpretation of the old QH comparison. It is not a complete repaired baseline or a successful new method experiment.
