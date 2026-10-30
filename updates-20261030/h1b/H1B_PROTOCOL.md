# H1b protocol: cross-tier OOD phase-transfer stress test on ds62 SiO₂

Frozen before any DPA embedding of a target structure is scored, and before any target free-energy label is used. Approved by the project owner on 2026-10-30.

## 0. What H1b is and is not

| Evidence tier | Status | What it can show |
|---|---|---|
| DFT-level development evidence (PBE-D3(BJ) SiO₂, task checkpoint) | Exists; all four SiO₂ folds previously inspected | The mechanism candidates |
| **H1b: ds62 assessed-tier external stress test** | This protocol | Whether the pattern appears under a different reference tier and out-of-distribution targets |
| Independent DFT-level confirmation of H1 | **Not done** | — |

**H1b cannot confirm or refute the DFT-level H1.** It must never be reported as an independent confirmation.

## 1. Budget

0 SU on Expanse and Bohrium. No DFT, TI or MD job. Local CPU only: DPA descriptor inference, BurnMan evaluation of published parameters, and ridge readouts. No running calculation is touched. No further material search.

## 2. Representation

- **Checkpoint.** Si/O task checkpoint `model.ckpt-5000.pt`, SHA-256 `9247692687a6cac04298499646fce7b2cd49c31bedfb97eead12d79edc2337dd`.
- **Load-compatibility copy.** For deepmd-kit 3.1.3, two config keys unknown to that runtime are removed, both at inert values: `descriptor.default_chg_spin = None` and `descriptor.repflow.sequential_update = False`. Weights are unchanged. Patched SHA-256 `9191c2f63350dde02bc229a00998bd6a728fb92478dcdfc4a3c2a97e89b4630c`; record in `stage0/CHECKPOINT_PATCH.json`.
- **Pooling.** `DeepPot.eval_descriptor` (fp64), unweighted mean over atoms, as in development.

**Gate E1.** Re-embedded source means must reproduce `task_checkpoint_embeddings.npz` to max abs difference ≤ 1×10⁻⁵. Otherwise **STOP**.

**Gate E2.** Target descriptors must be finite. Mean descriptors for each target's (3×3×3 or 4×4×6) supercell must agree with the unit cell to ≤ 1×10⁻⁶. Otherwise **STOP**.

## 3. Structures and phase mapping (audited before this freeze)

Records: `audit/STRUCTURE_AUDIT.json`, `audit/TOPOLOGY_AUDIT.json`.

| ds62 end-member | Role | Structure | Identity evidence | Documented deviation |
|---|---|---|---|---|
| q | source | development `quartz_beta` (243 atoms) | Si–Si coordination sequence 4,12,30,52,80 = qtz net, every Si | spglib C2 (not P6₂22); V/atom 12.00 vs ds62 V₀ 12.56 (−4.5%) |
| trd | source | development `tridymite_p63mmc` (216 atoms) | 4,12,25,44,67 = lon net, every Si | P1; V/atom 14.00 vs 15.50 (−9.7%) |
| crst | source | development `cristobalite_beta` (192 atoms) | 4,12,24,42,64 = dia net, every Si | P1 or P2₁/c; V/atom 14.29 vs 15.19 (−6.0%) |
| coe | target | COD 9007161, Angel et al. 2003, 0.1 MPa | C2/c; Si 4-coordinated; V/atom 11.392 vs V₀ 11.425 | — |
| stv | target | COD 9007150, Yamanaka et al. 2002, 0.1 MPa | P4₂/mnm; Si 6-coordinated; V/atom 7.769 vs V₀ 7.755 | — |

- **Excluded.** The C222₁ tridymite structure is excluded because ds62 has a single tridymite end-member.
- **Pressure sensitivity only (never used in any verdict).** COD 9007164 (coesite, 2.84 GPa) and COD 9007152 (stishovite, 9.26 GPa), from the same study series.
- **State mismatch (declared limitation).** Source representatives are near-zero-pressure, high-temperature, volume-scan configurations. Targets are ambient-pressure experimental structures. Labels are at 3 GPa and 10 GPa. The readout uses one structure per phase as phase identity, as in development.

## 4. Reference labels

- **Source.** ds62 through BurnMan 2.1.0 (`HP_2011_ds62.py` SHA-256 `0b14e147fa599a0adbae3adc6e956300d98bcf152bc633536fcd145eb10bede5`). G and S per atom (formula ÷ 3).
- **Temperature grid.** 500–2000 K, 1 K steps (1,501 points), fixed on validity grounds:
  - 500 K avoids the low-temperature extrapolation of the high-temperature trd/crst heat-capacity forms (the round-2 audit showed anomalous ΔS at 300 K).
  - 2000 K stays below SiO₂ melting at 3–10 GPa.
  - The window will not change after scoring.
- **Isobars.**
  - Primary: 3 GPa; target coe; primary pair coe − q.
  - Secondary: 10 GPa; target stv; primary pair stv − coe.
- **Reference entropy.** BurnMan's analytic S, not finite differences.

**Gate R1.** Re-evaluate the reference at 3 GPa and 10 GPa on the round-2 100 K grid and compare with `../fes-h1-candidate-audit-20261029/round2/SIO2_ASSESSED_TIERS_DS62_SLB.json`. Agreement must be ≤ 1×10⁻⁶ meV/atom. Otherwise **STOP**.

**T₀ rule.** T₀ is the temperature where the reference primary pair ΔG crosses zero on the 1 K grid, by linear interpolation. **Gate R2:** if there is not exactly one zero in the window, **STOP**.

**Reference uncertainty.** From the ds62 enthalpy covariance:
- σ(coe − q) = 0.1246 meV/atom.
- σ(stv − coe) = 1.3192 meV/atom.

Handling:
- Treated as a temperature-independent offset uncertainty on ΔG. S, V and Cp uncertainty is unavailable: unknown, not zero.
- States with |ΔG_ref| < 2σ count as reference-indeterminate. Stable accuracy is reported both unfiltered and with those states excluded.
- Reference Tc uncertainty is reported as 2σ / |ΔS_ref(T₀)|.

## 5. Readout (unchanged from the frozen A/B/C protocol unless stated)

- **Basis.** [1, x, x ln x, x²] with x = T/1000 K.
- **Alpha.** The same 13-value multiplier grid; alpha scale = trace(XᵀX)/dim; the same exact reduced-rank ridge solve.
- **Selection.** Source-validation pair ΔG MAE; ties go to the smaller multiplier. Target labels never enter fitting, scaling or selection.
- **Source label split on the grid index.** Index mod 10 = 5 is source validation; mod 10 = 0 is the source ID-test (reported only); all others are fit.
- **Arms.**
  - A: structure volume per atom. *Change from development:* computed from the representative structure rather than a reference table, so no target thermodynamic quantity enters.
  - B: task-checkpoint DPA.
  - C: DPA plus volume.
- **Scalings.** Coordinate (primary) and equal-block (secondary), as in development.

**Configurations**

| Test | Target | Primary two-source | Primary three-source | Also fitted and reported |
|---|---|---|---|---|
| Primary, 3 GPa | coe | {q, trd} | {q, trd, crst} | {q, crst}, {trd, crst} |
| Secondary, 10 GPa | stv | {q, coe} | {q, trd, coe} | all other 2-source subsets of {q, trd, crst, coe}; the 4-source set |

**Candidate pools for envelope scoring**

| Test | Primary pool | Secondary pool |
|---|---|---|
| Primary | {q, coe} | {q, trd, crst, coe} |
| Secondary | {coe, stv} | {q, trd, crst, coe, stv} |

Events are scored identity-aware: same from/to phases and direction, within 100 K.

## 6. Metrics

- Primary-pair gap MAE, bias and zero-gap skill (1 − MAE / MAE of ΔG = 0).
- Relative entropy MAE.
- Stable accuracy, unfiltered and σ-filtered.
- Hits, misses and false events; Tc error.
- **Anchor decomposition at T₀:** e_G(T) = e_G(T₀) − ∫_{T₀}^{T} δΔS dt, reported at T₀, 500 K and 2000 K.
- Out-of-distribution diagnostics from Stage 0, which are descriptive only and use no thresholds:
  - standardised distance from each target to the sources;
  - atomic nearest-neighbour distance quantiles relative to source leave-one-structure-out distances;
  - affine-hull weights and residual.

## 7. Pre-registered pattern criteria (B arm, coordinate scaling, primary pool)

| ID | Criterion |
|---|---|
| R1b | three-source primary-pair gap MAE < two-source |
| R2b | three-source stable accuracy < two-source, **and** three-source has 0 identity-correct target events |
| R3b | for the three-source model, \|e_G(T₀)\| > 100 K × \|ΔS_ref(T₀)\|, i.e. the anchor error alone moves the crossing beyond the 100 K tolerance to first order |

**Verdict rule** (applied separately to the primary and secondary tests; the primary test governs the headline):
- **PATTERN_PRESENT** if R1b, R2b and R3b all hold.
- **PATTERN_ABSENT** if R2b fails.
- **PARTIAL** otherwise.

Separately, and descriptively: whether any arm recovers the reference target event, and the zero-gap skill of each arm.

## 8. Stopping rules

- **Before scoring:** stop on failure of E1, E2, R1 or R2, or if the hash of this protocol does not match `H1B_PROTOCOL.sha256`. Report the blocker; do not score.
- **After scoring:** one run only. No change of readout, window, isobar, scaling, structures, criteria or checkpoint. Every result is retained, including failures and the sensitivity structures. No material-search extension follows from any outcome.
