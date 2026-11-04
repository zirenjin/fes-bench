# Four correction routes: result (development validation)

- **Protocol.** `PROTOCOL.md`, SHA-256 0b382dba…, frozen 2026-11-02T17:44Z.
- **Deviations.** `DEVIATIONS.md`. D1: phonopy integrates every mode above 0 THz, so the bins now cover 1e-4–60 THz. D2: the gate is checked against the same-run phonopy F, and route C is applied as a difference on the stored F. Both were recorded before any full-data route result. The 20-formula smoke run was used only to test the code and is not reported.
- **Gate.** Re-integration matches the same-run phonopy F for 100% of entries (maximum deviation 0.006 meV/atom for F0 and 0.013 for MS). Recompute against stored values: median 0.007 meV/atom, p99 0.06–0.08.
- **Outputs.** `results/ROUTES.json`, with the outer out-of-fold G of every run in `results/OOF_G.npz` on GenSi; `scripts/routes.py`; `scripts/summarize.py`.
- **Data.** 521 entries, 239 formulas, 351 pairs; 101 pairs have a reference crossing. MLIP weights are frozen and there is no new DFT. All data were seen before, so every result is development evidence.

Numbers are from repeat 1 unless marked. "Δ" is route minus "none" for the same model and setting, with a 95% formula-cluster bootstrap interval. M1 is the formula-weighted RMS of e_Δ(T_c^ref) in meV/atom. "Wrong" is the wrong-ordering rate averaged over 600/1000/1500 K.

## S0: no reference information

| Model | Route | M1 | Δ M1 | Within ε₅₀ | Wrong | Δ wrong (points) | Pair RMS 1000 K | δΔS 1000 K | p95 \|e₁₀₀₀\| | Pairs worse by > 10 | False crossings | Verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| F0 | none | 7.09 | — | 0.12 | 11.4% | — | 9.71 | 8.2 | 17.2 | — | 87 | — |
| F0 | A | 5.91 | **−1.18 [−2.04, −0.27]** | 0.06 | 12.3% | +1.0 [−1.7, 3.6] | 8.74 | 8.2 | 16.8 | 6 | 109 | partial |
| F0 | B | 6.65 | −0.44 [−1.23, 0.58] | 0.10 | 15.0% | **+3.7 [0.9, 6.6]** | 8.85 | 7.0 | 18.8 | 22 | 85 | failure |
| F0 | C | 7.14 | +0.06 [−0.86, 0.84] | 0.09 | 12.9% | +1.5 [−0.7, 3.7] | 10.95 | 10.6 | 19.7 | 28 | 79 | failure |
| MS | none | 9.45 | — | 0.07 | 10.6% | — | 10.70 | 10.5 | 18.6 | — | 83 | — |
| MS | A | 7.26 | −2.19 [−4.08, 0.03] | 0.13 | 9.7% | −0.9 [−2.9, 1.1] | 8.90 | 10.5 | 17.8 | 17 | 105 | partial |
| MS | B | 7.39 | −2.06 [−4.67, 0.62] | 0.07 | 11.2% | +0.6 [−2.0, 3.1] | 9.74 | 8.1 | 23.3 | 22 | 84 | failure |
| MS | C | 7.66 | −1.79 [−3.34, 0.03] | 0.08 | 12.3% | +1.6 [−0.7, 4.1] | 10.06 | 11.5 | 16.5 | 26 | 82 | failure |
| D1 (w = 0.5) | against F0 | 6.92 | −0.17 [−2.16, 1.89] | 0.11 | **8.7%** | **−2.7 [−5.1, −0.5]** | 7.49 | 6.1 | **12.8** | 12 | 76 | partial |
| D1 (w = 0.5) | against MS | 6.92 | **−2.53 [−4.50, −0.35]** | 0.11 | 8.7% | −2.0 [−5.4, 1.3] | 7.49 | 6.1 | 12.8 | 7 | 76 | success |
| D2 (global w, chosen 0.5–0.6) | against F0 / MS | 7.01 | −0.08 / −2.44 | 0.12 | 8.9% | −2.5 / −1.7 | 7.77 | 6.3 | 13.6 | 11 / 9 | 84 | partial / partial |
| D3 (gated w) | against F0 / MS | 7.18 | +0.10 / −2.26 | 0.04 | 8.9% | −2.5 / −1.7 | 9.14 | 7.8 | 13.2 | 11 / 10 | 90 | partial / partial |

## S1: static ΔE₀ provided to every route

| Model | Route | M1 | Δ M1 | Within ε₅₀ | Wrong | Δ wrong (points) | Pair RMS 1000 K | δΔS 1000 K | p95 \|e₁₀₀₀\| | Pairs worse by > 10 | Verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|
| F0 | none | 4.95 | — | 0.16 | 7.9% | — | 7.69 | 8.2 | 14.9 | — | — |
| F0 | A / B / C | 4.94 / 4.62 / 4.89 | −0.01 / −0.34 / −0.06 (all intervals include 0) | 0.03–0.07 | 9.0 / 8.9 / 9.8% | +1.1 / +0.9 / +1.9 | 6.92 / 7.02 / 8.10 | 8.2 / 6.9 / 10.6 | 15.8 / 15.5 / 15.1 | 9 / 22 / 23 | all failure |
| MS | none | 7.37 | — | 0.18 | 5.3% | — | 10.03 | 10.5 | 13.9 | — | — |
| MS | A / B / C | 6.14 / 6.14 / 5.82 | −1.23 / −1.23 / −1.55 (all intervals include 0) | 0.12–0.13 | 5.8 / 7.9 / 8.6% | +0.5 / **+2.7** / **+3.3** | 7.15 / 8.05 / 8.95 | 10.5 / 8.0 / 10.8 | 11.9 / 18.2 / 16.2 | 6 / 19 / 30 | partial / failure / failure |
| D1 | against F0 | 4.95 | −0.00 [−2.20, 2.12] | 0.19 | **3.9%** | **−4.1 [−6.4, −2.1]** | 6.04 | 6.1 | **9.1** | 10 | partial |
| D1 | against MS | 4.95 | **−2.42 [−3.96, −0.45]** | 0.19 | 3.9% | −1.4 [−3.4, 0.4] | 6.04 | 6.1 | 9.1 | 3 | success |
| D3 | against MS | 5.53 | **−1.83 [−3.11, −0.22]** | 0.21 | 4.8% | −0.5 | 7.42 | 7.4 | 10.4 | 5 | success |

## Fit against transfer (repeat 1; R² of the stage targets)

| Run | Stage-2 location α: train / out-of-fold | Stage-1 response: train / out-of-fold | Pair-curve loss: train / test |
|---|---|---|---|
| F0 S0 A | 0.76 / 0.21 | — | 30 / 71 |
| MS S0 A | 0.54 / 0.34 | — | 57 / 81 |
| MS S1 A | 0.81 / 0.55 | — | 33 / 52 |
| F0 S0 B | 0.55 / 0.12 | δ: 0.82 / 0.19 | 41 / 71 |
| MS S1 B | 0.57 / **−0.19** | δ: 0.78 / 0.43 | 23 / 55 |
| F0 S0 C | 0.61 / 0.18 | θ: 0.86 / 0.29 | 61 / 105 |
| MS S0 C | 0.62 / 0.18 | θ: 0.30 / 0.11 | 73 / 93 |

**Stability over the three repeats** (range of Δ):
- Route A, M1: F0 S0 −1.18 to −0.21; MS S0 −3.23 to −2.19.
- Route C, wrong ordering: always worse, +1.2 to +4.0 points.
- D1 has no fitted parameter, so its result is identical in every repeat.

## Reading (development evidence)

1. **Location (A) is partly learnable, but not to the needed precision, and it does not improve ordering.**
   - Without reference information, A lowers the near-transition error by 1.2–2.2 meV/atom.
   - That is still about 15 times the ±50 K budget (median ε₅₀ about 0.4 meV/atom). The fraction of pairs within budget does not rise.
   - False crossings increase (87 → 109 for F0; 83 → 105 for MS), and wrong ordering does not improve.
   - With the static ΔE₀ already provided (S1), A adds almost nothing for F0. Much of what A learns in S0 is therefore the static part.
   - The training fit is good (R² 0.54–0.81) while out-of-fold R² is 0.21–0.55. With these features and learners, the location correction is predicted only partly across formulas.
   - A's target is the residual averaged over 50–1500 K, not the local residual at the transition. A lower full-range residual is therefore not the same as better signs or crossings near the decision region.
2. **Adding a uniform frequency scale (B) does not help.**
   - It lowers δΔS somewhat (F0 8.2 → 7.0; MS 10.5 → 8.0–8.1).
   - But it worsens ordering (significantly for F0 S0, +3.7 points, and MS S1, +2.7), worsens 19–22 pairs by more than 10 meV/atom, and its location stage stops transferring (out-of-fold R² ≤ 0.12, down to −0.19).
   - Learning more response information made transfer worse.
3. **Spectral correction (C) fails.**
   - The band-matched frequency errors fit in training but transfer poorly: θ out-of-fold R² 0.11–0.29.
   - Re-integrating the predicted correction worsens δΔS for F0 (8.2 → 10.6), ordering (+1.5 to +3.3 points) and the tail (23–30 pairs worse by more than 10 meV/atom).
   - With these features, the phase-specific softening seen in the tail audit cannot be predicted across formulas.
4. **Model combination (D) gives the clearest gain, and the plain average does best.**
   - Wrong ordering falls from 11.4% (F0) and 10.6% (MS) to 8.7% (S0), and from 7.9% / 5.3% to 3.9% (S1). The gain against F0 has an interval excluding 0 in both settings.
   - Tails shrink: p95 of \|e₁₀₀₀\| falls from 17–19 to 12.8 meV/atom, and pairs with \|e\| > 20 from 22–25 to 14.
   - Conformal coverage rises from 0.31–0.32 to 0.46 at 1000 K.
   - Learned weights (D2, D3) do not beat w = 0.5.
   - Against F0, the near-transition error does not improve, and 10–12 pairs get worse by more than 10 meV/atom. That is above the protocol's 7-pair limit, so D is "partial" against F0 and "success" against MatterSim.
   - **Interpretation:** the two models' errors are complementary in some cases, so averaging dilutes part of each model's tail. This does not require the errors to be independent.
   - **It also imports the other model's tail.** The average is worse than F0 alone, by more than 10 meV/atom, on 12 pairs in S0 and 10 in S1. Almost all are the halide pairs where MatterSim softens one phase: AgI, LiI, CsCl, CsBr, CuCl, Ba₂SrI₆.
   - **No validation from the repeats.** D1 gives identical results in every repeat only because it has no fitted parameter. That is not extra validation.
5. **None of the routes approaches the root-finding budget.** At most 21% of transition pairs are within ±50 K, against 12–18% for "none".

## Attribution (as the protocol requires)

| Route | Label | Cause |
|---|---|---|
| A | partial | Limited out-of-fold accuracy with these features. The target is a full-range mean, not the decision region. It adds false crossings. |
| B | failure | Limited out-of-fold accuracy, plus tails: the added response flexibility overfits. The δΔS norm improves while ordering worsens. |
| C | failure | Limited out-of-fold accuracy, plus tails. The spectral target is not predicted from these model and structure features, and it is not the decision quantity. |
| D | success against MS, partial against F0 | The combination helps ordering and tails without fitting. It does not reduce the near-transition error relative to F0, and it brings MatterSim's tail into 10–12 pairs. |

**Answer to the routing question (qualified).** The current features and learners learn part of each correction, but their prediction across formulas is not accurate enough for the transition requirement. This does **not** show that the corrections are inherently non-transferable. The results mix three causes, and they cannot be assigned to one:
- **features** that do not carry the needed information;
- **target mismatch**: A targets the full-range mean residual and C a spectral target, neither of which is the decision-region quantity;
- **limited generalisation.**

Each route improves a different capability:
- A lowers the near-transition RMS but not the ordering, and it adds false crossings.
- B lowers δΔS but worsens ordering.
- D improves ordering and tails but not the near-transition error relative to F0.

No single success/failure label summarises them.

## Closure for this data set (2026-11-02)

**Method search on these data stops here.** A–C are reported as diagnostic interventions. They show that a lower aggregate residual is not enough for the transition task; they are not proposed methods.

**Benchmark default combination baseline (frozen): D1, equal-weight averaging of the two models' G(T)**, one weight of 0.5 for every formula.
- It is reported **alongside** the single models. It is not presented as a new potential or as a method guaranteed to be more accurate.
- Development standing, under the frozen criteria: **success against MatterSim, partial against F0.** It improves ordering in both settings, with intervals excluding 0 against F0. It does not improve the near-transition error against F0, and it worsens 10–12 individual pairs against F0 (listed above) and 3–7 against MatterSim.

| Capability (S0 / S1) | F0 | MatterSim | D1 average |
|---|---|---|---|
| Near-transition RMS e(T_c), meV/atom | 7.09 / 4.95 | 9.45 / 7.37 | 6.92 / 4.95 |
| Within the ±50 K budget | 0.12 / 0.16 | 0.07 / 0.18 | 0.11 / 0.19 |
| Wrong ordering, mean over 600/1000/1500 K | 11.4 / 7.9% | 10.6 / 5.3% | 8.7 / 3.9% |
| p95 \|e(1000 K)\|, meV/atom | 17.2 / 14.9 | 18.6 / 13.9 | 12.8 / 9.1 |
| Pairs with max \|e\| > 20 meV/atom | 25 / 26 | 22 / 18 | 14 / 10 |
| Pairs worse than the single model by > 10 meV/atom | — | — | 12 / 10 against F0; 7 / 3 against MatterSim |
| Conformal coverage at 1000 K (α = 0.05) | 0.32 / 0.36 | 0.31 / 0.40 | 0.46 / 0.57 |
| Inference cost on the 521 entries (stored per-entry wall time: forces plus phonopy) | 1.33 h | 1.60 h | 2.93 h (2.2× F0, 1.8× MatterSim) |

The wall times come from the original runs on the shared A100, with different worker counts, so they are approximate.

**Entry condition for the next phase: independent confirmation.**
- Neither the D1 baseline nor any route claim has been confirmed. Confirmation needs a fresh, tier-compatible formula set, with the baseline and metrics frozen in advance.
- None was found among the sources searched in this round (`fes-confirmation-screen-20261102/INVENTORY.md`), so this gap is left open and stated.
- Adding models, changing losses or searching further on these data are not substitutes for confirmation.
- No paid DFT is started.

**Limits.**
- All four routes were evaluated on data already analysed.
- Only one structure dataset and two models were used.
- D needs two model evaluations.
- Any confirmation needs a fresh formula set, with the routes and their variants frozen in advance.
