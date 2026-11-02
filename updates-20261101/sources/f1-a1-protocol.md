# Addendum A1: 2×2 anchor-calibration design and DFT-side mechanism strata

Frozen at 2026-11-01 ~22:20 UTC, while F1 training was at step 14,400 of 20,000. No F1 evaluation (`eval/`) existed.

F0 results on these pairs have been seen in the exploratory audit analyses. All F0-only quantities below are therefore descriptive. **The pre-specified part is the F1-vs-F0 contrast within each cell and stratum.** PROTOCOL.md and its primary endpoint are unchanged.

## A1.1 2×2 design (fixed same-formula pairs, not the stable-phase envelope)

For each pair, e(T) = ΔF_MLIP(T) − ΔF_DFT(T), with ΔF = E₀ + F_vib on the 0–1500 K grid (50 K steps).

| | No calibration | Constant (anchor) calibration at T₀ |
|---|---|---|
| F0 | e(T) | r(T) = e(T) − e(T₀) |
| F1 | e(T) | r(T) = e(T) − e(T₀) |

- **Anchor.** Primary T₀ = 300 K. Secondary T₀ = 0 K, where r(T) = δΔF_th(T), the main-protocol quantity.
- **Information budget.** The calibration uses the target pair's reference ΔF at T₀. That is one target thermodynamic label per pair, so these arms are anchored predictions, not zero-shot. F0 and F1 use the same anchor and the same budget.
- **Relative-entropy error.** δΔS(T) = −d e/dT, by central differences on the grid. It is reported in µeV/(atom K) at 600, 1000 and 1500 K. Constant calibration leaves it unchanged.

Metrics, on test formulas, at 600, 1000 and 1500 K:
- RMS of e and of r;
- RMS of δΔS;
- non-degenerate wrong-ordering rate (|ΔF_DFT| ≥ 1 meV/atom) for e and for r;
- F1/F0 RMS ratios with the paired formula-cluster bootstrap (2,000 draws).

**Pre-specified reading.** "Fine-tuning improves temperature dependence" is supported only if two things hold at 1000 K: r_F1 has a lower RMS than r_F0, and δΔS_F1 has a lower RMS than δΔS_F0, with both upper 95% ratio bounds below 1. A gain only in uncalibrated e, with no gain in r, is read as an anchor effect.

## A1.2 Mechanism strata (DFT-only features; no MLIP information)

Per pair:
- **|ΔV/V|.** Relative difference in DFT volume per atom: |V_i − V_j| / mean(V).
- **|ΔCN|.** Absolute difference in mean coordination number. CN of an atom = number of neighbours within 1.2 × that atom's nearest-neighbour distance, from the stored DFT unit cell under periodic boundaries.
- **|ΔS_DFT(1000 K)|.** Magnitude of the DFT vibrational-entropy difference, a proxy for "entropy-driven competition".

Each feature splits pairs at its median over all primary pairs (both splits), into small and large.

Per stratum and arm (F0, F1), at 1000 K:
- single-phase δF_th RMS;
- cancellation fraction c = 1 − ⟨δΔ²⟩ / ⟨δ_i² + δ_j²⟩ for δF_th;
- RMS of r (T₀ = 300 K) and of δΔS.

**Pre-specified F1 question.** Is the F1/F0 ratio for r and for δΔS different between the small and large strata? This is answered by the bootstrap CI of the difference of log ratios between strata.

The F0 stratum comparisons (does structural difference weaken cancellation?) are descriptive, because F0 has been seen.

## Not in scope

- Envelope (stable-phase) metrics stay those of PROTOCOL.md.
- Harmonic-validity and reference-discriminability strata are not available as DFT-only features here. The harmonic class is fixed by selection (all phases dynamically stable), and discriminability enters only through the |ΔF_DFT| ≥ 1 meV/atom ordering filter.
