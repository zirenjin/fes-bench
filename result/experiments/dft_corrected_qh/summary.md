# DFT-corrected QH static-energy replacement

The canonical `F_QH` static term was audited phase-by-phase before applying the correction. For SiO2 it is `F_QH - F_vib`; for metal QH files the row-wise `E_static_eV_per_atom` column is used. The correction therefore subtracts exactly the static term contained in `F_QH`, avoiding double counting.
The corrected curve is `G_corr(T) = F_QH(T) - E_DPA + E_DFT`, with both sides of every pair treated identically. The five-pair summaries exclude cristobalite--tridymite because it has no reference crossing; pair-level output retains the five crossing pairs with fold provenance.

| Split | Predictor | ΔG MAE (eV/atom) | balanced sign accuracy | Tc error | mean |Tc error| (K) | false | missed |
|---|---|---:|---:|---|---:|---:|---:|
| temp_extrap | T1_E_plus_F_QH | 0.14943348952911875 | 0.5 | n/a:missed_crossing | n/a:missed_crossing | 0 | 5 |
| temp_extrap | DFT_corrected_QH | 0.1609811916510002 | 0.5 | n/a:missed_crossing | n/a:missed_crossing | 0 | 5 |
| phase_lopo | T1_E_plus_F_QH | 0.113585483749105 | 0.5 | n/a:missed_crossing | n/a:missed_crossing | 0 | 10 |
| phase_lopo | DFT_corrected_QH | 0.12254590437368973 | 0.5 | n/a:missed_crossing | n/a:missed_crossing | 0 | 10 |

Maximum checked static-term difference relative to metadata E_DPA = 3.194e-03 eV/atom. The correction uses the exact QH-contained static term, not metadata E_DPA.

DFT convergence status: all nine adopted static outputs are `final-adopted-converged`; the Hf neighboring-k-mesh deviation remains documented in the DFT static convergence report.
