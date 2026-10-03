# DFT-corrected QH static-energy replacement

The canonical `F_QH` contains a constant static term plus `F_vib`; this was checked phase-by-phase before applying the correction.
The corrected curve is `G_corr(T) = F_QH(T) - E_DPA_in_F_QH + E_DFT`, which is the requested `F_QH - E_DPA + E_DFT` with the actual static term represented by the QH file. The provenance `E_DPA` is reported separately where it differs by about 1 meV/atom.

| Split | Predictor | ΔG MAE (eV/atom) | sign accuracy | mean |Tc error| (K) | false | missed |
|---|---|---:|---:|---:|---:|---:|
| temp_extrap | T1_E_plus_F_QH | 0.14943348952911875 | 0.4557863953791643 | n/a | 0 | 5 |
| temp_extrap | DFT_corrected_QH | 0.1609811916510002 | 0.146609479579732 | n/a | 0 | 5 |
| phase_lopo | T1_E_plus_F_QH | 0.113585483749105 | 0.5293352836014952 | n/a | 0 | 10 |
| phase_lopo | DFT_corrected_QH | 0.12254590437368973 | 0.45620022753128553 | n/a | 0 | 10 |

Maximum checked |F_QH − F_vib − E_DPA| = 3.194e-03 eV/atom.
