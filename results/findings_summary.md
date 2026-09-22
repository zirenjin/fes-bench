# Findings summary

- Bartel 2018 skill relative to `global_mean_delta_g` is −2.848 on `temp_extrap`, −2.950/−4.992/−3.384 on the three SiO₂ `phase_lopo` folds, and −1.326 on SiO₂ `system_loso`; evidence: `results/external_baselines/skill_scores.md`.
- Across all ten E1 checkpoints, polynomial all-pair skill is −3.412 to −3.440 and `tlog_polynomial` all-pair skill is −1.245 to −1.251, with every aggregate negative; evidence: `results/trivial_floor/skill_scores.md`.
- `interp_const` has temp-extrap matched-pair MAE 0.001969 eV/atom versus the 0.00638272 eV/atom global floor and skill +0.692; evidence: `results/external_baselines/skill_scores.md`.
- Using the held-out calibration window reduces held-out scalar G MAE by 51.95–52.00% for polynomial and 59.55–59.59% for `tlog_polynomial`, while ΔG MAE changes remain at approximately 10⁻¹⁷ eV/atom and Tc changes at approximately 10⁻¹⁰ K; evidence: `results/phase5/e2_leakage_audit.md`.
- SiO₂ QH minimum frequencies are −2.63652/−1.70034/−1.65362 THz for quartz/cristobalite/tridymite, with negative-mode fractions 0.005467/0.013695/0.014559 and `qh_reliable=false` for all three; evidence: `results/phase2/sio2_qh_report.md`.
- Reference ΔG standard deviations are 3.607 meV/atom for cristobalite–quartz, 0.059 meV/atom for cristobalite–tridymite, 3.576 meV/atom for quartz–tridymite, and 11.987 meV/atom for Hf bcc–hcp; evidence: `results/table2_delta_g_amplitude.md`.
