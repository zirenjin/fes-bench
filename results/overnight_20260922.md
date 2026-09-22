# Overnight batch status — 2026-09-22

This is a line-by-line queue report. Lines A–D were started independently; the live QH line does not block the completed baseline and survey work.

## Line A — external baselines: partial, core results complete

Commands run:

```text
PYTHONPATH=. python3 scripts/run_external_baselines.py --data-root data --output-root results/external_baselines --splits splits/temp_extrap.json splits/phase_lopo.json splits/system_loso.json
```

Artifacts: `results/external_baselines/metrics.json`, `skill_scores.md`, `README.md`; implementation in `fes_bench/baselines/` and `scripts/run_external_baselines.py`.

Key result: Bartel skill is negative wherever a matched global-mean pair floor exists: −2.848 on temp extrap, −2.950/−4.992/−3.384 on the SiO₂ LOPO folds, and −1.326 on SiO₂ LOSO. Interpolation is positive on temp extrap (+0.692). `global_mean_delta_g` is pair-level and uses all available training pair labels, so LOPO/LOSO retain a nontrivial floor. A mat-agent rerun is live to obtain phase-ID MLP metrics with the working Torch environment. The separate Hf DPA E0/volume helper hit a PyTorch/e3nn safe-load compatibility error and was terminated; Bartel Hf rows remain explicitly unavailable.

Bartel coefficient audit: original and implementation values match for all four coefficients. The fixed-composition limitation (mass identical, volume-only phase discrimination) is in the docstring and report.

## Line B — SiO₂ QH: running

Command launched in isolated `/share/jzr/codex_fes_overnight_20260922_1790052610710037915`:

```text
/root/miniconda3/envs/mat-agent/bin/python -m fes_bench.qh.run --config configs/qh/sio2_qh_domains_sse_pbe.yaml
```

The process was verified live (PID 1996888; it has `/dev/nvidia4` and `/dev/nvidia-uvm` open; last poll 2026-09-22 06:20 UTC) with log at `results/phase2_sio2_sse_pbe.log`. The config requests `Domains_SSE_PBE`, three SiO₂ phases, five volume scales, 0.01 Å displacements, and 16³ mesh. Final `fqh.csv`, imaginary-mode reliability report, comparison plot, and failure-mode entries are pending process completion.

## Line C — environment: partial

The curated code/data were synchronized to the isolated remote workspace; the unrelated pre-existing `/share/jzr/fes-bench` directory was not used. The direct V100 SSH probe timed out at `dcwq1547908.bohrium.tech:22`, and thu-GenSi cannot resolve `bohrium-v100`, so V100 synchronization remains blocked by DNS/network state. `compileall` passed remotely. The targeted baseline smoke process was launched and then terminated after it continued consuming CPU; the complete baseline result is already persisted locally. Pytest was initially absent, so a vendored pytest installation was made only under the isolated remote workspace; the corrected full run is live as PID 2148593 with log `results/phase6_readiness_pytest_vendor_full.log`. See `results/phase6_readiness/`.

## Line D — complete for lightweight items

Artifacts: `results/STATUS.md`, `results/phase5/sio2_second_truth_source_survey.md`, `results/table2_delta_g_amplitude.md`. The Forslund/DaRUS source is available as a 184.2 MB record with direct-upsampling thermodynamic tables and a separate functional-evaluation ZIP; it was not downloaded or processed.

## Deviations from design

No new training, no head-policy change, and no frozen-split modification. The remote QH runner uses the existing QH implementation; its final reliability/imaginary-mode status is not asserted until the process exits. The external MLP and full pytest suite remain environment holds with concrete missing library errors.
