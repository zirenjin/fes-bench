# Overnight batch status — 2026-09-22

This is a line-by-line queue report. Lines A–D were started independently; the live QH line does not block the completed baseline and survey work.

## Line A — external baselines: complete (structural N/A folds recorded)

Commands run:

```text
PYTHONPATH=. python3 scripts/run_external_baselines.py --data-root data --output-root results/external_baselines --splits splits/temp_extrap.json splits/phase_lopo.json splits/system_loso.json
```

Artifacts: `results/external_baselines/metrics.json`, `skill_scores.md`, `README.md`; implementation in `fes_bench/baselines/` and `scripts/run_external_baselines.py`.

Key result: Bartel skill is negative wherever a matched global-mean pair floor exists: −2.848 on temp extrap, −2.950/−4.992/−3.384 on the SiO₂ LOPO folds, and −1.326 on SiO₂ LOSO. Interpolation is positive on temp extrap (+0.692). The mat-agent rerun completed; phase-ID MLP skill is also negative on temp extrap (−0.858), while LOPO/LOSO are structurally unavailable because the held-out target has no training rows. `global_mean_delta_g` is pair-level and uses all available training pair labels, so LOPO/LOSO retain a nontrivial floor. The separate Hf DPA E0/volume helper hit a PyTorch/e3nn safe-load compatibility error and was terminated; Bartel Hf rows remain explicitly unavailable.

Bartel coefficient audit: original and implementation values match for all four coefficients. The fixed-composition limitation (mass identical, volume-only phase discrimination) is in the docstring and report.

## Line B — SiO₂ QH: complete

Command launched in isolated `/share/jzr/codex_fes_overnight_20260922_1790052610710037915`:

```text
/root/miniconda3/envs/mat-agent/bin/python -m fes_bench.qh.run --config configs/qh/sio2_qh_domains_sse_pbe.yaml
```

The raw QH process completed after roughly seven hours. It produced all three 1,649-point `fqh.csv` files and `qh_summary.json` under the isolated workspace; `qh.compare` produced `qh_comparison.json` and `fqh_vs_reference.{png,pdf}` with all-phase test MAE 21.8657 meV/atom. The serial same-protocol diagnoses also completed: all three phases have non-Γ-local negative modes and `qh_reliable=false` (details in `results/phase2/sio2_qh_report.md` and the remote JSON evidence).

## Line C — environment: partial

The curated code/data were synchronized to the isolated remote workspace; the unrelated pre-existing `/share/jzr/fes-bench` directory was not used. The direct V100 SSH probe timed out again at `dcwq1547908.bohrium.tech:22` at 15:35 UTC, and thu-GenSi cannot resolve `bohrium-v100`. Authenticated Bohr resource discovery lists V100 SKUs, but reports a 0 CNY balance, so no paid replacement job was submitted; V100 synchronization remains blocked by network/billing state. `compileall` passed remotely. The targeted baseline smoke process was launched and then terminated after it continued consuming CPU; the complete baseline result is already persisted locally. Pytest was initially absent, so a vendored pytest installation was made only under the isolated remote workspace; the corrected full run passed all 21 tests in 4750.31 s (exit 0), recorded in `results/phase6_readiness/pytest_vendor_full.log`.

## Line D — complete for lightweight items

Artifacts: `results/STATUS.md`, `results/phase5/sio2_second_truth_source_survey.md`, `results/table2_delta_g_amplitude.md`. The Forslund/DaRUS source is available as a 184.2 MB record with direct-upsampling thermodynamic tables and a separate functional-evaluation ZIP; it was not downloaded or processed.

## Deviations from design

No new training, no head-policy change, and no frozen-split modification. The remote QH runner uses the existing QH implementation; its final reliability/imaginary-mode status is not asserted until the process exits. The external MLP and full pytest suite remain environment holds with concrete missing library errors.
