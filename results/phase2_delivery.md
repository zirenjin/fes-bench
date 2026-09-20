# Phase 2 delivery — Hf QH gate

Status: the hcp gate is unlocked for `Domains_Alloy`; bcc remains held because
its spectrum is unstable.  No QH was expanded to other systems.

## Files changed or added

- `fes_bench/qh/run.py`: configuration-driven finite-displacement QH runner.
- `fes_bench/qh/compare.py`: one-system calibration and reference comparison.
- `fes_bench/qh/materialize.py`: provenance-preserving conversion of the
  completed diagnostic files into the canonical `data/<system>/<phase>/fqh.csv`
  schema.  It requires an explicit reliability decision and never upgrades a
  diagnostic curve.
- `configs/qh/hf_qh.yaml`, `hf_compare.yaml`,
  `hf_qh_domains_alloy.yaml`, `hf_compare_domains_alloy.yaml`, and
  `hf_materialize_domains_alloy.yaml`.
- `results/phase2/hf/` and `results/phase2_domains_alloy/hf/`: QH tables,
  phonon summaries, comparison JSON, PDF/PNG figures, and run logs.
- `results/phase2/hcp_imaginary_diagnosis.md`: Gamma-centered reciprocal-space
  diagnosis and explicit head decision.
- `results/phase2/hcp_domains_alloy_rerun_20260921.log`: single-point sanity
  plus the fresh isolated hcp QH rerun.

## Actual commands and key results

Run in the isolated thu-GenSi benchmark directory with the `mat-agent`
environment:

```bash
PYTHONPATH=. conda run -n mat-agent python -m fes_bench.qh.run \
  --config configs/qh/hf_qh_domains_alloy.yaml
PYTHONPATH=. conda run -n mat-agent python -m fes_bench.qh.compare \
  --config configs/qh/hf_compare_domains_alloy.yaml
```

The formal frozen-representation head, `DPA-3.1-3M.pt` + `Domains_Alloy`,
gives a stable Hf hcp harmonic spectrum (minimum frequency 0.230--0.242 THz)
but an unstable bcc spectrum (minimum frequency -1.989 to -1.906 THz).  On the
T > 1800 K evaluation region after the single system-level calibration,
the Hf QH baseline MAE is 208.642 meV/atom (hcp 92.595; bcc 324.690), with
maximum absolute error 408.710 meV/atom.

For diagnosis only, the generic `Domains_SSE_PBE` result is retained in
`results/phase2/hf/`: both hcp and bcc have larger imaginary modes and its
test MAE is 192.159 meV/atom.  It is not selected as the formal FES head.

## Artifacts

- `results/phase2_domains_alloy/hf/{hcp,bcc}_fqh.csv`
- `results/phase2_domains_alloy/hf/qh_summary.json`
- `results/phase2_domains_alloy/hf/qh_comparison.json`
- `results/phase2_domains_alloy/hf/fqh_vs_reference.{png,pdf}`
- `results/phase2/hf/` contains the same diagnostic artifacts for
  `Domains_SSE_PBE`.
- `data/hf/{hcp,bcc}/fqh.csv` and `data/hf/{hcp,bcc}/phonon_report.json`:
  canonicalized copies of the selected `Domains_Alloy` diagnostic.  Both
  reports set `qh_reliable: true` for hcp and `false` for bcc, with
  `diagnostic_only: true`; only hcp may pass the current reliability gate.
  Local and thu-GenSi execution
  logs are `results/phase2_materialization_20260919.local.log`,
  `results/phase2_materialization_20260920_thu.log`, and
  `results/phase2_materialization_20260921_v100_r17.log`. The V100 repeat
  passed the focused and full regression checks.

## Acceptance check

- Hf two-phase QH files and a visual `fqh_vs_reference` comparison: **pass**.
- Smooth, tens-of-meV `E_DPA + F_QH + c_system` agreement: **fail**; residual
  errors are about 0.1--0.4 eV/atom.
- Harmonic stability gate: **hcp pass, bcc fail**; the selected hcp head has
  no mode below -0.05 THz, while bcc retains imaginary modes.
- Expansion to SiO2, CaSiO3, Ti, and Zr: **not run**, as prescribed by the
  failed gate rather than consuming additional QH compute.

## Deviations from design

- The implementation currently uses three volume scales (0.98, 1.00, 1.02),
  a 2x2x2 supercell, and a 12x12x12 mesh.  This is a diagnostic prototype,
  not the final five-volume, >=10-A-supercell, >=15/A-mesh protocol required
  by the design.
- The runner records minimum frequencies but does not yet write the required
  imaginary-mode fraction, Vinet-fit parameters, cached force constants, or
  per-phase timing report.  These are intentionally not represented as
  complete Phase 2 artifacts.
- The canonical Hf `fqh.csv` files remain diagnostic three-volume materialized
  curves (no Vinet fit is claimed).  The hcp report is now reliable only after
  the recorded scattered-mode diagnosis, single-point sanity, fresh rerun, and
  smooth-deviation check; bcc remains unavailable to production predictors.
- `Domains_SSE_PBE` was run only as a rejected diagnostic comparison; the
  selected FES representation is `Domains_Alloy`.  The observed bcc
  instability and large residual slope mismatch mean neither diagnostic
  result is valid as the final physical QH baseline.
