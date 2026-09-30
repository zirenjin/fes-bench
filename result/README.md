# Result index

`result/` is the sole generated-evidence layer. Every table has a CSV for
review, a sibling `.meta.json` with commit and input hashes, and a producer in
`src/tables/`. Every experiment directory has `findings.csv`,
`findings.meta.json`, and a top-level `metrics.json` provenance index; raw
metrics remain below its `raw_runs/` directory when an experiment has them.

| Conclusion or output | Command | Source evidence |
| --- | --- | --- |
| system / phase / split / metric inventories | `python3 src/tables/{system,phase,split,metric}_inventory.py` | `experiments/data_prep/raw_inventory/` |
| predictor comparison and Tc errors | `python3 src/tables/{predictor_comparison,crossing_errors}.py --split <split>` | `experiments/external_baselines/raw_<split>/` |
| E1 / E2 / E3 | `python3 src/experiments/{crossing_reevaluation,calibration_window,synthetic_recovery}.py` | matching experiment directory |
| floor and skill audit | `python3 src/experiments/skill_floor.py` | `experiments/skill_floor/` |
| QH and imaginary-mode records | `python3 src/experiments/{quasi_harmonic,imaginary_modes}.py` | `experiments/quasi_harmonic/` |
| Ideal SiO₂ representative rebuild and QH comparison | `PYTHONPATH=src/lib python src/experiments/data_prep/rebuild_sio2_representatives.py --config configs/representatives/sio2_ideal_rebuild.yaml --relax --checkpoint <external>`; `python src/experiments/data_prep/compare_qh_rebuild_all.py` | `experiments/quasi_harmonic_rebuild/` and `data/processed/sio2/_archive_md_snapshot/` |
| Metal ≥10 Å QH comparison | `python src/experiments/data_prep/compare_qh_rebuild_all.py` (QH production commands are the three `configs/qh/*_qh_10a.yaml` configs) | `experiments/quasi_harmonic_10a/` |
| MAE, split-sign, representative, reference, leakage audits | `python3 src/experiments/<name>.py` | matching experiment directory |

All paths are relative to the repository root. Rebuild the complete derived
layer with `bash src/reproduce_all.sh`.

The Tc-error value for `interp_const` in
`tables/crossing_errors_temp_extrap.csv` comes from
`experiments/external_baselines/raw_temp_extrap/interp_const/seed_none/metrics.json`.
That file's `provenance` block records the imported source hash and frozen split
hash; `src/tables/crossing_errors.py --split temp_extrap` performs the read.

## N/A vocabulary

- `n/a:no_reference_crossing`: the frozen reference grid has no sign change.
- `n/a:pair_only_predictor`: the predictor does not produce phase-level `G`.
- `n/a:no_training_phase`: a held-out fold lacks a required training phase or system.
- `n/a:degenerate_prediction`: predicted ΔG is within `1e-12 eV/atom` of zero on the entire evaluation grid, so crossings and their derived errors are undefined.
- `n/a:input_unavailable`: a predictor's required immutable input is absent.

The matching `.meta.json` files retain field-specific explanations. The current
tables have 33, 15, and 32 cells of these respective types.

## E1 reproduction tolerance

The fixed policy is in `configs/reproduction_tolerance.yaml`; it is derived
from resolution, not fitted to a reproduction result. Crossing counts, false/
missed crossings, list cardinalities, `sign_accuracy`, and other ratios must
match exactly. `T_c` and `T_c`-error fields use an absolute tolerance of
0.05 K (one twentieth of the 1 K grid); G/ΔG and calibration-energy fields use
1e-6 eV/atom (one thousandth of the reported meV/atom precision). GPU
floating-point evaluation is not required to be bitwise identical across
runtime environments.

The canonical E1 checkpoint output is
`experiments/crossing_reevaluation/sio2_reanalysis.json`, produced by
`src/experiments/crossing_reevaluation_run.py`. Its field-level audit is
`experiments/crossing_reevaluation/reproduction_check.csv`; the canonical and
historical SHA-256 values, and the pass result, are recorded in its sibling
`sio2_reanalysis.meta.json`. The original input remains unchanged at
`experiments/legacy_support/phase5/sio2_reanalysis.json`.

E1 status: complete; acceptance criterion three did not pass because the
cristobalite–tridymite pair has one false crossing for each of 10/10
checkpoints. This is a model conclusion, not a reproduction failure.

## Known issue recorded, not resolved

Skill scores now use the zero floor uniformly on every split and subset:
`1 − MAE / MAE_zero`, where `MAE_zero` is the mean reference `|ΔG|` on the
evaluated points. The zero floor is label-free, available for every split, and
avoids making the denominator depend on whether a fold happens to contain a
training phase.

Bartel covers Hf from the Domains_Alloy-computed Hf E0. Its Hf error is much
larger than its SiO₂ error; the SiO₂-only table gives 28.2 meV/atom.

无训练 predictor 在 phase_lopo 与 system_loso 上指标相同是预期的：两者的测试点并集都是全部相对的全温区，权重比例一致。

temp_extrap 的无训练 predictor 只在共享的 `T>T*` 测试点评测；zero 排除 crist–trid 后为 3.4692 meV/atom（全网格为 5.4372），按测试点评测后与参考口径一致。

`_legacy/` contains six cleaned historical delivery records and `CLEANING.md`.
They are tracked for provenance but not read by table scripts.

## Table 2b structure/QH audit

The SiO₂ representatives were rebuilt from AFLOW prototype CIFs (the DaRUS-4999/
Forslund phase models agree: β-quartz P6₄22, β-cristobalite Fd-3m, and
β-tridymite P6₃/mmc) and hydrostatic variable-cell relaxed with the
`Domains_SSE_PBE` head under `FixSymmetry`; all three final space
groups match their target groups. The previous MD snapshots remain in
`data/processed/sio2/_archive_md_snapshot/`. The structured old/new comparison
is `experiments/quasi_harmonic_rebuild/findings.csv`; Bartel, downstream QH
comparisons, imaginary-mode summaries, and reference statistics are explicitly
marked in each phase `meta.json` as requiring recomputation and were not
silently regenerated from the new representatives.

QH now uses the same ≥10 Å shortest-edge protocol for SiO₂ and metals: SiO₂
uses `[2,2,2]` for all three phases, while Hf/Ti/Zr use phase-specific cells.
The canonical SiO₂ QH outputs and imaginary fractions are in
`experiments/quasi_harmonic_sio2_domains_alloy/`; superseded SSE-PBE versions
are archived under `experiments/legacy_support/qh_superseded_20260930/`; metal old/new comparisons are in
`experiments/quasi_harmonic_10a/{hf,ti,zr}/findings.csv`. The superseded
fixed-cell SiO₂ QH outputs are retained under the explicitly invalid
`legacy_invalid_fixed_cell_20260929/` directory and are not used by tables.

The new-metal imaginary-mode fractions (old `[2,2,2]` → new ≥10 Å cell) are:
Hf hcp `0 → 0`, Hf bcc `0.199074 → 0.210784`; Ti hcp `0 → 0`, Ti bcc
`0.273810 → 0.264881`; Zr hcp `0 → 0`, Zr bcc `0.182540 → 0.199405`.
The reliability decisions are unchanged (hcp true, bcc false).

Split v2 is frozen under `data/processed/splits_v2/`: `temp_extrap.json` and
`phase_lopo.json` include Hf/SiO₂/Ti/Zr, while `system_loso.json` is the three
metal folds (`hf`, `ti`, `zr`) only. Its generated SHA-256 is stored in each
JSON; v1 files under `data/processed/splits/` remain unchanged.

Formal SiO₂ E_DPA and F_QH provenance now uses `DPA-3.1-3M.pt` with
`Domains_Alloy`; the non-formal `Domains_SSE_PBE` relaxation/QH diagnostic is
retained separately. For β-quartz, the dedicated soft-mode run found a minimum
of −1.153 THz at Γ (20 negative modes across 19 q points; only 1/20 negative
modes inside the configured Γ neighborhood). It is annotated as a possible
physical soft mode, but is not Γ-exclusive; the qh_reliable rule is unchanged.
The 1% imaginary-mode fraction threshold may be insensitive to soft-mode-driven
phase transitions.

`phase_inventory.csv` takes every SiO₂ imaginary fraction from
`experiments/quasi_harmonic_sio2_domains_alloy/raw_runs/sio2/qh_summary.json`,
the `1.0` equilibrium-volume entry with `[2,2,2]` supercell; the superseded
SSE-PBE QH outputs are archived under
`experiments/legacy_support/qh_superseded_20260930/`. The dedicated β-quartz
soft-mode diagnostic reports a minimum of −1.1533986 THz at q = `[0,0,0]`
(Γ), recorded in the phase table without replacing the formal fraction.
The previously reported 0.161% / 4.75% / 5.57% values came from the
SSE-PBE diagnostic and are now archive-only, so they no longer compete with
the canonical table source.

Checkpoint provenance for E_DPA/F_QH artifacts is indexed by
`experiments/checkpoint_consistency/findings.csv`. Formal fields use SHA-256
`86dd3a804d78ca5d203ebf98747e8f16dff9713ba8950097ceb760b161e19907` and
`Domains_Alloy`; `Domains_SSE_PBE` is retained only for explicitly labelled PES
relaxation/soft-mode diagnostics. Bartel, QH, imaginary-mode, and E4
before/after changes are under `tables/_changes/`.
