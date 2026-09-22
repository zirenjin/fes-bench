"""Create E1 constant-floor skill scores and ΔG-MAE convention audit."""

from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path

import numpy as np

from fes_bench.data.load import load


PHASES = ("cristobalite_beta", "quartz_beta", "tridymite_p63mmc")
PAIRS = (("cristobalite_beta", "quartz_beta"), ("cristobalite_beta", "tridymite_p63mmc"), ("quartz_beta", "tridymite_p63mmc"))
PAIR_NAMES = ("crist–quartz", "crist–trid", "quartz–trid")
TRAIN_END = 592  # frozen temp_extrap: indices 0..591 / 851–1442 K


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--curves", required=True, help="Directory containing the ten formal-window curve JSON files")
    parser.add_argument("--data-root", default="data")
    parser.add_argument("--output-root", default="results/trivial_floor")
    parser.add_argument("--audit-output", default="results/phase5/mae_convention_audit.md")
    return parser


def _mean_std(values: list[float]) -> str:
    return f"{np.mean(values):.7f} ± {np.std(values):.7f}"


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    data_root = Path(args.data_root)
    reference = {phase: np.array([point.G_eV_per_atom for point in load("sio2", phase, data_root).G_table]) for phase in PHASES}
    ref_delta = [reference[left] - reference[right] for left, right in PAIRS]
    constants = [float(delta[:TRAIN_END].mean()) for delta in ref_delta]
    baseline_full = [float(np.mean(np.abs(delta - constant))) for delta, constant in zip(ref_delta, constants)]
    rows: list[dict[str, object]] = []
    for path_string in sorted(glob.glob(str(Path(args.curves) / "*_formal_train_window.json"))):
        path = Path(path_string)
        basis, seed = path.stem.removesuffix("_formal_train_window").rsplit("_seed", 1)
        payload = json.loads(path.read_text(encoding="utf-8"))["predictions"]
        prediction = {phase: np.array(payload[f"sio2:{phase}"]["G_eV_per_atom"]) for phase in PHASES}
        model_delta = [prediction[left] - prediction[right] for left, right in PAIRS]
        pair_mae = [float(np.mean(np.abs(predicted - observed))) for predicted, observed in zip(model_delta, ref_delta)]
        train_mae = [float(np.mean(np.abs(predicted[:TRAIN_END] - observed[:TRAIN_END]))) for predicted, observed in zip(model_delta, ref_delta)]
        test_mae = [float(np.mean(np.abs(predicted[TRAIN_END:] - observed[TRAIN_END:]))) for predicted, observed in zip(model_delta, ref_delta)]
        rows.append({"basis": basis, "seed": int(seed), "pair_mae": pair_mae, "train_mae": train_mae, "test_mae": test_mae})
    out = Path(args.output_root); out.mkdir(parents=True, exist_ok=True)
    lines = ["# E1 skill scores against the training-mean constant ΔG floor", "", "For each SiO2 pair, `constant_delta_g` is the mean reference ΔG over frozen temp-extrap training indices 0–591 (851–1442 K), then evaluated on the full 851–2499 K grid. `skill = 1 − MAE_model / MAE_constant`. The all-pair aggregate applies that formula to the mean pair MAE; it is not the arithmetic mean of pairwise skills.", "", "| Pair | training mean ΔG (eV/atom) | constant ΔG MAE (eV/atom) |", "|---|---:|---:|"]
    for name, constant, mae in zip(PAIR_NAMES, constants, baseline_full): lines.append(f"| {name} | {constant:.8f} | {mae:.8f} |")
    lines.extend(["", "## Per-checkpoint scores", "", "`crist–trid` is a no-reference-crossing pair with a 0.0727 meV/atom floor; it is shown, but its scale makes individual-pair skills highly unstable. **Every aggregate score below is negative**, meaning the checkpoint has larger full-grid ΔG MAE than the corresponding constant floor.", "", "| basis | seed | all pairs: model / constant MAE (meV), skill | excluding crist–trid: model / constant MAE (meV), skill | crist–quartz skill | crist–trid skill | quartz–trid skill |", "|---|---:|---|---|---:|---:|---:|"])
    for row in rows:
        mae = row["pair_mae"]
        skill_all = 1 - float(np.mean(mae)) / float(np.mean(baseline_full))
        skill_no = 1 - float(np.mean([mae[0], mae[2]])) / float(np.mean([baseline_full[0], baseline_full[2]]))
        individual = [1 - value / base for value, base in zip(mae, baseline_full)]
        lines.append(f"| {row['basis']} | {row['seed']} | {np.mean(mae)*1000:.3f} / {np.mean(baseline_full)*1000:.3f}, **{skill_all:.3f}** | {np.mean([mae[0], mae[2]])*1000:.3f} / {np.mean([baseline_full[0], baseline_full[2]])*1000:.3f}, **{skill_no:.3f}** | {individual[0]:.3f} | {individual[1]:.3f} | {individual[2]:.3f} |")
    lines.extend(["", "## Aggregation sensitivity", "", "| basis | full-grid model MAE, all pairs (meV/atom) | excluding crist–trid (meV/atom) |", "|---|---:|---:|"])
    for basis in ("polynomial", "tlog_polynomial"):
        selected = [row for row in rows if row["basis"] == basis]
        all_values = [float(np.mean(row["pair_mae"])) * 1000 for row in selected]
        no_values = [float(np.mean([row["pair_mae"][0], row["pair_mae"][2]])) * 1000 for row in selected]
        lines.append(f"| {basis} | {_mean_std(all_values)} | {_mean_std(no_values)} |")
    lines.extend(["", "Recommendation: report the no-reference-crossing-pair aggregate as the primary ΔG MAE, with the all-pair value in a footnote/supplement. The crist–trid pair remains useful for false-crossing diagnostics, but its 0.059 meV/atom reference standard deviation (and 0.073 meV/atom constant MAE) otherwise gives it disproportionate leverage.", "", "## Deviations from design", "", "No new model was trained. The E1 historical checkpoints are evaluated on their original full grid; the only fitted quantity is the explicitly requested pairwise training-mean constant.", ""])
    (out / "skill_scores.md").write_text("\n".join(lines), encoding="utf-8")

    def aggregate(basis: str, field: str) -> float:
        selected = [row for row in rows if row["basis"] == basis]
        return float(np.mean([np.mean(row[field]) for row in selected]))
    audit = ["# ΔG MAE convention audit", "", "## Resolution", "", "The reported **2.56 meV/atom is not an in-domain aggregate**. It is the single `quartz_beta__tridymite_c2221` cell (`0.002562650 eV/atom`) in the historical four-phase file `fes_static_only_4phase_extensive_gauge_c_a100/all_pair_metrics.json`, evaluated on 1649 full-grid points. It therefore cannot be compared to the current three-phase E1 aggregate as an ‘ID MAE’. The apparent contradiction is a mislabeled/detached table cell plus changed phase set/checkpoint, not a train-vs-test, pair-vs-phase, or grid convention difference.", "", "## Matched re-evaluation", "", "Current E1 uses the same ten historic four-phase checkpoints but restricts evaluation to quartz, cristobalite, and P6₃/mmc tridymite, on the original 851–2499 K / 1649-point grid. Its reported aggregate is the equal-weight mean of the three pairwise full-grid MAEs.", "", "| basis | full grid (meV/atom) | frozen train window 851–1442 K (meV/atom) | test window 1443–2499 K (meV/atom) |", "|---|---:|---:|---:|"]
    for basis in ("polynomial", "tlog_polynomial"):
        audit.append(f"| {basis} | {_mean_std([float(np.mean(row['pair_mae']))*1000 for row in rows if row['basis']==basis])} | {_mean_std([float(np.mean(row['train_mae']))*1000 for row in rows if row['basis']==basis])} | {_mean_std([float(np.mean(row['test_mae']))*1000 for row in rows if row['basis']==basis])} |")
    audit.extend(["", "These recomputed full-grid values reproduce E1 (polynomial ≈12.826 meV/atom; T-log ≈6.515 meV/atom). Neither matched train-window value (≈18.13 and ≈9.83 meV/atom) is 2.56 meV/atom. ΔG cancels the system gauge, so the formal-train versus heldout calibration convention also cannot change these ΔG MAEs.", "", "## Pair-set effect", "", "The historic 2.56265 meV value belongs to a pair containing C2221, whereas E1 deliberately excludes C2221 because its representative structure was not comparable. It also uses a different historical gauge-C run, rather than either of the polynomial/T-log checkpoint families. Thus the values answer different experimental questions and must not appear together in a common ΔG-MAE comparison.", "", "## Required reporting convention", "", "Use: ‘equal-weight mean of pair-level ΔG MAE, full 851–2499 K grid, three-phase subset, with no-crossing crist–trid reported separately (and excluded from the primary aggregate).’ Never label the isolated C2221-pair number as ID aggregate MAE.", "", "## Deviations from design", "", "No new model was trained; this is a read-only re-evaluation of frozen checkpoint curves and a read-only inspection of the archived historical metrics file.", ""])
    audit_path = Path(args.audit_output); audit_path.parent.mkdir(parents=True, exist_ok=True); audit_path.write_text("\n".join(audit), encoding="utf-8")
    print(out / "skill_scores.md"); print(audit_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
