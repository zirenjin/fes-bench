"""Evaluate a predictor specification against a frozen benchmark split."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

from fes_bench.data.load import ReferencePoint, load


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--predictor",
        required=True,
        help="reference, constant_offset:<eV>, or reference_iid_noise:<eV>",
    )
    parser.add_argument("--split", required=True, help="frozen split JSON")
    parser.add_argument("--data-root", required=True)
    parser.add_argument("--seeds", default="11,23,37,51,67")
    parser.add_argument("--output-root", default="results")
    return parser


def _root(delta: np.ndarray, temperatures: np.ndarray, merge_within_K: float = 0.0) -> list[float]:
    roots: list[float] = []
    index = 0
    while index < len(delta) - 1:
        left, right = delta[index], delta[index + 1]
        if left == 0:
            end = index
            while end + 1 < len(delta) and delta[end + 1] == 0:
                end += 1
            roots.append(float((temperatures[index] + temperatures[end]) / 2.0))
            index = end + 1
            continue
        if left * right < 0:
            fraction = abs(left) / (abs(left) + abs(right))
            roots.append(float(temperatures[index] + fraction * (temperatures[index + 1] - temperatures[index])))
        index += 1
    if merge_within_K <= 0 or len(roots) < 2:
        return roots
    clusters: list[list[float]] = [[roots[0]]]
    for value in roots[1:]:
        if value - clusters[-1][-1] <= merge_within_K:
            clusters[-1].append(value)
        else:
            clusters.append([value])
    return [float(np.mean(cluster)) for cluster in clusters]


def _predictor(spec: str, seed: int, key: str, reference: np.ndarray) -> np.ndarray:
    if spec == "reference":
        return reference.copy()
    prefix = "reference_noise:"
    if spec.startswith(prefix):
        sigma = float(spec[len(prefix) :])
        digest = hashlib.sha256(f"{seed}:{key}".encode()).digest()
        rng = np.random.default_rng(int.from_bytes(digest[:8], "little"))
        # One Gaussian energy-gauge perturbation per phase/seed preserves a
        # physically smooth G(T) curve, so Tc displacement can be compared to
        # MAE divided by the local crossing slope.
        return reference + rng.normal(0.0, sigma)
    prefix = "constant_offset:"
    if spec.startswith(prefix):
        return reference + float(spec[len(prefix) :])
    prefix = "reference_iid_noise:"
    if spec.startswith(prefix):
        sigma = float(spec[len(prefix) :])
        digest = hashlib.sha256(f"iid:{seed}:{key}".encode()).digest()
        rng = np.random.default_rng(int.from_bytes(digest[:8], "little"))
        return reference + rng.normal(0.0, sigma, size=reference.shape)
    raise ValueError(
        "predictor must be reference, constant_offset:<eV>, "
        "reference_noise:<eV>, or reference_iid_noise:<eV>"
    )


def _output_slug(spec: str) -> str:
    if spec.startswith("reference_iid_noise:"):
        return "root_stability_iid_" + spec.split(":", 1)[1]
    return spec.replace(":", "_")


def _table_cache(data_root: Path, rows: list[dict[str, object]]) -> dict[tuple[str, str], tuple[ReferencePoint, ...]]:
    keys = {(str(row["system"]), str(row["phase"])) for row in rows}
    return {key: load(key[0], key[1], data_root).G_table for key in keys}


def evaluate(spec: str, split: dict[str, object], data_root: Path, seeds: list[int]) -> dict[str, object]:
    if "folds" in split:
        raw_folds = split["folds"]
        if not isinstance(raw_folds, dict):
            raise ValueError("split folds must be an object")
        folds = {
            str(name): evaluate(spec, fold, data_root, seeds)
            for name, fold in raw_folds.items()
            if isinstance(fold, dict) and isinstance(fold.get("test"), list)
        }
        total_frames = sum(int(result["n_test_frames"]) for result in folds.values())
        weighted_mae = (
            sum(float(result["G_MAE_eV_per_atom"]) * int(result["n_test_frames"]) for result in folds.values())
            / total_frames
            if total_frames
            else math.nan
        )
        weighted_coverage = (
            sum(float(result["coverage_2sigma"]) * int(result["n_test_frames"]) for result in folds.values())
            / total_frames
            if total_frames
            else math.nan
        )
        return {
            "predictor": spec,
            "n_test_frames": total_frames,
            "seeds": seeds,
            "G_MAE_eV_per_atom": weighted_mae,
            "coverage_2sigma": weighted_coverage,
            "folds": folds,
        }
    test = split.get("test")
    if not isinstance(test, list):
        raise ValueError("split must provide a test list")
    train = split.get("train", [])
    if not isinstance(train, list):
        raise ValueError("split train entries must be a list")
    # A phase-LOPO fold evaluates only its held phase.  To calculate a physical
    # ΔG and Tc it must still be paired with predictions for the known phases
    # at the same state points.  Those partner frames remain outside scalar G
    # test metrics, but are explicitly included in pairwise curve metrics.
    pair_rows = test + train
    cache = _table_cache(data_root, pair_rows)
    by_system_phase: dict[tuple[str, str], list[int]] = defaultdict(list)
    test_by_system_phase: dict[tuple[str, str], list[int]] = defaultdict(list)
    for row in pair_rows:
        by_system_phase[(str(row["system"]), str(row["phase"]))].append(int(row["T_index"]))
    for row in test:
        test_by_system_phase[(str(row["system"]), str(row["phase"]))].append(int(row["T_index"]))
    point_errors: list[float] = []
    pair_metrics: dict[str, dict[str, object]] = {}
    coverage_values: list[bool] = []
    systems = sorted({system for system, _ in test_by_system_phase})
    for system in systems:
        phases = sorted(phase for item_system, phase in by_system_phase if item_system == system)
        predictions: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
        for phase in phases:
            indexes = sorted(by_system_phase[(system, phase)])
            reference = np.array([cache[(system, phase)][index].G_eV_per_atom for index in indexes])
            temperatures = np.array([cache[(system, phase)][index].T_K for index in indexes])
            if spec == "reference":
                seed_predictions = np.stack([reference.copy() for _ in seeds])
                mean, std = reference.copy(), np.zeros_like(reference)
            else:
                seed_predictions = np.stack(
                    [_predictor(spec, seed, f"{system}:{phase}", reference) for seed in seeds]
                )
                mean, std = seed_predictions.mean(axis=0), seed_predictions.std(axis=0)
            predictions[phase] = (temperatures, mean, seed_predictions)
            test_indexes = sorted(test_by_system_phase.get((system, phase), []))
            if test_indexes:
                test_temperatures = np.array([cache[(system, phase)][index].T_K for index in test_indexes])
                test_reference = np.array([cache[(system, phase)][index].G_eV_per_atom for index in test_indexes])
                prediction_map = dict(zip(temperatures.tolist(), mean.tolist()))
                std_map = dict(zip(temperatures.tolist(), std.tolist()))
                test_mean = np.array([prediction_map[value] for value in test_temperatures])
                test_std = np.array([std_map[value] for value in test_temperatures])
                point_errors.extend(np.abs(test_mean - test_reference).tolist())
                coverage_values.extend(np.abs(test_mean - test_reference) <= 2.0 * test_std + 1.0e-14)
        for left_index, left in enumerate(phases):
            for right in phases[left_index + 1 :]:
                if (system, left) not in test_by_system_phase and (system, right) not in test_by_system_phase:
                    continue
                left_t, left_g, left_seed_g = predictions[left]
                right_t, right_g, right_seed_g = predictions[right]
                shared = sorted(set(left_t.tolist()) & set(right_t.tolist()))
                if not shared:
                    continue
                left_map = dict(zip(left_t.tolist(), left_g.tolist()))
                right_map = dict(zip(right_t.tolist(), right_g.tolist()))
                temperatures = np.array(shared)
                observed = np.array([left_map[t] - right_map[t] for t in shared])
                ref_left = {(cache[(system, left)][i].T_K): cache[(system, left)][i].G_eV_per_atom for i in by_system_phase[(system, left)]}
                ref_right = {(cache[(system, right)][i].T_K): cache[(system, right)][i].G_eV_per_atom for i in by_system_phase[(system, right)]}
                reference = np.array([ref_left[t] - ref_right[t] for t in shared])
                roots_reference = _root(reference, temperatures)
                # Some source tables are quantized and contain a flat zero
                # plateau.  A smooth gauge perturbation can split that same
                # plateau into adjacent numerical roots; merge only these
                # near-neighbours, leaving all other crossings intact.
                roots_predicted = _root(
                    observed,
                    temperatures,
                    merge_within_K=3.0 if np.any(reference == 0) else 0.0,
                )
                seed_roots = []
                for seed_index in range(len(seeds)):
                    seed_observed = np.array(
                        [
                            dict(zip(left_t.tolist(), left_seed_g[seed_index].tolist()))[t]
                            - dict(zip(right_t.tolist(), right_seed_g[seed_index].tolist()))[t]
                            for t in shared
                        ]
                    )
                    seed_roots.append(
                        _root(
                            seed_observed,
                            temperatures,
                            merge_within_K=3.0 if np.any(reference == 0) else 0.0,
                        )
                    )
                delta_mae = float(np.mean(np.abs(observed - reference)))
                slope_by_root: list[float] = []
                crossing_mae: list[float] = []
                for root in roots_reference:
                    # The design calls for the local ΔG error at each
                    # reference crossing, not a whole-curve average.  Linear
                    # interpolation keeps this well-defined on nonuniform
                    # grids and makes the reported slope conversion auditable.
                    crossing_error = float(
                        np.interp(root, temperatures, observed - reference)
                    )
                    crossing_mae.append(abs(crossing_error))
                    nearest = int(np.argmin(np.abs(temperatures - root)))
                    if reference[nearest] == 0:
                        lower, upper = nearest, nearest
                        while lower > 0 and reference[lower] == 0:
                            lower -= 1
                        while upper < len(reference) - 1 and reference[upper] == 0:
                            upper += 1
                        slope = (reference[upper] - reference[lower]) / (temperatures[upper] - temperatures[lower])
                    elif nearest == 0:
                        slope = (reference[1] - reference[0]) / (temperatures[1] - temperatures[0])
                    elif nearest == len(temperatures) - 1:
                        slope = (reference[-1] - reference[-2]) / (temperatures[-1] - temperatures[-2])
                    else:
                        slope = (reference[nearest + 1] - reference[nearest - 1]) / (
                            temperatures[nearest + 1] - temperatures[nearest - 1]
                        )
                    slope_by_root.append(float(abs(slope)))
                pair_metrics[f"{system}:{left}_minus_{right}"] = {
                    "pair_support": "test_and_train_partner" if train else "test_only",
                    "delta_G_MAE_eV_per_atom": delta_mae,
                    "delta_G_RMSE_eV_per_atom": float(np.sqrt(np.mean((observed - reference) ** 2))),
                    "delta_G_MAE_at_crossing_eV_per_atom": crossing_mae,
                    "sign_accuracy": float(np.mean(np.sign(observed) == np.sign(reference))),
                    "reference_Tc_K": roots_reference,
                    "predicted_Tc_K": roots_predicted,
                    "predicted_Tc_by_seed_K": seed_roots,
                    "Tc_scatter_K": [
                        float(np.nanstd([roots[kk] for roots in seed_roots if len(roots) > kk]))
                        if any(len(roots) > kk for roots in seed_roots)
                        else math.nan
                        for kk in range(len(roots_reference))
                    ],
                    "Tc_error_K": [pred - ref for pred, ref in zip(roots_predicted, roots_reference)],
                    "crossing_slope_eV_per_atom_per_K": slope_by_root,
                    "Tc_err_from_dG_K": [
                        error / slope if slope else math.nan
                        for error, slope in zip(crossing_mae, slope_by_root)
                    ],
                    "false_crossings": max(0, len(roots_predicted) - len(roots_reference)),
                    "missed_crossings": max(0, len(roots_reference) - len(roots_predicted)),
                }
    return {
        "predictor": spec,
        "n_test_frames": len(test),
        "seeds": seeds,
        "G_MAE_eV_per_atom": float(np.mean(point_errors)) if point_errors else math.nan,
        "coverage_2sigma": float(np.mean(coverage_values)) if coverage_values else math.nan,
        "pairs": pair_metrics,
    }


def _number(value: object) -> str:
    if isinstance(value, (float, int)):
        return f"{value:.6g}"
    return "—"


def _summary_markdown(metrics: dict[str, object], title: str) -> str:
    """Render scalar and pairwise metrics as a portable Markdown table."""

    lines = [f"# {title}", "", f"Frames: {metrics['n_test_frames']}", "", f"G MAE: {_number(metrics['G_MAE_eV_per_atom'])} eV/atom", ""]
    folds = metrics.get("folds")
    metric_views: list[tuple[str, dict[str, object]]]
    if isinstance(folds, dict):
        metric_views = [(str(name), value) for name, value in folds.items() if isinstance(value, dict)]
    else:
        metric_views = [("all test frames", metrics)]
    for fold_name, view in metric_views:
        pairs = view.get("pairs", {})
        if not isinstance(pairs, dict) or not pairs:
            lines.extend([f"## {fold_name}", "", "No comparable phase pair is defined for this fold.", ""])
            continue
        lines.extend(
            [
                f"## {fold_name}",
                "",
                "| Pair | support | ΔG MAE (eV/atom) | ΔG RMSE (eV/atom) | ΔG MAE@Tc (eV/atom) | sign accuracy | reference Tc (K) | predicted Tc (K) | Tc error (K) | Tc from ΔG (K) | Tc scatter (K) | false / missed |",
                "|---|---|---:|---:|---|---:|---|---|---|---|---|---:|",
            ]
        )
        for pair_name, pair in sorted(pairs.items()):
            if not isinstance(pair, dict):
                continue
            crossing = lambda key: ", ".join(_number(value) for value in pair.get(key, []) if isinstance(value, (float, int))) or "—"
            lines.append(
                "| "
                + " | ".join(
                    (
                        str(pair_name),
                        str(pair.get("pair_support", "test_only")),
                        _number(pair.get("delta_G_MAE_eV_per_atom")),
                        _number(pair.get("delta_G_RMSE_eV_per_atom")),
                        crossing("delta_G_MAE_at_crossing_eV_per_atom"),
                        _number(pair.get("sign_accuracy")),
                        crossing("reference_Tc_K"),
                        crossing("predicted_Tc_K"),
                        crossing("Tc_error_K"),
                        crossing("Tc_err_from_dG_K"),
                        crossing("Tc_scatter_K"),
                        f"{pair.get('false_crossings', 0)} / {pair.get('missed_crossings', 0)}",
                    )
                )
                + " |"
            )
        lines.append("")
    return "\n".join(lines)


def run(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    split_path, data_root = Path(args.split).resolve(), Path(args.data_root).resolve()
    split = json.loads(split_path.read_text(encoding="utf-8"))
    seeds = [int(item) for item in args.seeds.split(",") if item]
    metrics = evaluate(args.predictor, split, data_root, seeds)
    out = Path(args.output_root).resolve() / _output_slug(args.predictor) / split_path.stem
    out.mkdir(parents=True, exist_ok=True)
    (out / "metrics.json").write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (out / "summary.md").write_text(
        _summary_markdown(metrics, f"{args.predictor} on {split_path.stem}"), encoding="utf-8"
    )
    print(out / "metrics.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(run())
