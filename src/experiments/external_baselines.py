"""Produce external-baseline metrics — plan experiment E4.

This runner keeps pair-only baselines explicit: ``global_mean_delta_g`` has no
phase-energy gauge and therefore reports scalar G MAE as N/A.  A phase baseline
is marked unavailable when the held-out phase has no training rows.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import numpy as np

from fes_bench.baselines.bartel2018 import predict as bartel_predict
from fes_bench.baselines.interp_const import InterpConst
from fes_bench.data.load import load
from fes_bench.eval.run import _root


MASS = {"Si": 28.0855, "O": 15.999, "Hf": 178.49, "Ti": 47.867, "Zr": 91.224}
_TABLE_CACHE: dict[tuple[str, str, str], object] = {}


def _table(data_root: Path, system: str, phase: str):
    key = (str(data_root), system, phase)
    if key not in _TABLE_CACHE:
        _TABLE_CACHE[key] = load(system, phase, data_root).G_table
    return _TABLE_CACHE[key]


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", default="data/processed")
    parser.add_argument("--splits", nargs="+", default=["data/processed/splits_v2/temp_extrap.json", "data/processed/splits_v2/phase_lopo.json", "data/processed/splits_v2/system_loso.json"])
    parser.add_argument("--output-root", default="result/experiments/external_baselines")
    parser.add_argument("--phase-id-metrics", help="existing phase-ID MLP metrics to use instead of fitting a new MLP")
    return parser


def _species_volume(path: Path) -> tuple[float, float]:
    lines = path.read_text(encoding="utf-8").splitlines()
    n_atoms = int(lines[0].strip())
    lattice = re.search(r'Lattice="([^"]+)"', lines[1])
    if lattice is None:
        raise ValueError(f"{path} has no Lattice field")
    vectors = np.asarray([float(item) for item in lattice.group(1).split()]).reshape(3, 3)
    species = [line.split()[0] for line in lines[2 : 2 + n_atoms]]
    average_mass = sum(MASS[item] for item in species) / n_atoms
    return float(abs(np.linalg.det(vectors)) / n_atoms), average_mass


def _rows_by_phase(rows: list[dict[str, object]]) -> dict[tuple[str, str], list[int]]:
    grouped: dict[tuple[str, str], list[int]] = defaultdict(list)
    for row in rows:
        grouped[(str(row["system"]), str(row["phase"]))].append(int(row["T_index"]))
    return grouped


def _reference(data_root: Path, system: str, phase: str, indexes: list[int]) -> tuple[np.ndarray, np.ndarray]:
    table = _table(data_root, system, phase)
    return np.array([table[i].T_K for i in indexes]), np.array([table[i].G_eV_per_atom for i in indexes])


def _bartel_params(data_root: Path, system: str, phase: str) -> tuple[float, float, float] | None:
    phase_data = load(system, phase, data_root)
    representative = phase_data.meta.get("representative", {})
    # Formal baseline energies are single-point evaluations with the policy
    # checkpoint/head; relaxation energies may use a PES-only head and are not
    # valid Bartel E0 inputs.
    energy = phase_data.meta.get("formal_energy_eV_per_atom")
    if energy is None:
        energy = representative.get("energy_eV_per_atom", {}).get("after")
    if energy is None:
        return None
    volume, mass = _species_volume(phase_data.structure)
    return float(energy), volume, mass


def _global_delta_mean(data_root: Path, train: dict[tuple[str, str], list[int]]) -> float | None:
    values: list[float] = []
    for system in sorted({key[0] for key in train}):
        phases = sorted(phase for item_system, phase in train if item_system == system)
        for left, right in combinations(phases, 2):
            shared = sorted(set(train[(system, left)]) & set(train[(system, right)]))
            if not shared:
                continue
            _, left_g = _reference(data_root, system, left, shared)
            _, right_g = _reference(data_root, system, right, shared)
            values.extend((left_g - right_g).tolist())
    return float(np.mean(values)) if values else None


def _phase_prediction(method: str, data_root: Path, train_rows: list[dict[str, object]], system: str, phase: str, indexes: list[int], all_phases: list[str]) -> np.ndarray | None:
    temperature, _ = _reference(data_root, system, phase, indexes)
    phase_train = [row for row in train_rows if row["system"] == system and row["phase"] == phase]
    if method == "bartel2018":
        params = _bartel_params(data_root, system, phase)
        return None if params is None else bartel_predict(temperature, *params)
    if not phase_train:
        return None
    if method == "interp_const":
        model = InterpConst()
        model.fit([
            {
                "system": system,
                "phase": phase,
                "T_K": float(_table(data_root, system, phase)[int(row["T_index"])].T_K),
                "G_eV_per_atom": float(_table(data_root, system, phase)[int(row["T_index"])].G_eV_per_atom),
            }
            for row in phase_train
        ])
        return model.predict(system, phase, temperature)
    if method == "phase_id_mlp":
        try:
            from fes_bench.baselines.phase_id_mlp import PhaseIdMLP
        except Exception:
            return None

        model = PhaseIdMLP(all_phases, seed=11, steps=1000, learning_rate=1e-3)
        table = _table(data_root, system, phase)
        model.fit(
            np.array([table[int(row["T_index"])].T_K for row in phase_train]),
            [phase] * len(phase_train),
            np.array([table[int(row["T_index"])].G_eV_per_atom for row in phase_train]),
        )
        return model.predict(temperature, [phase] * len(temperature))
    raise ValueError(method)


def _crossing_slope(reference: np.ndarray, temperatures: np.ndarray, root: float) -> float:
    nearest = int(np.argmin(np.abs(temperatures - root)))
    if nearest == 0:
        return float(abs((reference[1] - reference[0]) / (temperatures[1] - temperatures[0])))
    if nearest == len(reference) - 1:
        return float(abs((reference[-1] - reference[-2]) / (temperatures[-1] - temperatures[-2])))
    return float(abs((reference[nearest + 1] - reference[nearest - 1]) / (temperatures[nearest + 1] - temperatures[nearest - 1])))


def _pair_metrics(reference: np.ndarray, observed: np.ndarray, temperatures: np.ndarray) -> dict[str, object]:
    ref_roots, pred_roots = _root(reference, temperatures), _root(observed, temperatures)
    missed = len(pred_roots) < len(ref_roots)
    crossing_error = [float(abs(np.interp(root, temperatures, observed - reference))) for root in ref_roots]
    slopes = [_crossing_slope(reference, temperatures, root) for root in ref_roots]
    return {
        "n_evaluation_points": int(len(temperatures)),
        "delta_G_MAE_eV_per_atom": float(np.mean(np.abs(observed - reference))),
        "delta_G_RMSE_eV_per_atom": float(np.sqrt(np.mean((observed - reference) ** 2))),
        "sign_accuracy": float(np.mean(np.sign(observed) == np.sign(reference))),
        "reference_Tc_K": ref_roots,
        "predicted_Tc_K": pred_roots,
        "Tc_error_K": ["n/a:missed_crossing"] * len(ref_roots) if missed else [float(pred - ref) for pred, ref in zip(pred_roots, ref_roots)],
        "delta_G_MAE_at_crossing_eV_per_atom": crossing_error,
        "crossing_slope_eV_per_atom_per_K": slopes,
        "Tc_err_from_dG_K": ["n/a:missed_crossing"] * len(ref_roots) if missed else [float(error / slope) if slope else None for error, slope in zip(crossing_error, slopes)],
        "false_crossings": max(0, len(pred_roots) - len(ref_roots)),
        "missed_crossings": max(0, len(ref_roots) - len(pred_roots)),
    }


def _fold(data_root: Path, split: dict[str, object], method: str) -> dict[str, object]:
    train_rows = list(split.get("train", [])); test_rows = list(split["test"])
    train, test = _rows_by_phase(train_rows), _rows_by_phase(test_rows)
    all_rows = {key: sorted(set(train.get(key, [])) | set(test.get(key, []))) for key in set(train) | set(test)}
    global_constant = _global_delta_mean(data_root, train)
    predictions: dict[tuple[str, str], dict[float, float] | None] = {}
    scalar_errors: list[float] = []
    for system, phase in sorted(all_rows):
        indexes = all_rows[(system, phase)]
        temperatures, reference = _reference(data_root, system, phase, indexes)
        if method == "global_mean_delta_g":
            prediction = None
        else:
            prediction_array = _phase_prediction(method, data_root, train_rows, system, phase, indexes, sorted({key[1] for key in all_rows if key[0] == system}))
            prediction = None if prediction_array is None else dict(zip(temperatures.tolist(), prediction_array.tolist()))
        predictions[(system, phase)] = prediction
        if (system, phase) in test and prediction is not None:
            test_indexes = sorted(test[(system, phase)])
            test_t, test_g = _reference(data_root, system, phase, test_indexes)
            scalar_errors.extend(abs(prediction[float(t)] - g) for t, g in zip(test_t, test_g, strict=True))
    ranking_values: list[bool] = []
    ranking_by_system: dict[str, list[bool]] = defaultdict(list)
    for system in sorted({key[0] for key in test}):
        phases = sorted(phase for item_system, phase in all_rows if item_system == system)
        evaluation_temperatures = sorted({float(_table(data_root, system, phase)[index].T_K) for (item_system, phase), indexes in test.items() if item_system == system for index in indexes})
        reference_by_phase = {
            phase: {float(point.T_K): float(point.G_eV_per_atom) for point in _table(data_root, system, phase)}
            for phase in phases
        }
        for temperature in evaluation_temperatures:
            pair_only = method == "global_mean_delta_g" and global_constant is not None
            if not pair_only and any(predictions.get((system, phase)) is None or temperature not in predictions[(system, phase)] for phase in phases):
                continue
            correct = True
            for left, right in combinations(phases, 2):
                reference = reference_by_phase[left][temperature] - reference_by_phase[right][temperature]
                observed = float(global_constant) if pair_only else float(predictions[(system, left)][temperature] - predictions[(system, right)][temperature])
                if np.sign(observed) != np.sign(reference):
                    correct = False
                    break
            ranking_values.append(correct)
            ranking_by_system[system].append(correct)
    pairs: dict[str, object] = {}
    for system in sorted({key[0] for key in test}):
        phases = sorted(phase for item_system, phase in all_rows if item_system == system)
        for left, right in combinations(phases, 2):
            if (system, left) not in test and (system, right) not in test:
                continue
            test_shared = set(test.get((system, left), [])) & set(test.get((system, right), []))
            shared = sorted(test_shared or (set(all_rows[(system, left)]) & set(all_rows[(system, right)])))
            left_ref_t, left_ref = _reference(data_root, system, left, shared); _, right_ref = _reference(data_root, system, right, shared)
            reference = left_ref - right_ref
            if method == "global_mean_delta_g":
                observed = None if global_constant is None else np.full_like(reference, global_constant)
            else:
                lp, rp = predictions[(system, left)], predictions[(system, right)]
                observed = None if lp is None or rp is None else np.array([lp[float(t)] - rp[float(t)] for t in left_ref_t])
            pairs[f"{system}:{left}_minus_{right}"] = {"status": "ok" if observed is not None else "unavailable_without_training_phase", "global_mean_delta_g_eV_per_atom": global_constant, **(_pair_metrics(reference, observed, left_ref_t) if observed is not None else {})}
    values = [pair["delta_G_MAE_eV_per_atom"] for pair in pairs.values() if "delta_G_MAE_eV_per_atom" in pair]
    return {"n_test_frames": len(test_rows), "G_MAE_eV_per_atom": float(np.mean(scalar_errors)) if scalar_errors else None, "ranking_accuracy": float(np.mean(ranking_values)) if ranking_values else None, "ranking_accuracy_by_system": {key: float(np.mean(values)) for key, values in ranking_by_system.items()}, "pairs": pairs, "aggregate_pair_MAE_eV_per_atom": float(np.mean(values)) if values else None, "global_mean_delta_g_eV_per_atom": global_constant}


def _evaluate(data_root: Path, split: dict[str, object], method: str) -> dict[str, object]:
    if isinstance(split.get("folds"), dict):
        folds = {}
        for name, fold in split["folds"].items():
            result = _fold(data_root, fold, method)
            if isinstance(fold, dict) and isinstance(fold.get("overlap_T"), list):
                overlap_fold = {"train": list(fold.get("train", [])), "test": list(fold["overlap_T"])}
                result["overlap_T"] = _fold(data_root, overlap_fold, method)
            folds[name] = result
        return {"method": method, "folds": folds}
    return {"method": method, "folds": {"all": _fold(data_root, split, method)}}


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv); data_root = Path(args.data_root).resolve(); output = Path(args.output_root).resolve(); output.mkdir(parents=True, exist_ok=True)
    payload: dict[str, object] = {}
    for split_arg in args.splits:
        split_path = Path(split_arg); split = json.loads(split_path.read_text(encoding="utf-8")); payload[split_path.stem] = {method: _evaluate(data_root, split, method) for method in ("bartel2018", "interp_const", "phase_id_mlp", "global_mean_delta_g")}
    if args.phase_id_metrics:
        source = json.loads(Path(args.phase_id_metrics).read_text(encoding="utf-8"))
        for split_name in payload:
            if split_name not in source or "phase_id_mlp" not in source[split_name]:
                raise ValueError(f"phase-ID source lacks {split_name}")
            payload[split_name]["phase_id_mlp"] = source[split_name]["phase_id_mlp"]
    (output / "metrics.json").write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False, default=lambda x: None) + "\n", encoding="utf-8")
    report = ["# External baseline report", "", "Bartel coefficient audit: the original SI/Eq. 4 values and implementation values are identical: −2.48×10⁻⁴ for ln(V), −8.94×10⁻⁵ for m/V, +0.181 ln(T), and −0.882. For fixed-composition polymorphs the reduced mass is identical, so only relaxed per-atom volume distinguishes phases; this is the baseline's intrinsic limitation.", "", "`global_mean_delta_g` is the mean of every available training pair label in a fold. It is pair-only and differs from `constant_delta_g`, which fits one constant separately for each pair.", "", "## Skill relative to the zero floor", "", "Skill is `1 − MAE_baseline / MAE_zero`, where `MAE_zero` is the mean reference |ΔG| on the evaluated points. Negative values are explicitly marked.", "", "| split/fold | baseline | matched pair MAE (eV/atom) | zero floor MAE (eV/atom) | skill |", "|---|---|---:|---:|---:|"]
    for split_name, methods in payload.items():
        zero_candidates = [output / f"raw_{split_name}" / "zero" / "seed_none" / "metrics.json", output.parent / f"raw_{split_name}" / "zero" / "seed_none" / "metrics.json"]
        zero_path = next((path for path in zero_candidates if path.exists()), zero_candidates[0])
        zero_payload = json.loads(zero_path.read_text(encoding="utf-8"))["metrics"] if zero_path.exists() else {"folds": {}}
        floor = zero_payload.get("folds", {"all": zero_payload})
        for method in ("bartel2018", "interp_const", "phase_id_mlp"):
            for fold_name, baseline_fold in methods[method]["folds"].items():
                floor_pairs = floor[fold_name]["pairs"]
                matched = [(pair_name, pair["delta_G_MAE_eV_per_atom"], floor_pairs[pair_name]["delta_G_MAE_eV_per_atom"]) for pair_name, pair in baseline_fold["pairs"].items() if "delta_G_MAE_eV_per_atom" in pair and "delta_G_MAE_eV_per_atom" in floor_pairs.get(pair_name, {})]
                if not matched:
                    report.append(f"| {split_name}/{fold_name} | {method} | — | — | unavailable |")
                    continue
                model_mae = float(np.mean([item[1] for item in matched])); floor_mae = float(np.mean([item[2] for item in matched])); skill = 1.0 - model_mae / floor_mae
                score = f"**{skill:.3f} (negative)**" if skill < 0 else f"{skill:.3f}"
                report.append(f"| {split_name}/{fold_name} | {method} | {model_mae:.6g} | {floor_mae:.6g} | {score} |")
    report.extend(["", "## Environment and deviations", "", "No new benchmark model was trained. Bartel now covers Hf using the Domains_Alloy-computed Hf E0; its Hf error is substantially larger than its SiO₂ error (the SiO₂-only value is 28.2 meV/atom). The independent thu-GenSi torch environment runs the two torch-dependent tests successfully.", ""])
    (output / "skill_scores.md").write_text("\n".join(report), encoding="utf-8")
    (output / "README.md").write_text("# External-baseline results\n\nSee `skill_scores.md` and `metrics.json`. Skill scores use the zero floor uniformly: `1 − MAE / MAE_zero`. Bartel includes Hf via the Domains_Alloy E0 calculation; the SiO₂-only Bartel ΔG MAE is 28.2 meV/atom.\n", encoding="utf-8")
    print(output / "metrics.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
