"""Re-evaluate historic four-phase SiO2 FES checkpoints on a named three-phase subset.

The script is intentionally an adapter, not a trainer.  It evaluates frozen
TorchScript heads on their original static structures, writes temperature-keyed
curves, and then invokes ``fes_bench.eval.run`` with its ``precomputed:``
predictor mode.  Thus the benchmark evaluator, rather than this adapter,
owns root finding, slope conversion, and false/missed-crossing accounting.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np


CANONICAL_PHASES = ("quartz_beta", "cristobalite_beta", "tridymite_p63mmc")
LEGACY_PHASES = CANONICAL_PHASES + ("tridymite_c2221",)


def read_csv(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    return (
        np.array([float(row["T_K"]) for row in rows]),
        np.array([float(row["G_eV_per_atom"]) for row in rows]),
    )


def raw_phase(legacy_root: Path, phase: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    phase_root = legacy_root / phase
    atype = np.loadtxt(phase_root / "type.raw", dtype=np.int64)
    parts = []
    for directory in sorted(phase_root.glob("set.*")):
        parts.append(
            (
                np.load(directory / "fparam.npy"),
                np.load(directory / "coord.npy"),
                np.load(directory / "box.npy"),
                np.load(directory / "free_energy.npy"),
            )
        )
    fparam, coord, box, free_energy = (np.concatenate([part[i] for part in parts]) for i in range(4))
    order = np.argsort(fparam[:, 0])
    return fparam[order], coord[order], box[order], free_energy[order], atype, len(atype)


def reference_rows(data_root: Path, phase: str) -> tuple[np.ndarray, np.ndarray]:
    return read_csv(data_root / "sio2" / phase / "reference_G.csv")


def evaluate_model(model_path: Path, legacy_root: Path, device: str, batch_size: int) -> dict[str, dict[str, list[float]]]:
    import torch

    torch.serialization.add_safe_globals([slice])
    from deepmd.infer import DeepPot

    wrapper = DeepPot(str(model_path))
    model = wrapper.deep_eval.dp.model["Default"].eval()
    output: dict[str, dict[str, list[float]]] = {}
    for phase in CANONICAL_PHASES:
        fparam, coord, box, _labels, atype, natoms = raw_phase(legacy_root, phase)
        predictions = []
        for start in range(0, len(fparam), batch_size):
            stop = min(start + batch_size, len(fparam))
            coord_t = torch.tensor(coord[start:stop].reshape(stop - start, natoms, 3), dtype=torch.float64, device=device)
            box_t = torch.tensor(box[start:stop].reshape(stop - start, 3, 3), dtype=torch.float64, device=device)
            atype_t = torch.tensor(np.tile(atype, (stop - start, 1)), dtype=torch.long, device=device)
            fparam_t = torch.tensor(fparam[start:stop], dtype=torch.float64, device=device)
            with torch.no_grad():
                values = model(coord_t, atype_t, box_t, fparam_t)["free_energy"]
            predictions.extend((values.detach().cpu().numpy().reshape(-1) / natoms).tolist())
        output[f"sio2:{phase}"] = {"T_K": fparam[:, 0].astype(float).tolist(), "G_eV_per_atom": predictions}
    return output


def checkpoint_type_map(model_path: Path) -> list[str]:
    """Read the element mapping embedded in the frozen TorchScript model."""

    import torch

    torch.serialization.add_safe_globals([slice])
    from deepmd.infer import DeepPot

    return list(DeepPot(str(model_path)).get_type_map())


def rows_for_temperatures(data_root: Path, temperatures: set[float]) -> list[dict[str, object]]:
    rows = []
    for phase in CANONICAL_PHASES:
        source_t, _ = reference_rows(data_root, phase)
        for index, temperature in enumerate(source_t):
            if float(temperature) in temperatures:
                rows.append({"system": "sio2", "phase": phase, "T_index": index})
    return rows


def calibration(predictions: dict[str, dict[str, list[float]]], data_root: Path, temperatures: set[float]) -> float:
    residuals: list[float] = []
    for phase in CANONICAL_PHASES:
        reference_t, reference_g = reference_rows(data_root, phase)
        predicted = predictions[f"sio2:{phase}"]
        by_t = dict(zip(predicted["T_K"], predicted["G_eV_per_atom"]))
        for temperature, value in zip(reference_t, reference_g):
            if float(temperature) in temperatures:
                residuals.append(float(value) - float(by_t[float(temperature)]))
    if not residuals:
        raise ValueError("calibration window has no shared frames")
    return float(np.mean(residuals))


def shifted(predictions: dict[str, dict[str, list[float]]], constant: float) -> dict[str, dict[str, list[float]]]:
    return {
        key: {"T_K": curve["T_K"], "G_eV_per_atom": [float(value) + constant for value in curve["G_eV_per_atom"]]}
        for key, curve in predictions.items()
    }


def invoke_eval(repo: Path, curves: Path, split: Path, data_root: Path, output_root: Path, seed: int) -> dict[str, object]:
    spec = f"precomputed:{curves}"
    command = [
        sys.executable,
        "-m",
        "fes_bench.eval.run",
        "--predictor",
        spec,
        "--split",
        str(split),
        "--data-root",
        str(data_root),
        "--seeds",
        str(seed),
        "--output-root",
        str(output_root),
    ]
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(repo) + os.pathsep + environment.get("PYTHONPATH", "")
    subprocess.run(command, cwd=repo, env=environment, check=True, capture_output=True, text=True)
    slug = spec.replace(":", "_").replace("/", "_")
    metrics_path = output_root / slug / split.stem / "metrics.json"
    return json.loads(metrics_path.read_text(encoding="utf-8"))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint-root", required=True, type=Path)
    parser.add_argument("--legacy-data-root", required=True, type=Path)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument(
        "--train-max-index",
        type=int,
        default=591,
        help="inclusive canonical temp_extrap training index (current frozen split: 591)",
    )
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--batch-size", type=int, default=16)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    args.output_root.mkdir(parents=True, exist_ok=True)

    protocol = json.loads((args.legacy_data_root / "static_protocol.json").read_text(encoding="utf-8"))
    legacy_input_type_map = (args.legacy_data_root / "type_map.raw").read_text(encoding="utf-8").split()
    first_checkpoint = next((args.checkpoint_root / "polynomial").glob("seed*/model.pth"))
    type_map = checkpoint_type_map(first_checkpoint)
    if type_map != legacy_input_type_map:
        raise ValueError(f"checkpoint type_map {type_map!r} disagrees with legacy inputs {legacy_input_type_map!r}")
    map_rows = []
    for index, legacy_phase in enumerate(LEGACY_PHASES):
        fparam, _coord, _box, labels, _atype, natoms = raw_phase(args.legacy_data_root, legacy_phase)
        item: dict[str, object] = {
            "checkpoint_phase_index": index,
            "legacy_phase": legacy_phase,
            "canonical_phase": legacy_phase if legacy_phase in CANONICAL_PHASES else None,
            "n_atoms": natoms,
            "n_frames": len(fparam),
            "temperature_range_K": [float(fparam[0, 0]), float(fparam[-1, 0])],
            "type_map": type_map,
            "legacy_protocol_n_atoms": protocol["phases"][legacy_phase]["n_atoms"],
        }
        if legacy_phase in CANONICAL_PHASES:
            canonical_t, canonical_g = reference_rows(args.data_root, legacy_phase)
            raw_per_atom = labels.reshape(-1) / natoms
            raw_by_t = dict(zip(fparam[:, 0].astype(float), raw_per_atom))
            canonical_by_t = dict(zip(canonical_t.astype(float), canonical_g))
            common = sorted(set(raw_by_t) & set(canonical_by_t))
            item["reference_G_validation"] = {
                "n_common_temperatures": len(common),
                "max_abs_difference_eV_per_atom": float(max(abs(raw_by_t[t] - canonical_by_t[t]) for t in common)),
                "raw_reference_range_eV_per_atom": [float(min(raw_by_t.values())), float(max(raw_by_t.values()))],
                "canonical_reference_range_eV_per_atom": [float(min(canonical_by_t.values())), float(max(canonical_by_t.values()))],
            }
        map_rows.append(item)
    phase_map = {
        "scope": "historic four-phase checkpoints evaluated only on the named canonical three-phase subset",
        "checkpoint_type_map": type_map,
        "legacy_input_type_map": legacy_input_type_map,
        "checkpoint_phase_order_source": "training_data.systems order in adjacent config.json; model type_map contains species, not phase IDs",
        "entries": map_rows,
        "excluded_from_reporting": "tridymite_c2221",
    }
    (args.output_root / "sio2_checkpoint_phase_map.json").write_text(json.dumps(phase_map, indent=2) + "\n", encoding="utf-8")

    canonical_temperatures, _ = reference_rows(args.data_root, CANONICAL_PHASES[0])
    for phase in CANONICAL_PHASES[1:]:
        phase_temperatures, _ = reference_rows(args.data_root, phase)
        if not np.array_equal(phase_temperatures, canonical_temperatures):
            raise ValueError("canonical three-phase tables do not share a temperature grid")
    train_temperatures = {float(value) for value in canonical_temperatures[: args.train_max_index + 1]}
    held_temperatures = {float(value) for value in canonical_temperatures[args.train_max_index + 1 :]}
    all_temperatures = train_temperatures | held_temperatures
    full_split = {"name": "sio2_three_phase_full_range", "train": [], "test": rows_for_temperatures(args.data_root, all_temperatures)}
    subset_split = {
        "name": "sio2_temp_extrap_subset",
        "train": rows_for_temperatures(args.data_root, train_temperatures),
        "test": rows_for_temperatures(args.data_root, held_temperatures),
    }
    full_split_path = args.output_root / "sio2_three_phase_full_range.json"
    subset_split_path = args.output_root / "sio2_temp_extrap_subset.json"
    full_split_path.write_text(json.dumps(full_split, indent=2) + "\n", encoding="utf-8")
    subset_split_path.write_text(json.dumps(subset_split, indent=2) + "\n", encoding="utf-8")

    records: list[dict[str, object]] = []
    for basis in ("polynomial", "tlog_polynomial"):
        for seed_dir in sorted((args.checkpoint_root / basis).glob("seed*")):
            seed = int(seed_dir.name.removeprefix("seed"))
            config = json.loads((seed_dir / "config.json").read_text(encoding="utf-8"))
            raw = evaluate_model(seed_dir / "model.pth", args.legacy_data_root, args.device, args.batch_size)
            formal_c = calibration(raw, args.data_root, train_temperatures)
            held_c = calibration(raw, args.data_root, held_temperatures)
            variants = {"formal_train_window": formal_c, "heldout_window_audit": held_c}
            result: dict[str, object] = {
                "basis": basis,
                "seed": seed,
                "checkpoint": str(seed_dir / "model.pth"),
                "config_training_systems": config["training"]["training_data"]["systems"],
                "model_internal_c_system_field": config["model"]["fitting_net"].get("c_system"),
                "calibration_constants_eV_per_atom": variants,
                "metrics": {},
            }
            for label, constant in variants.items():
                curve_path = args.output_root / "curves" / f"{basis}_seed{seed}_{label}.json"
                curve_path.parent.mkdir(parents=True, exist_ok=True)
                curve_path.write_text(json.dumps({"predictions": shifted(raw, constant)}, indent=2) + "\n", encoding="utf-8")
                result["metrics"][label] = {
                    "full_range": invoke_eval(repo, curve_path, full_split_path, args.data_root, args.output_root / "eval", seed),
                    "temp_extrap": invoke_eval(repo, curve_path, subset_split_path, args.data_root, args.output_root / "eval", seed),
                }
            records.append(result)
            print(f"completed {basis} seed {seed}", flush=True)
    payload = {
        "scope": "four-phase historic SiO2 FES checkpoints; reported metrics restricted to canonical quartz/cristobalite/p63mmc subset",
        "calibration_windows_K": {
            "formal_train": [min(train_temperatures), max(train_temperatures)],
            "heldout_audit": [min(held_temperatures), max(held_temperatures)],
            "frozen_temp_extrap_train_max_index": args.train_max_index,
        },
        "checkpoint_records": records,
    }
    (args.output_root / "sio2_reanalysis.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
