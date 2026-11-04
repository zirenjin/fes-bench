"""Regression checks for the frozen-split evaluator."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from fes_bench.eval.run import evaluate


def _write_linear_system(root: Path) -> dict[str, object]:
    system = root / "toy"
    phases = ("alpha", "beta")
    system.mkdir()
    (system / "system.json").write_text(
        json.dumps(
            {
                "phases": list(phases),
                "type_map": ["X"],
                "reference_grid": {"T_K": list(range(11)), "P_GPa": [0.0]},
                "truth_level": "unit-test",
            }
        ),
        encoding="utf-8",
    )
    meta = {
        "functional": "test", "dispersion": "none", "supercell": "1", "kpoints": "gamma",
        "convergence": "test", "size_error_eV_per_atom": 0.0, "method": "test",
        "doi": "test", "notes": "linear fixture",
    }
    test: list[dict[str, object]] = []
    for phase in phases:
        phase_dir = system / phase
        phase_dir.mkdir()
        (phase_dir / "structure.extxyz").write_text("1\nProperties=species:S:1:pos:R:3\nX 0 0 0\n", encoding="utf-8")
        (phase_dir / "meta.json").write_text(json.dumps(meta), encoding="utf-8")
        with (phase_dir / "reference_G.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(["T_K", "P_GPa", "G_eV_per_atom", "level", "source"])
            for index, temperature in enumerate(range(11)):
                # alpha - beta crosses once at 5 K with a 10 meV/K slope.
                energy = 0.0 if phase == "alpha" else 0.01 * (temperature - 5.0)
                writer.writerow([temperature, 0.0, energy, "test", "fixture"])
                test.append({"system": "toy", "phase": phase, "T_index": index})
    return {"name": "temp_extrap", "test": test}


def test_reference_passthrough_is_exact(tmp_path: Path) -> None:
    split = _write_linear_system(tmp_path)
    metrics = evaluate("reference", split, tmp_path, [11, 23, 37, 51, 67])

    assert metrics["G_MAE_eV_per_atom"] == 0.0
    pair = metrics["pairs"]["toy:alpha_minus_beta"]
    assert pair["Tc_error_K"] == [0.0]
    assert pair["delta_G_MAE_at_crossing_eV_per_atom"] == [0.0]
    assert pair["false_crossings"] == pair["missed_crossings"] == 0


def test_constant_offset_tc_matches_slope_conversion(tmp_path: Path) -> None:
    split = _write_linear_system(tmp_path)
    metrics = evaluate("constant_offset:0.005", split, tmp_path, [11])
    pair = metrics["pairs"]["toy:alpha_minus_beta"]

    direct = abs(pair["Tc_error_K"][0])
    converted = pair["Tc_err_from_dG_K"][0]
    assert direct == pytest.approx(converted, abs=1.0e-12)
    assert pair["Tc_scatter_K"] == [0.0]
    assert pair["false_crossings"] == pair["missed_crossings"] == 0


def test_root_stability_iid_noise_reports_scatter_without_gate(tmp_path: Path) -> None:
    split = _write_linear_system(tmp_path)
    metrics = evaluate("reference_iid_noise:0.005", split, tmp_path, [11, 23, 37])
    pair = metrics["pairs"]["toy:alpha_minus_beta"]

    assert len(pair["predicted_Tc_by_seed_K"]) == 3
    assert len(pair["Tc_scatter_K"]) == 1
    assert pair["Tc_scatter_K"][0] >= 0.0
    assert pair["false_crossings"] == pair["missed_crossings"] == 0


def test_phase_lopo_uses_training_phase_as_explicit_pair_partner(tmp_path: Path) -> None:
    full = _write_linear_system(tmp_path)
    test = [row for row in full["test"] if row["phase"] == "alpha"]
    train = [row for row in full["test"] if row["phase"] == "beta"]
    metrics = evaluate("reference", {"name": "phase_lopo", "train": train, "test": test}, tmp_path, [11])

    assert metrics["G_MAE_eV_per_atom"] == 0.0
    pair = metrics["pairs"]["toy:alpha_minus_beta"]
    assert pair["pair_support"] == "test_and_train_partner"
    assert pair["Tc_error_K"] == [0.0]


def test_precomputed_checkpoint_curve_uses_temperature_keys(tmp_path: Path) -> None:
    split = _write_linear_system(tmp_path)
    curves = {"predictions": {}}
    for phase in ("alpha", "beta"):
        with (tmp_path / "toy" / phase / "reference_G.csv").open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        curves["predictions"][f"toy:{phase}"] = {
            "T_K": [float(row["T_K"]) for row in rows],
            "G_eV_per_atom": [float(row["G_eV_per_atom"]) for row in rows],
        }
    path = tmp_path / "curves.json"
    path.write_text(json.dumps(curves), encoding="utf-8")

    metrics = evaluate(f"precomputed:{path}", split, tmp_path, [11])
    pair = metrics["pairs"]["toy:alpha_minus_beta"]
    assert metrics["G_MAE_eV_per_atom"] == 0.0
    assert pair["Tc_error_K"] == [0.0]


def test_multi_pressure_grid_pairs_only_matching_states(tmp_path: Path) -> None:
    system = tmp_path / "toy"
    phases = ("alpha", "beta")
    system.mkdir()
    (system / "system.json").write_text(
        json.dumps({
            "phases": list(phases), "type_map": ["X"],
            "reference_grid": {"T_K": [300, 400, 500], "P_GPa": [0, 5]},
            "truth_level": "multi-pressure-fixture",
        }), encoding="utf-8",
    )
    split: dict[str, object] = {"test": []}
    curves: dict[str, object] = {"predictions": {}}
    for phase in phases:
        phase_dir = system / phase
        phase_dir.mkdir()
        (phase_dir / "structure.extxyz").write_text(
            "1\nProperties=species:S:1:pos:R:3\nX 0 0 0\n", encoding="utf-8"
        )
        (phase_dir / "meta.json").write_text(json.dumps({
            "functional": "test", "dispersion": "none", "supercell": "1",
            "kpoints": "gamma", "convergence": "test", "size_error_eV_per_atom": 0.0,
            "method": "test", "doi": "test", "notes": "multi-pressure fixture",
        }), encoding="utf-8")
        rows: list[list[object]] = []
        pred_t: list[float] = []
        pred_p: list[float] = []
        pred_g: list[float] = []
        for pressure in (0.0, 5.0):
            for temperature in (300.0, 400.0, 500.0):
                value = 0.0 if phase == "alpha" else 0.001 * (temperature - (400.0 if pressure == 0.0 else 500.0))
                rows.append([temperature, pressure, value, "test", "fixture"])
                pred_t.append(temperature)
                pred_p.append(pressure)
                pred_g.append(value)
                split["test"].append({"system": "toy", "phase": phase, "T_index": len(rows) - 1})
        with (phase_dir / "reference_G.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(["T_K", "P_GPa", "G_eV_per_atom", "level", "source"])
            writer.writerows(rows)
        curves["predictions"][f"toy:{phase}"] = {
            "T_K": pred_t, "P_GPa": pred_p, "G_eV_per_atom": pred_g,
        }

    metrics = evaluate("reference", split, tmp_path, [11, 23])
    assert metrics["G_MAE_eV_per_atom"] == 0.0
    assert set(metrics["pairs"]) == {
        "toy:alpha_minus_beta@P=0GPa", "toy:alpha_minus_beta@P=5GPa",
    }
    assert metrics["pairs"]["toy:alpha_minus_beta@P=0GPa"]["reference_Tc_K"] == [400.0]
    assert metrics["pairs"]["toy:alpha_minus_beta@P=5GPa"]["reference_Tc_K"] == [500.0]

    curves_path = tmp_path / "multi_pressure_curves.json"
    curves_path.write_text(json.dumps(curves), encoding="utf-8")
    precomputed = evaluate(f"precomputed:{curves_path}", split, tmp_path, [11, 23])
    assert precomputed["G_MAE_eV_per_atom"] == 0.0
    assert precomputed["pairs"]["toy:alpha_minus_beta@P=0GPa"]["delta_G_MAE_eV_per_atom"] == 0.0
    assert precomputed["pairs"]["toy:alpha_minus_beta@P=5GPa"]["delta_G_MAE_eV_per_atom"] == 0.0

    temperature_only = {"predictions": {}}
    for key, curve in curves["predictions"].items():
        temperature_only["predictions"][key] = {
            "T_K": curve["T_K"], "G_eV_per_atom": curve["G_eV_per_atom"],
        }
    ambiguous_path = tmp_path / "ambiguous_temperature_only_curves.json"
    ambiguous_path.write_text(json.dumps(temperature_only), encoding="utf-8")
    with pytest.raises(ValueError, match="duplicate state"):
        evaluate(f"precomputed:{ambiguous_path}", split, tmp_path, [11, 23])
