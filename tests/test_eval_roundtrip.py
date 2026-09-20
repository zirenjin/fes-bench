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


def test_smooth_noise_tc_matches_slope_conversion(tmp_path: Path) -> None:
    split = _write_linear_system(tmp_path)
    metrics = evaluate("reference_noise:0.005", split, tmp_path, [11])
    pair = metrics["pairs"]["toy:alpha_minus_beta"]

    direct = abs(pair["Tc_error_K"][0])
    converted = pair["Tc_err_from_dG_K"][0]
    assert direct == pytest.approx(converted, rel=0.20)
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
