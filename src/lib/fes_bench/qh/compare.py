"""Calibrate and inspect a QH baseline only on a configured training range."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from fes_bench.config import ConfigError, load_mapping


def _resolve(config_path: Path, value: str) -> Path:
    value_path = Path(value).expanduser()
    return value_path if value_path.is_absolute() else (config_path.parent / value_path).resolve()


def _rows(path: Path) -> dict[float, dict[str, float]]:
    records: dict[float, dict[str, float]] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            numeric: dict[str, float] = {}
            for key, value in row.items():
                if key == "T_K":
                    continue
                try:
                    numeric[key] = float(value)
                except (TypeError, ValueError):
                    continue
            records[float(row["T_K"])] = numeric
    return records

def run(config_path: str | Path) -> int:
    try:
        import numpy as np
    except ImportError as exc:
        raise ConfigError("numpy is required for QH comparison") from exc
    try:
        import matplotlib.pyplot as plt
        from matplotlib import gridspec
        import seaborn as sns
        have_plot = True
    except ImportError as exc:
        have_plot = False
    config_file = Path(config_path).expanduser().resolve()
    config = load_mapping(config_file)
    for key in ("data_root", "result_root", "system"):
        if not isinstance(config.get(key), str):
            raise ConfigError(f"{key} is required")
    phases = config.get("phases")
    train_max = float(config.get("train_T_max_K"))
    if not isinstance(phases, list) or not all(isinstance(phase, str) for phase in phases):
        raise ConfigError("phases must be a string list")
    system_dir = _resolve(config_file, config["data_root"]) / config["system"]
    result_dir = _resolve(config_file, config["result_root"]) / config["system"]
    merged: dict[str, list[dict[str, float]]] = {}
    train_offsets: list[float] = []
    for phase in phases:
        reference = _rows(system_dir / phase / "reference_G.csv")
        qh = _rows(result_dir / f"{phase}_fqh.csv")
        common = sorted(set(reference) & set(qh))
        if not common:
            raise ConfigError(f"{phase} has no common QH/reference temperatures")
        records = [
            {"T_K": temperature, "reference": reference[temperature]["G_eV_per_atom"], "qh": qh[temperature]["F_QH_eV_per_atom"]}
            for temperature in common
        ]
        selected = [record["reference"] - record["qh"] for record in records if record["T_K"] <= train_max]
        if not selected:
            raise ConfigError(f"{phase} has no train-range QH data")
        train_offsets.extend(selected)
        merged[phase] = records
    calibration = float(np.mean(train_offsets))
    summary = result_dir / "qh_summary.json"
    summary_payload = json.loads(summary.read_text(encoding="utf-8")) if summary.exists() else {}
    metrics: dict[str, object] = {"system": config["system"], "train_T_max_K": train_max, "calibration_eV_per_atom": calibration, "checkpoint": summary_payload.get("checkpoint", ""), "checkpoint_sha256": summary_payload.get("checkpoint_sha256", ""), "head": summary_payload.get("head", ""), "phases": {}}
    all_test_errors: list[float] = []
    for phase, records in merged.items():
        for record in records:
            record["prediction"] = record["qh"] + calibration
            record["residual_meV_per_atom"] = (record["prediction"] - record["reference"]) * 1000.0
        test = [abs(record["residual_meV_per_atom"]) for record in records if record["T_K"] > train_max]
        all_test_errors.extend(test)
        values = [record["prediction"] for record in records]
        slopes = np.diff(values)
        metrics["phases"][phase] = {
            "n_points": len(records), "test_mae_meV_per_atom": float(np.mean(test)),
            "test_max_abs_error_meV_per_atom": float(np.max(test)),
            "nonmonotonic_positive_steps": int(np.count_nonzero(slopes > 1e-10)),
        }
    metrics["test_mae_meV_per_atom"] = float(np.mean(all_test_errors))
    metrics["test_max_abs_error_meV_per_atom"] = float(np.max(all_test_errors))
    (result_dir / "qh_comparison.json").write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")

    if not have_plot:
        print(f"QH comparison complete (metrics only): test_MAE_meV_per_atom={metrics['test_mae_meV_per_atom']:.6g}")
        return 0
    sns.set_theme(font_scale=1.0, style="whitegrid", font="DejaVu Sans")
    palette = sns.color_palette("colorblind", len(phases))
    n_phases = len(phases)
    figure = plt.figure(figsize=(5.2 * n_phases, 7), dpi=150)
    grid = gridspec.GridSpec(2, n_phases)
    grid.update(wspace=0.08, hspace=0.30, left=0.08, right=0.99, top=0.91, bottom=0.10)
    axes = [plt.subplot(grid[row, column]) for row in range(2) for column in range(n_phases)]
    figure.suptitle(f"{config['system']} QH baseline calibrated only below {train_max:g} K", fontsize=14, color="dimgrey")
    for column, (phase, color) in enumerate(zip(phases, palette)):
        records = merged[phase]
        temperature = np.array([record["T_K"] for record in records])
        reference = np.array([record["reference"] for record in records])
        prediction = np.array([record["prediction"] for record in records])
        residual = np.array([record["residual_meV_per_atom"] for record in records])
        curve_axis, residual_axis = axes[column], axes[n_phases + column]
        curve_axis.plot(temperature, reference, color="dimgrey", linewidth=2, label="reference G")
        curve_axis.plot(temperature, prediction, color=color, linewidth=2, linestyle="--", label="E + F_QH + c")
        curve_axis.axvline(train_max, color="lightgrey", linewidth=1.2, linestyle=":")
        curve_axis.set_title(r"$\bf{(" + chr(ord("a") + column) + r")}\ $" + phase, loc="left", fontsize=11, pad=7)
        curve_axis.tick_params(axis="both", which="both", length=0, labelcolor="dimgrey")
        curve_axis.grid(False); curve_axis.patch.set_edgecolor("lightgrey"); curve_axis.patch.set_linewidth(0.8)
        residual_axis.axhline(0.0, color="lightgrey", linewidth=1.2)
        residual_axis.plot(temperature, residual, color=color, linewidth=1.8)
        residual_axis.axvline(train_max, color="lightgrey", linewidth=1.2, linestyle=":")
        residual_axis.set_title(r"$\bf{(" + chr(ord("c") + column) + r")}\ $" + phase + " residual", loc="left", fontsize=11, pad=7)
        residual_axis.tick_params(axis="both", which="both", length=0, labelcolor="dimgrey")
        residual_axis.grid(False); residual_axis.patch.set_edgecolor("lightgrey"); residual_axis.patch.set_linewidth(0.8)
        if column == 0:
            curve_axis.set_ylabel("Free energy (eV/atom)", color="dimgrey")
            residual_axis.set_ylabel("Residual (meV/atom)", color="dimgrey")
        else:
            curve_axis.tick_params(labelleft=False); residual_axis.tick_params(labelleft=False)
        residual_axis.set_xlabel("Temperature (K)", color="dimgrey")
    axes[0].legend(frameon=True, facecolor="white", framealpha=0.8, edgecolor="lightgrey", labelcolor="dimgrey", loc="best", fontsize=9)
    figure.text(0.99, 0.01, f"Test-range MAE: {metrics['test_mae_meV_per_atom']:.1f} meV/atom; dotted line = training cutoff", ha="right", va="bottom", fontsize=9, color="dimgrey", style="italic")
    sns.despine(left=True, bottom=True)
    figure.savefig(result_dir / "fqh_vs_reference.png", dpi=150, bbox_inches="tight")
    figure.savefig(result_dir / "fqh_vs_reference.pdf", dpi=150, bbox_inches="tight")
    print(f"QH comparison complete: test_MAE_meV_per_atom={metrics['test_mae_meV_per_atom']:.6g}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    args = parser.parse_args(argv)
    try:
        return run(args.config)
    except (ConfigError, OSError, ValueError, KeyError) as exc:
        print(f"QH comparison failed: {exc}")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
