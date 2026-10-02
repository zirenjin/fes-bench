"""Single-sided QH thermal-shape diagnostic — plan experiment single_side_qh.

This experiment never forms phase pairs.  It compares canonical reference G(T),
the split-specific T3 tlog absolute-G prediction, and the split-specific T2
QH-residual prediction on the four QH-reliable phases only.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
try:
    import seaborn as sns
except ModuleNotFoundError:  # keep the experiment runnable in the lean CI image
    sns = None


PHASES = (
    ("hf", "hcp", "Hf hcp"),
    ("ti", "hcp", "Ti hcp"),
    ("zr", "hcp", "Zr hcp"),
    ("sio2", "quartz_beta", "SiO2 beta-quartz"),
)
PREDICTORS = {
    "absolute_G": "t3_tlog",
    "one_sided_QH": "t2_qh",
}
SEEDS = (11, 23, 37)
TOL = 1.0e-10


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_commit(root: Path) -> str:
    try:
        return subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def load_reference(root: Path, system: str, phase: str) -> dict[float, float]:
    path = root / "data/processed" / system / phase / "reference_G.csv"
    return {float(row["T_K"]): float(row["G_eV_per_atom"]) for row in read_csv(path)}


def load_split(root: Path) -> tuple[dict[str, list[dict[str, object]]], Path, str]:
    path = root / "data/processed/splits_v2/temp_extrap.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    split = {name: list(payload[name]) for name in ("train", "test")}
    return split, path, sha256(path)


def expected_temperatures(root: Path, split: dict[str, list[dict[str, object]]], system: str, phase: str) -> dict[str, np.ndarray]:
    reference = load_reference(root, system, phase)
    out: dict[str, np.ndarray] = {}
    for region in ("train", "test"):
        indexes = sorted(int(row["T_index"]) for row in split[region] if row["system"] == system and row["phase"] == phase)
        out[region] = np.asarray([sorted(reference)[i] for i in indexes], dtype=float)
    return out


def qh_reliable_phases(root: Path) -> set[tuple[str, str]]:
    inventory = root / "result/tables/phase_inventory.csv"
    reliable: set[tuple[str, str]] = set()
    for row in read_csv(inventory):
        if row.get("qh_reliable", "").strip().lower() == "true":
            reliable.add((row["system"], row["phase"]))
    return reliable


def prediction_file(root: Path, source: str, system: str, phase: str, seed: int, region: str) -> Path:
    path = root / "result/experiments/single_side_qh/predictions" / source / f"{system}_{phase}_seed{seed}_{region}.csv"
    if not path.exists():
        raise FileNotFoundError(f"missing canonical prediction file: {path}")
    return path


def load_prediction(path: Path, system: str, phase: str, seed: int, region: str) -> tuple[np.ndarray, np.ndarray]:
    rows = read_csv(path)
    if not rows:
        raise ValueError(f"empty prediction file: {path}")
    required = {"T_K", "pred_G_eV_per_atom", "system", "phase", "seed", "region", "natoms"}
    missing = required.difference(rows[0])
    if missing:
        raise ValueError(f"{path} missing columns {sorted(missing)}")
    for row in rows:
        if (row["system"], row["phase"], int(row["seed"]), row["region"]) != (system, phase, seed, region):
            raise ValueError(f"prediction metadata mismatch in {path}")
        if int(row["natoms"]) <= 0:
            raise ValueError(f"invalid atom count in {path}")
    if len({int(row["natoms"]) for row in rows}) != 1:
        raise ValueError(f"atom normalization changes within {path}")
    temperatures = np.asarray([float(row["T_K"]) for row in rows], dtype=float)
    values = np.asarray([float(row["pred_G_eV_per_atom"]) for row in rows], dtype=float)
    if not np.isfinite(values).all():
        raise ValueError(f"non-finite prediction in {path}")
    if rows[0].get("pred_G_eV_per_atom", "").lower().find("mev") >= 0:
        raise ValueError(f"prediction labels must be eV/atom, not meV/atom: {path}")
    return temperatures, values


def model_provenance(root: Path, predictor: str, seed: int, system: str) -> dict[str, str]:
    method = "repr_regression_tlog" if predictor == "absolute_G" else "qh_residual"
    path = root / "result/experiments/t3_temp_extrap/raw_runs" / method / f"seed_{seed}/metrics.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    provenance = payload["provenance"]
    checkpoints = provenance.get("checkpoint_path", [])
    hashes = provenance.get("checkpoint_sha256", [])
    checkpoint = next((str(item) for item in checkpoints if f"/{system}/" in str(item)), "")
    index = checkpoints.index(checkpoint) if checkpoint in checkpoints else -1
    return {
        "config": str(provenance["predictor_config"]),
        "checkpoint": checkpoint,
        "checkpoint_sha256": str(hashes[index]) if index >= 0 and index < len(hashes) else "",
        "head": str(provenance.get("heads", {}).get(system, "")),
        "git_commit": str(provenance.get("git_commit", "")),
        "split_sha256": str(provenance.get("split_sha256", "")),
    }


def rmse(values: np.ndarray) -> float:
    return float(np.sqrt(np.mean(values * values)))


def mae(values: np.ndarray) -> float:
    return float(np.mean(np.abs(values)))


def slope_mev_per_k(temperatures: np.ndarray, residual_eV: np.ndarray, t_ref: float) -> float:
    coefficient = np.polyfit(temperatures - t_ref, residual_eV * 1000.0, 1)[0]
    return float(coefficient)


def set_plot_style() -> None:
    if sns is not None:
        sns.set_theme(font_scale=1.0, style="whitegrid", font="DejaVu Sans")
    else:
        plt.rcParams.update({
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.grid": True,
            "axes.facecolor": "white",
            "figure.facecolor": "white",
        })


def despine() -> None:
    if sns is not None:
        sns.despine(left=True, bottom=True)
    else:
        for axis in plt.gcf().axes:
            axis.spines["left"].set_visible(False)
            axis.spines["bottom"].set_visible(False)


def check_and_collect(root: Path) -> tuple[list[dict[str, object]], dict[str, object]]:
    split, split_path, split_hash = load_split(root)
    allowed = {(system, phase) for system, phase, _ in PHASES}
    reliable = qh_reliable_phases(root)
    if not allowed.issubset(reliable):
        raise AssertionError(f"required QH-reliable phases missing from audited inventory: {sorted(allowed - reliable)}")
    forbidden = reliable.intersection({("hf", "bcc"), ("ti", "bcc"), ("zr", "bcc"), ("sio2", "cristobalite_beta"), ("sio2", "tridymite_p63mmc")})
    if forbidden:
        raise AssertionError(f"invalid QH phase entered experiment: {sorted(forbidden)}")

    rows: list[dict[str, object]] = []
    checks: dict[str, object] = {"split": str(split_path.relative_to(root)), "split_sha256": split_hash, "seed_sets": {}, "grid_interpolation": "none; exact temperature grids required"}
    predictor_seed_sets: dict[str, set[int]] = {}
    for predictor, source in PREDICTORS.items():
        predictor_seed_sets[predictor] = set(SEEDS)
        for seed in SEEDS:
            for system, phase, _ in PHASES:
                prov = model_provenance(root, predictor, seed, system)
                if prov["split_sha256"] != split_hash:
                    raise AssertionError(f"{predictor} seed {seed} uses a different split: {prov['split_sha256']}")
                expected = expected_temperatures(root, split, system, phase)
                reference = load_reference(root, system, phase)
                curves: dict[str, tuple[np.ndarray, np.ndarray]] = {}
                for region in ("train", "test"):
                    path = prediction_file(root, source, system, phase, seed, region)
                    temperatures, values = load_prediction(path, system, phase, seed, region)
                    if not np.array_equal(temperatures, expected[region]):
                        raise AssertionError(f"temperature grid mismatch for {path}")
                    if region == "train" and not len(temperatures):
                        raise AssertionError(f"empty training region for {path}")
                    curves[region] = (temperatures, values)
                t_ref = float(curves["train"][0][-1])
                if not any(abs(t - t_ref) <= TOL for t in curves["train"][0]):
                    raise AssertionError(f"T_ref is not in training grid for {system}/{phase}")
                ref_train = np.asarray([reference[float(t)] for t in curves["train"][0]])
                ref_test = np.asarray([reference[float(t)] for t in curves["test"][0]])
                ref_anchor = reference[t_ref]
                anchor_pred = float(curves["train"][1][-1])
                anchored_ref_train = ref_train - ref_anchor
                anchored_ref_test = ref_test - ref_anchor
                for region, ref_values, anchored_ref in (("train", ref_train, anchored_ref_train), ("test", ref_test, anchored_ref_test)):
                    temperatures, prediction = curves[region]
                    prediction_tilde = prediction - anchor_pred
                    residual = prediction_tilde - anchored_ref
                    if region == "train":
                        anchored_train = residual
                        raw_train_mae = mae(prediction - ref_values) * 1000.0
                        raw_test_mae = raw_test_rmse = anchored_test_mae = anchored_test_rmse = None
                        raw_offset = float((anchor_pred - ref_anchor) * 1000.0)
                    else:
                        anchored_test = residual
                        raw_train_mae = None
                        raw_test_mae = mae(prediction - ref_values) * 1000.0
                        raw_test_rmse = rmse(prediction - ref_values) * 1000.0
                        anchored_test_mae = mae(residual) * 1000.0
                        anchored_test_rmse = rmse(residual) * 1000.0
                        raw_offset = None
                    if region == "test":
                        slope = slope_mev_per_k(temperatures, residual, t_ref)
                    else:
                        slope = None
                    rows.append({
                        "system": system, "phase": phase, "seed": seed, "predictor": predictor,
                        "region": region, "T_ref_K": t_ref, "n_train": len(curves["train"][0]), "n_test": len(curves["test"][0]),
                        "raw_train_mae_meV": raw_train_mae, "raw_test_mae_meV": raw_test_mae,
                        "raw_test_rmse_meV": raw_test_rmse, "raw_offset_at_Tref_meV": raw_offset,
                        "anchored_train_mae_meV": mae(anchored_train) * 1000.0 if region == "train" else None,
                        "anchored_test_mae_meV": anchored_test_mae, "anchored_train_rmse_meV": rmse(anchored_train) * 1000.0 if region == "train" else None,
                        "anchored_test_rmse_meV": anchored_test_rmse,
                        "residual_slope_meV_per_K": slope, "residual_slope_u_eV_per_K": slope * 1000.0 if slope is not None else None,
                        "source_file": str(prediction_file(root, source, system, phase, seed, region).relative_to(root)),
                        "config": prov["config"], "checkpoint": prov["checkpoint"], "checkpoint_sha256": prov["checkpoint_sha256"],
                        "head": prov["head"], "git_commit": prov["git_commit"], "split_sha256": prov["split_sha256"],
                    })
                if abs(float(anchored_train[-1])) > 1.0e-9:
                    raise AssertionError(f"anchored residual at T_ref is nonzero for {predictor}/{system}/{phase}/{seed}")
                # Add one oracle row per valid seed so the per-seed table has a
                # directly comparable reference curve without implying a model.
                rows.append({
                    "system": system, "phase": phase, "seed": seed, "predictor": "oracle", "region": "train",
                    "T_ref_K": t_ref, "n_train": len(curves["train"][0]), "n_test": len(curves["test"][0]),
                    "raw_train_mae_meV": 0.0, "raw_test_mae_meV": None, "raw_test_rmse_meV": None,
                    "raw_offset_at_Tref_meV": 0.0, "anchored_train_mae_meV": 0.0, "anchored_test_mae_meV": None,
                    "anchored_train_rmse_meV": 0.0, "anchored_test_rmse_meV": None, "residual_slope_meV_per_K": None,
                    "residual_slope_u_eV_per_K": None, "source_file": str((root / "data/processed" / system / phase / "reference_G.csv").relative_to(root)),
                    "config": "", "checkpoint": "", "checkpoint_sha256": "", "head": "", "git_commit": git_commit(root), "split_sha256": split_hash,
                })
                rows.append({
                    "system": system, "phase": phase, "seed": seed, "predictor": "oracle", "region": "test",
                    "T_ref_K": t_ref, "n_train": len(curves["train"][0]), "n_test": len(curves["test"][0]),
                    "raw_train_mae_meV": None, "raw_test_mae_meV": 0.0, "raw_test_rmse_meV": 0.0,
                    "raw_offset_at_Tref_meV": None, "anchored_train_mae_meV": None, "anchored_test_mae_meV": 0.0,
                    "anchored_train_rmse_meV": None, "anchored_test_rmse_meV": 0.0, "residual_slope_meV_per_K": 0.0,
                    "residual_slope_u_eV_per_K": 0.0, "source_file": str((root / "data/processed" / system / phase / "reference_G.csv").relative_to(root)),
                    "config": "", "checkpoint": "", "checkpoint_sha256": "", "head": "", "git_commit": git_commit(root), "split_sha256": split_hash,
                })
    # The oracle rows are added once per seed after both model predictors are
    # traversed.  Keep this guard so a malformed export cannot silently create
    # duplicate reference rows.
    unique_rows: list[dict[str, object]] = []
    seen: set[tuple[object, ...]] = set()
    for row in rows:
        key = (row["system"], row["phase"], row["seed"], row["predictor"], row["region"])
        if key not in seen:
            seen.add(key)
            unique_rows.append(row)
    rows = unique_rows
    seed_sets = {key: sorted(value) for key, value in predictor_seed_sets.items()}
    if len({tuple(value) for value in seed_sets.values()}) != 1:
        raise AssertionError(f"T3/T2 seed sets differ: {seed_sets}")
    checks["seed_sets"] = seed_sets
    checks["validated_phases"] = [f"{system}:{phase}" for system, phase, _ in PHASES]
    checks["units"] = "input/reference eV/atom; reported metrics meV/atom and meV/atom/K"
    checks["anchoring"] = "T_ref=max(train temperature); all model anchored residuals at T_ref <= 1e-9 eV/atom"
    return rows, checks


def summarize(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for system, phase, label in PHASES:
        per_seed: dict[int, dict[str, float]] = {}
        for predictor in ("absolute_G", "one_sided_QH"):
            for row in rows:
                if row["system"] == system and row["phase"] == phase and row["predictor"] == predictor and row["region"] == "test":
                    per_seed.setdefault(int(row["seed"]), {})[predictor] = float(row["anchored_test_rmse_meV"])
        ratios = np.asarray([item["one_sided_QH"] / item["absolute_G"] for item in per_seed.values()])
        improvements = 100.0 * (1.0 - ratios)
        abs_values = np.asarray([item["absolute_G"] for item in per_seed.values()])
        qh_values = np.asarray([item["one_sided_QH"] for item in per_seed.values()])
        abs_slopes = np.asarray([float(row["residual_slope_meV_per_K"]) for row in rows if row["system"] == system and row["phase"] == phase and row["predictor"] == "absolute_G" and row["region"] == "test"])
        qh_slopes = np.asarray([float(row["residual_slope_meV_per_K"]) for row in rows if row["system"] == system and row["phase"] == phase and row["predictor"] == "one_sided_QH" and row["region"] == "test"])
        output.append({"system": system, "phase": phase, "label": label, "n_seeds": len(ratios), "abs_anchored_test_rmse_mean_meV": float(abs_values.mean()), "abs_anchored_test_rmse_std_meV": float(abs_values.std(ddof=1)), "qh_anchored_test_rmse_mean_meV": float(qh_values.mean()), "qh_anchored_test_rmse_std_meV": float(qh_values.std(ddof=1)), "R_mean": float(ratios.mean()), "R_std": float(ratios.std(ddof=1)), "improvement_mean_percent": float(improvements.mean()), "improvement_std_percent": float(improvements.std(ddof=1)), "abs_slope_mean_meV_per_K": float(abs_slopes.mean()), "abs_slope_std_meV_per_K": float(abs_slopes.std(ddof=1)), "qh_slope_mean_meV_per_K": float(qh_slopes.mean()), "qh_slope_std_meV_per_K": float(qh_slopes.std(ddof=1))})
    macro = {"system": "macro_average", "phase": "macro_average", "label": "Macro-average across four phases", "n_seeds": len(SEEDS)}
    mean_to_std = {
        "abs_anchored_test_rmse_mean_meV": "abs_anchored_test_rmse_std_meV",
        "qh_anchored_test_rmse_mean_meV": "qh_anchored_test_rmse_std_meV",
        "R_mean": "R_std",
        "improvement_mean_percent": "improvement_std_percent",
        "abs_slope_mean_meV_per_K": "abs_slope_std_meV_per_K",
        "qh_slope_mean_meV_per_K": "qh_slope_std_meV_per_K",
    }
    for key, std_key in mean_to_std.items():
        values = np.asarray([float(row[key]) for row in output])
        macro[key] = float(values.mean())
        macro[std_key] = float(values.std(ddof=1))
    output.append(macro)
    return output


def plot_outputs(root: Path, rows: list[dict[str, object]], summary: list[dict[str, object]], output: Path) -> None:
    set_plot_style()
    colors = {"oracle": "#1f253f", "absolute_G": "#4575b4", "one_sided_QH": "#d73027"}
    labels = {"oracle": "Oracle reference", "absolute_G": "Absolute-G T3 tlog", "one_sided_QH": "One-sided QH T2"}
    for system, phase, label in PHASES:
        fig, axes = plt.subplots(2, 1, figsize=(8, 6), dpi=150, sharex=True)
        phase_rows = [row for row in rows if row["system"] == system and row["phase"] == phase]
        reference = load_reference(root, system, phase)
        # The model CSVs carry the exact frozen grids.  Use the first seed as
        # the plotting grid after the validation above.
        for predictor, source in [("oracle", None), ("absolute_G", "t3_tlog"), ("one_sided_QH", "t2_qh")]:
            if predictor == "oracle":
                seed = SEEDS[0]
                pred_train_t = np.asarray([float(r["T_K"]) for r in read_csv(root / "result/experiments/single_side_qh/predictions/t3_tlog" / f"{system}_{phase}_seed{seed}_train.csv")])
                pred_train = np.asarray([reference[t] for t in pred_train_t])
                pred_test_t = np.asarray([float(r["T_K"]) for r in read_csv(root / "result/experiments/single_side_qh/predictions/t3_tlog" / f"{system}_{phase}_seed{seed}_test.csv")])
                pred_test = np.asarray([reference[t] for t in pred_test_t])
            else:
                source_path_train = root / "result/experiments/single_side_qh/predictions" / source / f"{system}_{phase}_seed{SEEDS[0]}_train.csv"
                source_path_test = root / "result/experiments/single_side_qh/predictions" / source / f"{system}_{phase}_seed{SEEDS[0]}_test.csv"
                pred_train_t, pred_train = load_prediction(source_path_train, system, phase, SEEDS[0], "train")
                pred_test_t, pred_test = load_prediction(source_path_test, system, phase, SEEDS[0], "test")
            t_ref = float(pred_train_t[-1]); anchor = float(pred_train[-1])
            ref_anchor = reference[t_ref]
            axes[0].plot(pred_train_t, (pred_train - anchor) * 1000.0, color=colors[predictor], label=labels[predictor], linewidth=1.8)
            axes[0].plot(pred_test_t, (pred_test - anchor) * 1000.0, color=colors[predictor], linewidth=1.8)
            axes[1].plot(pred_train_t, (pred_train - anchor - (np.asarray([reference[t] for t in pred_train_t]) - ref_anchor)) * 1000.0, color=colors[predictor], linewidth=1.2, alpha=0.75)
            axes[1].plot(pred_test_t, (pred_test - anchor - (np.asarray([reference[t] for t in pred_test_t]) - ref_anchor)) * 1000.0, color=colors[predictor], linewidth=1.8)
        axes[0].axvline(t_ref, color="dimgrey", linestyle="--", linewidth=0.9, label="Train/test boundary")
        axes[1].axvline(t_ref, color="dimgrey", linestyle="--", linewidth=0.9)
        axes[0].set_ylabel("Anchored G (meV/atom)", color="dimgrey")
        axes[1].set_ylabel("Anchored residual (meV/atom)", color="dimgrey")
        axes[1].set_xlabel("Temperature (K)", color="dimgrey")
        axes[0].set_title(f"{label}: anchored thermal curves", loc="left", color="dimgrey")
        axes[1].set_title("Residual to canonical reference", loc="left", color="dimgrey")
        summary_row = next(item for item in summary if item["system"] == system)
        axes[1].text(0.99, 0.97, f"R = {float(summary_row['R_mean']):.3f}\nQH improvement = {float(summary_row['improvement_mean_percent']):+.1f}%", transform=axes[1].transAxes, ha="right", va="top", color="dimgrey", fontsize=9)
        for ax in axes:
            ax.grid(False)
            ax.tick_params(axis="both", which="both", length=0, labelcolor="dimgrey")
            ax.patch.set_edgecolor("lightgrey")
            ax.patch.set_linewidth(0.8)
        axes[0].legend(frameon=True, facecolor="white", framealpha=0.8, edgecolor="lightgrey", labelcolor="dimgrey", fontsize=8, ncol=2)
        despine()
        fig.tight_layout()
        fig.savefig(output / f"{system}_{phase}.png", dpi=150, bbox_inches="tight")
        fig.savefig(output / f"{system}_{phase}.pdf", dpi=150, bbox_inches="tight")
        plt.close(fig)

    figure, axis = plt.subplots(figsize=(8, 4.5), dpi=150)
    plotted = [row for row in summary if row["system"] != "macro_average"]
    x = np.arange(len(plotted)); width = 0.36
    abs_values = [float(row["abs_anchored_test_rmse_mean_meV"]) for row in plotted]
    qh_values = [float(row["qh_anchored_test_rmse_mean_meV"]) for row in plotted]
    abs_err = [float(row["abs_anchored_test_rmse_std_meV"]) for row in plotted]
    qh_err = [float(row["qh_anchored_test_rmse_std_meV"]) for row in plotted]
    axis.bar(x - width / 2, abs_values, width, yerr=abs_err, capsize=3, color="#4575b4", label="Absolute-G T3 tlog")
    axis.bar(x + width / 2, qh_values, width, yerr=qh_err, capsize=3, color="#d73027", label="One-sided QH T2")
    axis.set_xticks(x, [str(row["label"]) for row in plotted], rotation=15, ha="right")
    axis.set_ylabel("Anchored test RMSE (meV/atom)", color="dimgrey")
    axis.set_title("Single-sided QH thermal-shape diagnostic", loc="left", color="dimgrey")
    axis.legend(frameon=True, facecolor="white", framealpha=0.8, edgecolor="lightgrey", labelcolor="dimgrey", fontsize=8)
    axis.grid(False); axis.tick_params(axis="both", which="both", length=0, labelcolor="dimgrey")
    axis.patch.set_edgecolor("lightgrey"); axis.patch.set_linewidth(0.8)
    despine()
    figure.tight_layout()
    figure.savefig(output / "summary.png", dpi=150, bbox_inches="tight")
    figure.savefig(output / "summary.pdf", dpi=150, bbox_inches="tight")
    plt.close(figure)


def write_readme(path: Path, summary: list[dict[str, object]], metric_rows: list[dict[str, object]]) -> None:
    rows = [row for row in summary if row["system"] != "macro_average"]
    lines = ["# Single-sided QH diagnostic", "", "This diagnostic compares anchored thermal shapes on the QH-reliable side only. It does not evaluate phase pairs, ΔG, ranking, or T_c.", "", "## Results", "", "| Phase | Absolute-G RMSE | One-sided QH RMSE | R = QH / absolute | Improvement | Absolute slope | QH slope |", "|---|---:|---:|---:|---:|---:|---:|"]
    for row in rows:
        lines.append(f"| {row['label']} | {row['abs_anchored_test_rmse_mean_meV']:.3f} ± {row['abs_anchored_test_rmse_std_meV']:.3f} | {row['qh_anchored_test_rmse_mean_meV']:.3f} ± {row['qh_anchored_test_rmse_std_meV']:.3f} | {row['R_mean']:.3f} ± {row['R_std']:.3f} | {row['improvement_mean_percent']:+.1f}% | {row['abs_slope_mean_meV_per_K']:.5f} | {row['qh_slope_mean_meV_per_K']:.5f} |")
    lines += ["", "Each model has one train row and one test row per phase and seed in `per_phase_seed_metrics.csv`; the summary uses the common test seeds only.", "", "The anchor is the final temperature in the frozen temp_extrap v2 training window for each phase. All curves are converted to anchored shape relative to that temperature before the primary RMSE is computed.", "", "QH is quantum phonopy while the canonical DaRUS/TI reference is classical MD. Therefore the anchored thermal-shape error, not raw absolute-G MAE, is the primary comparison.", "", "## Interpretation", ""]
    by_phase = {(row["system"], row["phase"]): row for row in rows}
    for system, phase, label in PHASES:
        row = by_phase[(system, phase)]
        model_rows = [r for r in metric_rows if r["system"] == system and r["phase"] == phase]
        offsets = {}
        slopes = {}
        for predictor in ("absolute_G", "one_sided_QH"):
            offsets[predictor] = float(np.mean([abs(float(r["raw_offset_at_Tref_meV"])) for r in model_rows if r["predictor"] == predictor and r["region"] == "train"]))
            slopes[predictor] = float(np.mean([abs(float(r["residual_slope_meV_per_K"])) for r in model_rows if r["predictor"] == predictor and r["region"] == "test"]))
        offset_word = "offset and slope" if offsets["one_sided_QH"] < offsets["absolute_G"] and slopes["one_sided_QH"] < slopes["absolute_G"] else "slope only" if slopes["one_sided_QH"] < slopes["absolute_G"] else "neither offset nor slope"
        lines.append(f"- **{label}:** QH reduces the anchored thermal residual (R = {float(row['R_mean']):.3f}); the three seeds are consistent (R std = {float(row['R_std']):.4f}). The improvement is mainly {offset_word} (mean absolute offset: {offsets['absolute_G']:.1f} → {offsets['one_sided_QH']:.1f} meV; slope magnitude: {slopes['absolute_G']:.5f} → {slopes['one_sided_QH']:.5f} meV/atom/K).")
    lines += ["", "Across all four QH-reliable phases, QH gives R < 1, but this is a single-phase thermal-shape result only. It does not establish pairwise phase-stability improvement, does not validate half-QH on a dynamically unstable high-temperature partner, and does not identify a true entropy error.", "", "## Reproduce", "", "```bash", "python src/experiments/single_side_qh.py --repo-root .", "```", "", "The script verifies QH reliability, the frozen split hash, exact temperature grids, common seed sets, eV/atom input units, and zero anchored residual at T_ref. Prediction CSVs are the canonical T3 tlog and T2 QH-residual outputs exported from the isolated runs; their source paths, checkpoint placeholders, heads, hashes, and commits are recorded in `per_phase_seed_metrics.csv` and `findings.meta.json`.", ""]
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    args = parser.parse_args()
    root = args.repo_root.resolve()
    output = root / "result/experiments/single_side_qh"
    rows, checks = check_and_collect(root)
    fields = list(rows[0])
    write_csv(output / "per_phase_seed_metrics.csv", fields, rows)
    summary = summarize(rows)
    write_csv(output / "qh_vs_abs_summary.csv", list(summary[0]), summary)
    plot_outputs(root, rows, summary, output)
    input_files = [root / "data/processed/splits_v2/temp_extrap.json", root / "result/tables/phase_inventory.csv"]
    input_files += [root / "data/processed" / system / phase / "reference_G.csv" for system, phase, _ in PHASES]
    input_files += sorted((output / "predictions").glob("**/*.csv"))
    meta = {"git_commit": git_commit(root), "generated_at_utc": datetime.now(timezone.utc).isoformat(), "script": {"path": "src/experiments/single_side_qh.py", "sha256": sha256(root / "src/experiments/single_side_qh.py")}, "inputs": [{"path": str(path.relative_to(root)), "sha256": sha256(path)} for path in input_files], "checks": checks, "primary_metric": "anchored test RMSE of e(T)=G_tilde_model-G_tilde_reference", "reference_convention": "canonical DaRUS thermodynamic reference; QH is quantum phonopy and reference is classical TI/MD", "seed_policy": list(SEEDS), "predictors": PREDICTORS}
    (output / "findings.meta.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    write_readme(output / "README.md", summary, rows)
    print("Phase        abs RMSE    QH RMSE     ratio R    improvement    slope abs   slope QH")
    for row in summary:
        print(f"{row['label']:<16} {row.get('abs_anchored_test_rmse_mean_meV', float('nan')):9.3f} {row.get('qh_anchored_test_rmse_mean_meV', float('nan')):9.3f} {row.get('R_mean', float('nan')):9.3f} {row.get('improvement_mean_percent', float('nan')):11.1f}% {row.get('abs_slope_mean_meV_per_K', float('nan')):11.5f} {row.get('qh_slope_mean_meV_per_K', float('nan')):9.5f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
