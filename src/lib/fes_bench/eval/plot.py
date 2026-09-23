"""Render reference and predictor Delta-G curves for frozen benchmark splits."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np

try:
    import seaborn as sns
except ImportError:  # V100 benchmark environment deliberately has no seaborn.
    sns = None

from fes_bench.eval.run import _predictor, _table_cache


def _phase_label(phase: str) -> str:
    return {
        "quartz_beta": "β-quartz",
        "cristobalite_beta": "β-cristobalite",
        "tridymite_p63mmc": "P6₃/mmc tridymite",
    }.get(phase, phase)


def _views(split: dict[str, object], requested_fold: str) -> list[tuple[str, dict[str, object]]]:
    if "folds" not in split:
        if requested_fold != "all":
            raise ValueError("--fold is only valid for a nested split")
        return [(str(split.get("name", "split")), split)]
    raw_folds = split["folds"]
    if not isinstance(raw_folds, dict):
        raise ValueError("split folds must be an object")
    if requested_fold == "all":
        return [(str(name), fold) for name, fold in sorted(raw_folds.items()) if isinstance(fold, dict)]
    try:
        fold = raw_folds[requested_fold]
    except KeyError as exc:
        raise ValueError(f"unknown fold {requested_fold!r}") from exc
    if not isinstance(fold, dict):
        raise ValueError(f"fold {requested_fold!r} must be an object")
    return [(requested_fold, fold)]


def _curves_for_view(
    name: str,
    split: dict[str, object],
    data_root: Path,
    predictor: str,
    seeds: list[int],
) -> list[tuple[str, np.ndarray, np.ndarray, np.ndarray]]:
    test = split.get("test")
    train = split.get("train", [])
    if not isinstance(test, list) or not isinstance(train, list):
        raise ValueError("each view must provide list-valued train and test entries")
    rows = test + train
    cache = _table_cache(data_root, rows)
    grouped: dict[tuple[str, str], list[int]] = defaultdict(list)
    tested: set[tuple[str, str]] = set()
    for row in rows:
        grouped[(str(row["system"]), str(row["phase"]))].append(int(row["T_index"]))
    for row in test:
        tested.add((str(row["system"]), str(row["phase"])))
    curves: list[tuple[str, np.ndarray, np.ndarray, np.ndarray]] = []
    for system in sorted({key[0] for key in tested}):
        phases = sorted(phase for item_system, phase in grouped if item_system == system)
        values: dict[str, tuple[np.ndarray, np.ndarray]] = {}
        for phase in phases:
            indices = sorted(set(grouped[(system, phase)]))
            table = cache[(system, phase)]
            temperature = np.array([table[index].T_K for index in indices])
            reference = np.array([table[index].G_eV_per_atom for index in indices])
            predicted = (
                reference
                if predictor == "reference"
                else np.stack([_predictor(predictor, seed, f"{system}:{phase}", reference) for seed in seeds]).mean(axis=0)
            )
            values[phase] = temperature, predicted
        for index, left in enumerate(phases):
            for right in phases[index + 1 :]:
                if (system, left) not in tested and (system, right) not in tested:
                    continue
                left_t, left_g = values[left]
                right_t, right_g = values[right]
                shared = np.array(sorted(set(left_t.tolist()) & set(right_t.tolist())))
                left_map, right_map = dict(zip(left_t, left_g)), dict(zip(right_t, right_g))
                ref_left = {cache[(system, left)][i].T_K: cache[(system, left)][i].G_eV_per_atom for i in grouped[(system, left)]}
                ref_right = {cache[(system, right)][i].T_K: cache[(system, right)][i].G_eV_per_atom for i in grouped[(system, right)]}
                curves.append(
                    (
                        f"{name}: {_phase_label(left)} − {_phase_label(right)}",
                        shared,
                        np.array([ref_left[value] - ref_right[value] for value in shared]),
                        np.array([left_map[value] - right_map[value] for value in shared]),
                    )
                )
    return curves


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictor", required=True)
    parser.add_argument("--split", required=True)
    parser.add_argument("--data-root", required=True)
    parser.add_argument("--seeds", default="11,23,37,51,67")
    parser.add_argument("--fold", default="all", help="nested fold name, or all")
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    split = json.loads(Path(args.split).read_text(encoding="utf-8"))
    data_root = Path(args.data_root).resolve()
    seeds = [int(value) for value in args.seeds.split(",") if value]
    curves = [curve for name, view in _views(split, args.fold) for curve in _curves_for_view(name, view, data_root, args.predictor, seeds)]
    if not curves:
        raise ValueError("no comparable phase-pair curves in requested split")

    # --- Style Setup ---
    if sns is not None:
        sns.set_theme(font_scale=1.0, style="whitegrid", font="DejaVu Sans")
    else:
        plt.style.use("seaborn-v0_8-whitegrid")
        plt.rcParams["font.family"] = "DejaVu Sans"
    n_columns = 2 if len(curves) <= 6 else 3
    n_rows = int(np.ceil(len(curves) / n_columns))
    figure = plt.figure(figsize={(2, 2): (10, 8), (3, 2): (14, 9), (3, 3): (14, 12)}.get((n_rows, n_columns), (14, 12)), dpi=150)
    grid = gridspec.GridSpec(n_rows, n_columns, figure=figure)
    grid.update(wspace=0.10, hspace={2: 0.30, 3: 0.35}.get(n_rows, 0.40), left=0.08, right=0.99, top=0.92, bottom=0.08)
    axes = [figure.add_subplot(grid[row, column]) for row in range(n_rows) for column in range(n_columns)]
    all_values = np.concatenate([np.concatenate((reference, predicted)) for _, _, reference, predicted in curves]) * 1000.0
    pad = max(0.2, 0.08 * (float(all_values.max()) - float(all_values.min())))
    reference_color = "#1f253f"  # darkest documented cubehelix tone
    predicted_color = "#bd0c0c"
    max_difference = max(float(np.max(np.abs(reference - predicted))) for _, _, reference, predicted in curves) * 1000.0
    figure.suptitle("Frozen-split ΔG curves", fontsize=14, y=0.98, color="dimgrey")

    # --- Plot ---
    for panel, (axis, (title, temperature, reference, predicted)) in enumerate(zip(axes, curves)):
        row, column = divmod(panel, n_columns)
        axis.plot(temperature, reference * 1000.0, color=reference_color, lw=2.0, label="reference")
        axis.plot(
            temperature,
            predicted * 1000.0,
            color=predicted_color,
            lw=1.6,
            ls="--",
            label="predictor" if args.predictor == "reference" else args.predictor,
        )
        axis.axhline(0.0, color="0.35", lw=0.8, zorder=0)
        axis.set_ylim(float(all_values.min()) - pad, float(all_values.max()) + pad)
        axis.set_title(r"$\bf{(" + chr(ord("a") + panel) + r")}$  " + title, loc="left", fontsize=9, pad=7)
        axis.grid(False)
        axis.tick_params(axis="both", which="both", length=0, labelcolor="dimgrey")
        axis.patch.set_edgecolor("lightgrey")
        axis.patch.set_linewidth(0.8)
        if column == 0:
            axis.set_ylabel("ΔG (meV/atom)", color="dimgrey")
        else:
            axis.tick_params(labelleft=False)
        if row == n_rows - 1:
            axis.set_xlabel("Temperature (K)", color="dimgrey")
    for axis in axes[len(curves) :]:
        axis.set_visible(False)
    axes[0].legend(frameon=True, facecolor="white", framealpha=0.8, edgecolor="lightgrey", labelcolor="dimgrey", fontsize=8, ncol=2)
    figure.text(0.99, 0.01, f"Largest |ΔG prediction error| across panels: {max_difference:.3g} meV/atom", ha="right", va="bottom", fontsize=9, color="dimgrey", style="italic")
    for axis in axes[: len(curves)]:
        axis.spines[["top", "right", "left", "bottom"]].set_visible(False)

    # --- Save ---
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=150, bbox_inches="tight")
    figure.savefig(output.with_suffix(".pdf"), dpi=150, bbox_inches="tight")
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
