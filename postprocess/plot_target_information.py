"""S2 target-information intervention: how much target information repairs the held-out cristobalite prediction.

Archived S2 results only (sio2-information-intervention/data/S2_INFORMATION_LADDER_SUMMARY.csv); earlier
two-source LOPO readouts C-simple and C-DPA, PBE-D3(BJ) reference, held-out beta-cristobalite.
Categories are different information budgets, so they are drawn as grouped points, not a line.
Practical few-shot budgets (0, 1, 2 pre-declared target labels) are separated from the two oracle
diagnostics, which use reference information at or around the crossing.
Usage: plot_target_information.py <S2_INFORMATION_LADDER_SUMMARY.csv> <out_dir>
"""
import csv
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

SRC, OUT = Path(sys.argv[1]), Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)
INK, TEAL, SLATE, ORACLE_BG = "#203546", "#087f80", "#5f7d93", "#eef1f4"
ARMS = [("zero_shot", "Zero-shot\n0 labels"), ("few_shot_1_label", "One label\n1263 K"),
        ("few_shot_2_label", "Two labels\n1263, 2087 K"), ("oracle_gap", "Oracle gap\n1675 K"),
        ("oracle_gap_slope", "Oracle gap\n+ slope\n11 points")]
READOUT = {"C-simple": (SLATE, "o"), "C-DPA": (TEAL, "s")}
rows = {(r["method"], r["arm"]): r for r in csv.DictReader(open(SRC, encoding="utf-8"))}
val = {k: (float(r["target_quartz_gap_MAE_meV_atom"]), float(r["envelope_hit_Tc_MAE_K"]), int(r["target_label_budget"])) for k, r in rows.items()}

# Statements used in the page caption.
for m in READOUT:
    assert all(rows[(m, a)]["envelope_matched_events"] == "1" and rows[(m, a)]["envelope_false_events"] == "0" for a, _ in ARMS)
    assert val[(m, "few_shot_1_label")][0] > val[(m, "zero_shot")][0]          # one label does not guarantee a lower MAE
    assert val[(m, "few_shot_2_label")][0] < 0.1                               # two fixed labels repair the gap
    assert val[(m, "oracle_gap_slope")][1] < 1.3 and val[(m, "oracle_gap_slope")][2] == 11
assert [val[("C-simple", a)][2] for a, _ in ARMS] == [0, 1, 2, 1, 11]

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": "#8aa0ae", "axes.labelcolor": INK,
                     "xtick.color": INK, "ytick.color": INK})
fig, axes = plt.subplots(1, 2, figsize=(11.0, 3.9), dpi=200)
for ax, (j, ylab, fmt) in zip(axes, ((0, "Cristobalite − quartz gap MAE (meV/atom)", "{:.3f}"),
                                     (1, "Transition-temperature error (K)", "{:.1f}"))):
    ax.axvspan(2.5, 4.5, color=ORACLE_BG, lw=0, zorder=0)
    top = max(v[j] for v in val.values()) * 1.22
    ax.text(3.5, top * 0.985, "oracle diagnostics\n(use reference information)", ha="center", va="top", fontsize=7.6, color="#536b7b")
    ax.text(1.0, top * 0.985, "few-shot budgets", ha="center", va="top", fontsize=7.6, color="#536b7b")
    for k, (m, (col, mk)) in enumerate(READOUT.items()):
        dx = -0.13 if k == 0 else 0.13
        for i, (a, _) in enumerate(ARMS):
            v = val[(m, a)][j]; oracle = a.startswith("oracle")
            ax.plot(i + dx, v, mk, ms=7, color=col, mfc="white" if oracle else col, mew=1.6, zorder=3,
                    label=m if i == 0 else None)
            ax.text(i + dx + (-0.07 if k == 0 else 0.07), v + top * 0.025, (fmt if v >= 2 or j == 0 else "{:.2f}").format(v), ha="right" if k == 0 else "left",
                    va="bottom", fontsize=7, color=col)
    ax.set_xticks(range(len(ARMS)), [t for _, t in ARMS], fontsize=7.8)
    ax.set_xlim(-0.55, 4.55); ax.set_ylim(0, top)
    ax.set_ylabel(ylab, fontsize=8.8)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.legend(*axes[0].get_legend_handles_labels(), fontsize=8, frameon=False, ncol=2, loc="lower center", bbox_to_anchor=(0.5, -0.07))
fig.tight_layout(w_pad=2.5)
for ext in ("png", "svg", "pdf"):
    fig.savefig(OUT / f"TARGET_INFORMATION_S2.{ext}", bbox_inches="tight")
summary = {f"{m}|{a}": {"gap_MAE_meV_atom": val[(m, a)][0], "Tc_error_K": val[(m, a)][1], "target_labels": val[(m, a)][2]}
           for m in READOUT for a, _ in ARMS}
(OUT / "TARGET_INFORMATION_S2.json").write_text(json.dumps(summary, indent=1) + "\n")
print(json.dumps(summary, indent=1))
