"""Fig. 5 of the v0.1 paper: representational capacity (gates 1b/1c, complete labels) and reference efficiency (RE1).
Reads results/GATE1C_A1.json and results/RE1.json only. Usage: make_fig5.py <out_dir>
"""
import json, sys
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

R = Path(__file__).parents[1] / "results"; out = Path(sys.argv[1])
COL = {"F0": "#0072B2", "MS": "#E69F00"}; NAME = {"F0": "DPA-3.1 (F0)", "MS": "MatterSim-v1"}; MK = {"F0": "o", "MS": "s"}
INK, MUTED, GRID = "#1a1a1a", "#555555", "#dddddd"
plt.rcParams.update({"font.size": 8, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False, "pdf.fonttype": 42})
c = json.load(open(R / "GATE1C_A1.json")); r = json.load(open(R / "RE1.json"))
fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.6), gridspec_kw=dict(wspace=0.35, width_ratios=[1, 1.25]))

labels = ["Uncorrected", "R3\n(24 par.)", "Full tensor\n(189 par.)"]
for k, M in enumerate(("F0", "MS")):
    b = c[M]; base = b["A1_proj_vs_base_proj"]["success_base"]; x = np.arange(3) + (k - 0.5) * 0.28
    vals = [base, b["R3_proj_vs_base_proj"]["success_var"], b["A1_proj_vs_base_proj"]["success_var"]]
    lo = [0, base + b["R3_proj_vs_base_proj"]["d_success"]["ci95"][0], base + b["A1_proj_vs_base_proj"]["d_success"]["ci95"][0]]
    hi = [0, base + b["R3_proj_vs_base_proj"]["d_success"]["ci95"][1], base + b["A1_proj_vs_base_proj"]["d_success"]["ci95"][1]]
    ax[0].bar(x, vals, width=0.26, color=COL[M], label=NAME[M])
    ax[0].errorbar(x[1:], vals[1:], yerr=[np.array(vals[1:]) - lo[1:], np.array(hi[1:]) - vals[1:]], fmt="none", ecolor=INK, lw=0.8, capsize=2)
ax[0].axhline(0.5, color=MUTED, lw=0.8, ls="-."); ax[0].text(-0.4, 0.515, "pre-specified threshold", fontsize=6.3, color=MUTED, ha="left", va="bottom")
ax[0].set_xticks(range(3)); ax[0].set_xticklabels(labels, fontsize=6.8); ax[0].set_ylim(0, 0.8); ax[0].set_ylabel("Crossing success (93 crossings)")
ax[0].legend(fontsize=6.5, loc="upper left", bbox_to_anchor=(0, 0.93)); ax[0].grid(axis="y", color=GRID, lw=0.4)
ax[0].set_title("a  Complete reference labels", loc="left", fontsize=8.5)

cost = r["reference_cost"]; frac = [cost[str(b)]["1"]["fraction"] for b in (0.1, 0.25, 0.5, 1.0)]
for M in ("F0", "MS"):
    m2 = [r[M][str(b)]["M2"]["success"] for b in (0.1, 0.25, 0.5, 1.0)]
    ax[1].plot(frac, m2, color=COL[M], marker=MK[M], ms=4, lw=1.4, label=f"{NAME[M]}: MLIP + residual")
    ax[1].axhline(r[M]["M0"]["success"], color=COL[M], lw=0.8, ls="--", label=f"{NAME[M]}: uncorrected")
m1 = {M: [r[M][str(b)]["M1"]["success"] for b in (0.1, 0.25, 0.5, 1.0)] for M in ("F0", "MS")}
assert np.allclose(m1["F0"], m1["MS"]), "the no-MLIP fit is model-independent; one line is drawn for both"
ax[1].plot(frac, m1["F0"], color=MUTED, marker="^", ms=3.5, mfc="white", lw=1.0, ls=":", label="Residual without MLIP (both models)")
ax[1].axhline(0.5, color=MUTED, lw=0.8, ls="-."); ax[1].axvline(0.25, color=GRID, lw=3, zorder=0)
ax[1].set_xscale("log"); ax[1].set_xticks([0.1, 0.25, 0.5, 1.0]); ax[1].set_xticklabels(["10%", "25%\n(primary)", "50%", "100%"], fontsize=6.8)
ax[1].set_ylim(0, 0.8); ax[1].set_xlabel("Reference displacements used (random subset)"); ax[1].set_ylabel("Crossing success")
ax[1].legend(fontsize=5.8, loc="upper left", bbox_to_anchor=(0, 0.97)); ax[1].grid(axis="y", color=GRID, lw=0.4)
ax[1].set_title("b  Reduced reference labels", loc="left", fontsize=8.5)
for ext in ("pdf", "png"): fig.savefig(out / f"fig5_curvature_capacity_efficiency.{ext}", bbox_inches="tight", dpi=200)
print("wrote", out)
