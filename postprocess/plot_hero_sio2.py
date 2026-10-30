"""Hero figure: SiO2 held-out cristobalite, three-source readout, DPA (B) vs DPA+volume (C).

Real archived predictions only. Curves come from the post-hoc regeneration of the frozen
A/B/C run (fes-decision-metric-audit-20261030/curves/SIO2_ABC_TASKCKPT.npz), which matches
the archived gap MAE to 4.5e-11 relative and all archived event lists exactly. This script
re-checks the plotted model against the frozen PRIMARY_COMPARISON.csv and EVENTS.csv before
drawing, and writes the plotted values as source data.
Usage: plot_hero_sio2.py <curves.npz> <abc_run_dir> <sio2_reference_phase_dir> <out_dir>
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

NPZ, RUN, REF, OUT = (Path(a) for a in sys.argv[1:5])
OUT.mkdir(parents=True, exist_ok=True)
POOL = ["quartz_beta", "tridymite_p63mmc", "cristobalite_beta"]
KEY = "cristobalite_beta__high_all_three__equal_block__{arm}__{p}"
NAMES = {"quartz_beta": "β-quartz", "tridymite_p63mmc": "β-tridymite", "cristobalite_beta": "β-cristobalite"}
PHASE_COLOR = {"quartz_beta": "#c7d5df", "tridymite_p63mmc": "#9a7bb6", "cristobalite_beta": "#5aa469"}
INK, TEAL, AMBER = "#203546", "#087f80", "#a96613"

c = np.load(NPZ); T = c["T_K"]
ref = {}
for p in POOL:
    rows = list(csv.DictReader(open(REF / p / "reference_G.csv")))
    assert np.array_equal(np.array([float(r["T_K"]) for r in rows]), T)
    ref[p] = np.array([float(r["G_eV_per_atom"]) for r in rows])
pred = {arm: {p: c[KEY.format(arm=arm, p=p)] for p in POOL} for arm in ("B", "C")}


def envelope(curves):
    lab = np.array(POOL)[np.argmin(np.column_stack([curves[p] for p in POOL]), axis=1)]
    ev = []
    for i in np.flatnonzero(lab[:-1] != lab[1:]):  # same interpolation rule as the frozen scorer
        a, b = str(lab[i]), str(lab[i + 1]); d0 = curves[a][i] - curves[b][i]; d1 = curves[a][i + 1] - curves[b][i + 1]
        f = abs(d0) / (abs(d0) + abs(d1)) if d0 * d1 < 0 else .5
        ev.append((a, b, float(T[i] + f * (T[i + 1] - T[i]))))
    return lab, ev


ref_lab, ref_ev = envelope(ref)
lab = {arm: envelope(pred[arm]) for arm in ("B", "C")}

# ---- consistency with the frozen run ----
prim = {r["arm"]: r for r in csv.DictReader(open(RUN / "PRIMARY_COMPARISON.csv")) if r["n_sources"] == "3" and r["scaling"] == "equal_block"}
for arm in ("B", "C"):
    acc = 100 * np.mean(lab[arm][0] == ref_lab)
    assert abs(acc - float(prim[arm]["stable_accuracy_percent"])) < 1e-6, (arm, acc)
    gaps = []
    for anc in ("quartz_beta", "tridymite_p63mmc"):
        e = (pred[arm]["cristobalite_beta"] - pred[arm][anc]) - (ref["cristobalite_beta"] - ref[anc])
        gaps.append(np.mean(abs(e)) * 1000)
    assert abs(np.mean(gaps) - float(prim[arm]["gap_MAE_meV_atom"])) < 1e-6, (arm, gaps)
for r in csv.DictReader(open(RUN / "EVENTS.csv")):
    if r["target"] == "cristobalite_beta" and r["configuration"] == "high_all_three" and r["scaling"] == "equal_block" and r["arm"] in ("B", "C") and r["candidate_pool"] == ";".join(POOL):
        arch = [(e["from"], e["to"]) for e in json.loads(r["predicted"])]
        assert arch == [(a, b) for a, b, _ in lab[r["arm"]][1]], (r["arm"], arch)

# ---- figure ----
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.edgecolor": "#8aa0ae", "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK})
fig = plt.figure(figsize=(8.6, 5.4), dpi=200)
gs = fig.add_gridspec(2, 1, height_ratios=[3.0, 1.15], hspace=0.32)
ax = fig.add_subplot(gs[0]); sx = fig.add_subplot(gs[1], sharex=ax)
mev = lambda a, b, d: (d[a] - d[b]) * 1000
# Margin of the held-out phase to its best competitor; below zero = cristobalite is the stable phase.
margin = lambda d: (d["cristobalite_beta"] - np.minimum(d["quartz_beta"], d["tridymite_p63mmc"])) * 1000
ax.axhspan(-3, 0, color="#5aa469", alpha=0.07, lw=0)
ax.axhline(0, color="#8aa0ae", lw=0.8)
ax.text(T[0] + 15, -0.08, "below 0: cristobalite is the stable phase", ha="left", va="top", fontsize=8, color="#3d7a4a")
ax.plot(T, margin(ref), color=INK, lw=2.6, label="DFT reference (PBE-D3(BJ))", zorder=4)
ax.plot(T, margin(pred["B"]), color=TEAL, lw=2.0, ls="--", label="DPA readout", zorder=5)
ax.plot(T, margin(pred["C"]), color=AMBER, lw=2.2, ls=(0, (1.0, 1.3)), label="DPA + volume readout", zorder=5)
T0 = [t for a, b, t in ref_ev if b == "cristobalite_beta"][0]
TB = [t for a, b, t in lab["B"][1] if b == "cristobalite_beta"][0]
ax.axvline(T0, color=INK, lw=0.8, ls=":")
ax.vlines(TB, -0.8, 0.0, color=TEAL, lw=0.8, ls=":")
ax.annotate(f"reference\n{T0:.0f} K", xy=(T0, 0), xytext=(T0 - 190, -0.5), fontsize=8.5, color=INK, ha="center", va="center",
            arrowprops=dict(arrowstyle="->", color=INK, lw=0.7))
ax.annotate(f"DPA: {TB:.0f} K, after a\ntridymite interval", xy=(TB, 0), xytext=(TB + 190, 0.95), fontsize=8.5, color=TEAL, ha="center", va="center",
            arrowprops=dict(arrowstyle="->", color=TEAL, lw=0.7))
ax.annotate("DPA + volume:\nnever stable", xy=(2400, float(np.interp(2400, T, margin(pred["C"])))), xytext=(2330, 0.95), fontsize=8.5, color=AMBER, ha="center", va="center",
            arrowprops=dict(arrowstyle="->", color=AMBER, lw=0.7))
ax.set_ylabel("G(cristobalite) − min G(quartz, tridymite)\n(meV/atom)", fontsize=9)
ax.set_ylim(-0.8, 2.4)
ax.legend(fontsize=8, frameon=False, loc="upper center", bbox_to_anchor=(0.66, 1.0))
ax.tick_params(labelbottom=False)
for s in ("top", "right"): ax.spines[s].set_visible(False)
rowsy = {"DFT reference": 2, "DPA": 1, "DPA + volume": 0}
for name, labs in (("DFT reference", ref_lab), ("DPA", lab["B"][0]), ("DPA + volume", lab["C"][0])):
    y = rowsy[name]
    starts = np.r_[0, np.flatnonzero(labs[1:] != labs[:-1]) + 1]; ends = np.r_[starts[1:], len(T)]
    for s, e in zip(starts, ends):
        sx.add_patch(plt.Rectangle((T[s], y - 0.36), T[e - 1] - T[s] + 1, 0.72, color=PHASE_COLOR[labs[s]], lw=0))
sx.set_yticks([2, 1, 0], ["DFT reference", "DPA", "DPA + volume"], fontsize=8.5)
sx.set_ylim(-0.6, 2.6); sx.set_xlim(T[0], T[-1])
sx.axvline(T0, color=INK, lw=0.8, ls=":")
sx.set_xlabel("Temperature (K)")
sx.set_title("Stable phase among quartz, tridymite and cristobalite", fontsize=8.5, loc="left", color=INK, pad=4)
for s in ("top", "right", "left"): sx.spines[s].set_visible(False)
sx.tick_params(axis="y", length=0)
handles = [plt.Rectangle((0, 0), 1, 1, color=PHASE_COLOR[p]) for p in POOL]
fig.legend(handles, [NAMES[p] for p in POOL], fontsize=8, frameon=False, ncol=3, loc="lower center", bbox_to_anchor=(0.55, -0.045))
fig.align_ylabels([ax, sx])
for ext in ("png", "svg", "pdf"):
    fig.savefig(OUT / f"HERO_SIO2_DPA_vs_DPAV.{ext}", bbox_inches="tight")

with open(OUT / "HERO_SIO2_SOURCE_DATA.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["T_K", "ref_crst_margin_meV", "DPA_crst_margin_meV", "DPAV_crst_margin_meV",
                "ref_crst_minus_q_meV", "ref_trd_minus_q_meV", "DPA_crst_minus_q_meV", "DPA_trd_minus_q_meV", "DPAV_crst_minus_q_meV", "DPAV_trd_minus_q_meV",
                "ref_stable", "DPA_stable", "DPAV_stable"])
    for i in range(len(T)):
        w.writerow([T[i], margin(ref)[i], margin(pred["B"])[i], margin(pred["C"])[i],
                    mev("cristobalite_beta", "quartz_beta", ref)[i], mev("tridymite_p63mmc", "quartz_beta", ref)[i],
                    mev("cristobalite_beta", "quartz_beta", pred["B"])[i], mev("tridymite_p63mmc", "quartz_beta", pred["B"])[i],
                    mev("cristobalite_beta", "quartz_beta", pred["C"])[i], mev("tridymite_p63mmc", "quartz_beta", pred["C"])[i],
                    ref_lab[i], lab["B"][0][i], lab["C"][0][i]])
summary = {"model": "held-out cristobalite; sources quartz, P6_3/mmc tridymite, C222_1 tridymite; equal-block scaling; Si/O task checkpoint",
           "reference_events": ref_ev, "DPA_events": lab["B"][1], "DPAV_events": lab["C"][1],
           "frozen_metrics": {arm: {k: prim[arm][k] for k in ("gap_MAE_meV_atom", "relative_entropy_MAE_microeV_atom_K", "stable_accuracy_percent")} for arm in ("B", "C")}}
(OUT / "HERO_SIO2_SUMMARY.json").write_text(json.dumps(summary, indent=1) + "\n")
print(json.dumps(summary, indent=1))
