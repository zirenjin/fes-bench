"""Anchor / integrated-entropy decomposition of the cristobalite gap error, DPA (B) vs DPA+volume (C).

Same model as the hero figure: held-out cristobalite, three source phases, equal-block scaling.
Plots the archived post-fit decomposition (fes-sio2-supervision-physics-goal-20261028/abc_volume,
ANCHOR_ENTROPY_ERROR_CURVES.csv), in which

    e_dG(T) = e_dG(T0) - integral_{T0}^{T} e_dS(t) dt,   T0 = reference quartz->cristobalite crossing.

This is an error identity, not a causal model. Before drawing, the script checks
  1. the archived error curves against the regenerated predictions and the DFT reference,
  2. the anchor / 2499 K / MAE values against ANCHOR_ENTROPY_PAIR_SUMMARY.csv,
  3. that "gap error below the reference margin" reproduces where each readout predicts cristobalite stable.
Usage: plot_error_decomposition.py <abc_volume_dir> <curves.npz> <sio2_reference_phase_dir> <out_dir>
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ABC, NPZ, REF, OUT = (Path(a) for a in sys.argv[1:5])
OUT.mkdir(parents=True, exist_ok=True)
POOL = ["quartz_beta", "tridymite_p63mmc", "cristobalite_beta"]
KEY = "cristobalite_beta__high_all_three__equal_block__{arm}__{p}"
REP = {"B": "B_DPA", "C": "C_DPA_plus_volume_equal_group"}
LABEL = {"B": "DPA", "C": "DPA + volume (equal-block)"}  # names as in the phase-decision figure
INK, TEAL, AMBER, GREEN = "#203546", "#087f80", "#a96613", "#5aa469"
COLOR = {"B": TEAL, "C": AMBER}

# Archived curves. The directed gap is written as two fields (first, second), so rows have one
# more field than the header names: representation, first, second, T_K, then the header's columns 4-9.
rows = list(csv.reader(open(ABC / "ANCHOR_ENTROPY_ERROR_CURVES.csv")))[1:]
cur = {}
for arm, rep in REP.items():
    sel = [r for r in rows if r[0] == rep and r[1] == "cristobalite_beta" and r[2] == "quartz_beta"]
    a = np.array([[float(x) for x in r[3:]] for r in sel])
    cur[arm] = dict(T=a[:, 0], e=a[:, 1], dS=a[:, 2], anchor=a[:, 3], integ=a[:, 4], closure=a[:, 6])
summ = {r["representation"]: r for r in csv.DictReader(open(ABC / "ANCHOR_ENTROPY_PAIR_SUMMARY.csv"))
        if r["gap_first_minus_second"] == "cristobalite_beta-quartz_beta"}
T0 = float(summ["B_DPA"]["T0_K"])

c = np.load(NPZ); T = c["T_K"]
ref = {}
for p in POOL:
    rr = list(csv.DictReader(open(REF / p / "reference_G.csv")))
    assert np.array_equal(np.array([float(r["T_K"]) for r in rr]), T)
    ref[p] = np.array([float(r["G_eV_per_atom"]) for r in rr]) * 1000
pred = {arm: {p: c[KEY.format(arm=arm, p=p)] * 1000 for p in POOL} for arm in ("B", "C")}
margin = lambda d: d["cristobalite_beta"] - np.minimum(d["quartz_beta"], d["tridymite_p63mmc"])
m_ref = margin(ref)

checks = {}
for arm in ("B", "C"):
    d = cur[arm]
    assert np.array_equal(d["T"], T)
    # 1. archived error curve = regenerated prediction - reference
    e_regen = (pred[arm]["cristobalite_beta"] - pred[arm]["quartz_beta"]) - (ref["cristobalite_beta"] - ref["quartz_beta"])
    dev = float(np.max(np.abs(e_regen - d["e"])))
    assert dev < 1e-5, (arm, dev)  # floating-point level; the regeneration matches gap MAE to ~5e-11 relative
    assert np.max(np.abs(d["closure"])) < 1e-7 and np.allclose(d["anchor"] + d["integ"], d["e"], atol=1e-7)
    # 2. archived summary values
    s = summ[REP[arm]]
    i_end = len(T) - 1
    for got, key in ((d["anchor"][0], "anchor_error_at_T0_meV_atom"), (d["e"][i_end], "at_2499K_deltaG_error_meV_atom"),
                     (d["integ"][i_end], "at_2499K_integrated_entropy_contribution_meV_atom"),
                     (np.mean(np.abs(d["e"])), "deltaG_error_MAE_meV_atom"), (np.mean(np.abs(d["dS"])), "deltaS_error_MAE_micro_eV_atom_K")):
        assert abs(got - float(s[key])) < 1e-9, (arm, key, got, s[key])
    # 3. decision rule: cristobalite predicted stable <=> predicted margin < 0 <=> e_margin < -m_ref.
    #    The plotted gap is crist - quartz; it differs from the margin error only through the
    #    (small) tridymite - quartz error, so count where the simple rule disagrees with the labels.
    m_pred = margin(pred[arm])
    stable = m_pred < 0
    rule = d["e"] < -m_ref
    mism = int(np.sum(stable != rule))
    entry = T[np.flatnonzero(stable)[0]] if stable.any() else None
    checks[arm] = dict(regen_max_abs_dev_meV=dev, rule_label_mismatches=mism, first_stable_T_K=entry,
                       max_abs_margin_error_minus_gap_error_meV=float(np.max(np.abs((m_pred - m_ref) - d["e"]))))
    assert mism <= 3, checks[arm]

# ---- figure ----
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9.5, "axes.edgecolor": "#8aa0ae", "axes.labelcolor": INK,
                     "xtick.color": INK, "ytick.color": INK})
fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.1), dpi=200, sharey=True)
YLO, YHI = -0.85, 0.75
for ax, arm in zip(axes, ("B", "C")):
    d = cur[arm]; col = COLOR[arm]; e0 = d["anchor"][0]
    ax.fill_between(T, YLO, np.maximum(-m_ref, YLO), color=GREEN, alpha=0.10, lw=0)
    ax.plot(T, -m_ref, color=GREEN, lw=1.0)
    ax.axhline(0, color="#8aa0ae", lw=0.6)
    ax.axvline(T0, color=INK, lw=0.7, ls=":")
    ax.fill_between(T, e0, d["e"], color=col, alpha=0.22, lw=0)
    ax.axhline(e0, color=INK, lw=1.0, ls="--")
    ax.plot(T, d["e"], color=col, lw=2.3)
    ax.plot([T0], [e0], "o", color=INK, ms=4, zorder=6)
    s = summ[REP[arm]]
    ax.set_title(f"{LABEL[arm]}   ·   ΔS MAE {float(s['deltaS_error_MAE_micro_eV_atom_K']):.3f} µeV/(atom·K)",
                 fontsize=9.5, loc="left", color=col, fontweight="bold", pad=6)
    ax.annotate(f"anchor error at T₀\n{e0:+.3f}", xy=(T0, e0), xytext=(T0 - 330, 0.64), fontsize=8, color=INK, ha="center",
                va="center", arrowprops=dict(arrowstyle="->", color=INK, lw=0.6))
    TE, ee, ie = T[-1], d["e"][-1], d["integ"][-1]
    if abs(ie) > 0.15:
        ax.annotate("", xy=(TE + 18, ee), xytext=(TE + 18, e0), arrowprops=dict(arrowstyle="->", color=col, lw=1.1))
        ax.text(TE - 20, (e0 + ee) / 2, f"signed entropy-\nerror integral\n{ie:+.3f}", fontsize=8, color=col, ha="right", va="center")
        ax.text(TE - 20, ee - 0.09, f"final gap error at {TE:.0f} K: {ee:+.3f}", fontsize=8, color=col, ha="right", va="center", fontweight="bold")
    else:  # integrated term too small to draw as an arrow; label above the curve
        ax.text(TE - 20, 0.70, f"final gap error at {TE:.0f} K: {ee:+.3f}", fontsize=8, color=col, ha="right", va="top", fontweight="bold")
        ax.text(TE - 20, 0.62, f"signed entropy-error integral: {ie:+.3f}", fontsize=8, color=col, ha="right", va="top")
    ent = checks[arm]["first_stable_T_K"]
    if ent is not None:
        ax.axvline(ent, color=col, lw=0.7, ls=":")
        ax.text(ent + 15, YLO + 0.05, f"cristobalite enters\nat {ent:.0f} K", fontsize=7.5, color=col, ha="left", va="bottom")
    else:
        ax.text(TE - 20, YLO + 0.05, "error never falls into the band:\ncristobalite never predicted stable", fontsize=7.5,
                color=col, ha="right", va="bottom")
    ax.text(T0 + 12, YLO + 0.05, f"T₀ = {T0:.0f} K", fontsize=7.5, color=INK, ha="left", va="bottom", rotation=90)
    ax.set_xlim(T[0], T[-1] + 40); ax.set_ylim(YLO, YHI)
    ax.set_xlabel("Temperature (K)")
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
axes[0].set_ylabel("Error in ΔG(cristobalite − quartz)\n(meV/atom)", fontsize=9)
axes[0].text(T[0] + 15, -0.15, "green band: error below\nthe reference margin,\ncristobalite still\npredicted stable",
             fontsize=7.5, color="#3d7a4a", ha="left", va="top")
fig.tight_layout(w_pad=2.0)
for ext in ("png", "svg", "pdf"):
    fig.savefig(OUT / f"ERROR_DECOMPOSITION_SIO2.{ext}", bbox_inches="tight")

with open(OUT / "ERROR_DECOMPOSITION_SOURCE_DATA.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["T_K", "ref_crst_margin_to_nearest_competitor_meV",
                *[f"{LABEL[a].replace(' ', '')}_{k}" for a in ("B", "C") for k in
                  ("gap_error_meV", "anchor_error_meV", "integrated_entropy_term_meV", "dS_error_microeV_atom_K", "crst_predicted_stable")]])
    stab = {a: margin(pred[a]) < 0 for a in ("B", "C")}
    for i in range(len(T)):
        w.writerow([T[i], m_ref[i], *[v for a in ("B", "C") for v in
                    (cur[a]["e"][i], cur[a]["anchor"][i], cur[a]["integ"][i], cur[a]["dS"][i], bool(stab[a][i]))]])
summary = {"model": "held-out cristobalite; sources quartz, P6_3/mmc tridymite, C222_1 tridymite; equal-block scaling; Si/O task checkpoint",
           "gap": "cristobalite_beta - quartz_beta", "T0_K": T0,
           "archived_values_meV_atom": {LABEL[a]: {k: float(summ[REP[a]][k]) for k in
               ("anchor_error_at_T0_meV_atom", "deltaG_error_MAE_meV_atom", "deltaS_error_MAE_micro_eV_atom_K",
                "at_2499K_deltaG_error_meV_atom", "at_2499K_integrated_entropy_contribution_meV_atom")} for a in ("B", "C")},
           "checks": checks,
           "note": "Post-fit error diagnostic based on a thermodynamic identity; not an independently established causal mechanism."}
(OUT / "ERROR_DECOMPOSITION_SUMMARY.json").write_text(json.dumps(summary, indent=1) + "\n")
print(json.dumps(summary, indent=1))
