"""Context figures for the database pairs: error scales against the crossing budget, and example Delta G(T) curves.
Reads the manifest and stored curves only (no model inference). Examples are chosen by a fixed rule: for DPA-3.1 (S0) on
model-valid crossing pairs, the pair with the median |e(T_x)| in each class (hit within +-50 K; location error, correct
sign of the relative entropy; response error, wrong sign).
Usage: python database_pairs/make_context_figs.py
"""
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from score import Bench, GATE

HERE = Path(__file__).parent; B = Bench(); T = B.T; i1000 = B.it[1000]
INK, MUTED, GRID = "#1a1a1a", "#555555", "#dddddd"; BLUE, PINK = "#0072B2", "#CC79A7"
plt.rcParams.update({"font.size": 8, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False, "pdf.fonttype": 42})
load = lambda f: np.load(f)
F0 = load(HERE / "predictions/F0_DPA-3.1-3M_S0.npz"); BA = load(HERE.parent / "analyses/descriptor_baseline/predictions/bartel_S1.npz")
G0 = F0["E0"][:, None] + F0["F"]; GB = BA["E0"][:, None] + BA["F"]
mv0 = B.ref_valid & (F0["fmin"][B.pi] >= GATE) & (F0["fmin"][B.pj] >= GATE); mvB = B.ref_valid & (BA["fmin"][B.pi] >= GATE) & (BA["fmin"][B.pj] >= GATE)
ev0 = F0["fmin"] >= GATE
d0 = G0[B.pi] - G0[B.pj]; dB = GB[B.pi] - GB[B.pj]; dr = B.dref


def at_tc(d, k): return float(np.interp(B.tc[k], T, d[k]))


cross = np.flatnonzero(B.ref_valid & B.cross)
rows = [("Reference $|\\Delta G|$ at 1000 K", np.abs(dr[B.ref_valid, i1000]), INK),
        ("Crossing budget $\\varepsilon_{50}$", B.eps50[cross], INK),
        ("DPA-3.1 single phase, 1000 K (S0)", np.abs((F0["F"] - B.ref["F"])[ev0 & (B.ref["fmin"] >= GATE), i1000]), BLUE),
        ("DPA-3.1 pair, 1000 K (S0)", np.abs((d0 - dr)[mv0, i1000]), BLUE),
        ("DPA-3.1 pair at crossing (S0)", np.array([abs(at_tc(d0 - dr, k)) for k in cross if mv0[k]]), BLUE),
        ("Descriptor pair at crossing (S1)", np.array([abs(at_tc(dB - dr, k)) for k in cross if mvB[k]]), PINK)]
fig, ax = plt.subplots(figsize=(7.0, 2.6))
for n, (lab, v, c) in enumerate(rows):
    y = -n; v = v[np.isfinite(v) & (v > 0)]; q1, med, q3 = np.percentile(v, [25, 50, 75]); lo, hi = np.percentile(v, [5, 95])
    ax.plot([lo, hi], [y, y], color=c, lw=1, alpha=0.6); ax.plot([q1, q3], [y, y], color=c, lw=5, solid_capstyle="butt")
    ax.plot(med, y, "o", color="white", ms=4.5, mec=c, mew=1.4); ax.text(hi * 1.15, y, f"{med:.2g} (n={len(v)})", va="center", fontsize=6.8, color=INK)
ax.axvline(np.median(B.eps50[cross]), color=MUTED, lw=0.7, ls=":")
ax.set_xscale("log"); ax.set_xlim(0.01, 100); ax.set_yticks([-n for n in range(len(rows))]); ax.set_yticklabels([r[0] for r in rows], fontsize=7)
ax.set_xlabel("meV atom$^{-1}$ (bar: interquartile range; line: 5th--95th percentile; circle: median)"); ax.grid(axis="x", color=GRID, lw=0.4)
for ext in ("pdf", "png"): fig.savefig(HERE / f"figures/fig0_error_scales.{ext}", bbox_inches="tight", dpi=200)
plt.close(fig)

# example curves
cls = {}
for k in cross:
    if not mv0[k]: continue
    e = abs(at_tc(d0 - dr, k)); tm = B.crossings(d0[k]); hit = any(abs(t - B.tc[k]) <= 50 for t in tm)
    sr = np.sign(np.gradient(dr[k], T)[np.argmin(abs(T - B.tc[k]))]); sm = np.sign(np.gradient(d0[k], T)[np.argmin(abs(T - B.tc[k]))])
    cls.setdefault("hit" if hit else ("response" if sr != sm else "location"), []).append((e, k))
pick = {c: sorted(v)[len(v) // 2][1] for c, v in cls.items()}
TITLE = {"hit": "a  Crossing within $\\pm50$ K", "location": "b  Location error", "response": "c  Response error (wrong sign of $\\Delta S$)"}
fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.3), gridspec_kw=dict(wspace=0.35))
for a, c in zip(axs, ("hit", "location", "response")):
    k = pick[c]; p = B.P[k]; sel = T >= 50
    a.axhline(0, color=GRID, lw=0.8); a.axvspan(B.tc[k] - 50, B.tc[k] + 50, color=GRID, alpha=0.6, lw=0)
    a.plot(T[sel], dr[k][sel], color=INK, lw=1.6, label="Reference (PBE)"); a.plot(T[sel], d0[k][sel], color=BLUE, lw=1.4, label="DPA-3.1")
    a.set_title(TITLE[c] + f"\n{p['formula']}, $|e(T_\\times)|$ = {abs(at_tc(d0 - dr, k)):.1f} meV atom$^{{-1}}$", loc="left", fontsize=7.5)
    a.set_xlabel("T (K)"); a.grid(color=GRID, lw=0.3)
axs[0].set_ylabel("$\\Delta G$ (meV atom$^{-1}$)"); axs[0].legend(fontsize=6.5, loc="upper right")
for ext in ("pdf", "png"): fig.savefig(HERE / f"figures/fig_dG_examples.{ext}", bbox_inches="tight", dpi=200)
print("classes", {c: len(v) for c, v in cls.items()}, "picked", {c: B.P[k]["formula"] for c, k in pick.items()})
