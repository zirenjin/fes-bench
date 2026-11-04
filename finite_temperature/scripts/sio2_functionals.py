"""Reference-level sensitivity of the SiO2 crossings: pairwise crossing temperatures under the eight exchange-correlation
treatments released in DaRUS-4999 (reference data only; no model). Writes results/SIO2_FUNCTIONALS.json and a figure.
Usage: python finite_temperature/scripts/sio2_functionals.py   (from the repository root)
"""
import csv, json
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).parents[2]; SRC = ROOT / "data/processed/sio2/functionals"; OUT = ROOT / "finite_temperature"
FUN = [("lda", "LDA"), ("pbe", "PBE"), ("r2scan", "r2SCAN"), ("pbe_d3bj", "PBE-D3(BJ)"), ("dfcx", "DFCX"),
       ("r2scan_d4", "r2SCAN-D4"), ("r2scan_d3bj", "r2SCAN-D3(BJ)"), ("r2scan_rvv10", "r2SCAN-rVV10")]
PH = ["quartz_beta", "cristobalite_beta", "tridymite_beta"]; PAIRS = [(0, 1, "quartz–cristobalite"), (0, 2, "quartz–tridymite"), (1, 2, "cristobalite–tridymite")]


def load(phase, f):
    r = list(csv.DictReader(open(SRC / f"{phase}_{f}.csv"))); return np.array([float(x["T_K"]) for x in r]), np.array([float(x["G_eV_per_atom"]) for x in r]) * 1000


res = {}
for f, name in FUN:
    G = {p: load(p, f) for p in PH}; T = G[PH[0]][0]; blk = dict(T_range=[float(T[0]), float(T[-1])])
    for a, b, lab in PAIRS:
        d = np.interp(T, *G[PH[a]]) - np.interp(T, *G[PH[b]]); s = np.sign(d); k = np.flatnonzero(s[:-1] * s[1:] < 0)
        blk[lab] = [float(T[i] - d[i] * (T[i + 1] - T[i]) / (d[i + 1] - d[i])) for i in k]
    res[name] = blk
(OUT / "results/SIO2_FUNCTIONALS.json").write_text(json.dumps(res, indent=1) + "\n")

INK, MUTED, GRID = "#1a1a1a", "#555555", "#dddddd"
plt.rcParams.update({"font.size": 8, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False, "pdf.fonttype": 42})
fig, ax = plt.subplots(figsize=(5.2, 2.6)); MK = {"quartz–cristobalite": ("o", "#0072B2"), "quartz–tridymite": ("s", "#E69F00")}
for n, (_, name) in enumerate(FUN):
    y = -n
    for lab, (mk, c) in MK.items():
        v = res[name][lab]
        if v: ax.plot(v, [y] * len(v), mk, color=c, ms=5, label=lab if n == 3 else None)
    if not res[name]["quartz–cristobalite"] and not res[name]["quartz–tridymite"]:
        ax.text(2580, y, "no crossing in range", va="center", ha="right", fontsize=6.8, color=MUTED)
ax.axvspan(1543 - 50, 1543 + 50, color=GRID, alpha=0.7, lw=0); ax.text(1543, 0.65, "$\\pm50$ K around PBE-D3(BJ)", ha="center", fontsize=6.5, color=MUTED)
ax.set_yticks([-n for n in range(len(FUN))]); ax.set_yticklabels([n for _, n in FUN], fontsize=7); ax.set_ylim(-len(FUN) + 0.4, 1.0)
ax.set_xlim(850, 2600); ax.set_xlabel("Reference crossing temperature (K)"); ax.grid(axis="x", color=GRID, lw=0.4)
ax.legend(fontsize=6.5, loc="upper left", bbox_to_anchor=(0.0, 0.92))
for ext in ("pdf", "png"): fig.savefig(OUT / f"figures/figS_sio2_functionals.{ext}", bbox_inches="tight", dpi=200)
print(json.dumps(res, indent=0))
