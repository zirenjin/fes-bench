"""Carbon negative control: same-phase fixed-volume Svib(300 K), DPA readout vs quadratic V/V0 baseline.

Archived development results only (updates-20261030/carbon). DPA = arm B, coordinate-standardized,
within-phase leave-one-volume-out, pooled over the 5 interior + 2 endpoint folds (n = 7 per phase).
Quadratic = fixed OLS on [1, V/V0, (V/V0)^2], same folds. Values are asserted against the
FOLLOWUP_RESULTS.md table (0.00743 / 0.00920 / 0.00467 / 0.00372 µeV/(atom K)).
Usage: plot_carbon_control.py <carbon_dir> <out_dir>
"""
import csv
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

CARBON, OUT = Path(sys.argv[1]), Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)
INK, TEAL, GREY = "#203546", "#087f80", "#8aa0ae"
TABLE = {("graphite", "DPA"): 0.00743, ("diamond", "DPA"): 0.00920, ("graphite", "quadratic"): 0.00467, ("diamond", "quadratic"): 0.00372}

mae = {}
m = [r for r in csv.DictReader(open(CARBON / "METRICS.csv")) if r["scale_rule"] == "coordinate_standardized" and r["arm"] == "B"
     and r["experiment"] == "within_phase_leave_one_volume_out"]
for ph in ("graphite", "diamond"):
    rows = [r for r in m if r["phase"] == ph]
    n = sum(int(r["n_structures"]) for r in rows); assert n == 7
    mae[(ph, "DPA")] = sum(int(r["n_structures"]) * float(r["mae_ueV_atom_K"]) for r in rows) / n
for r in csv.DictReader(open(CARBON / "VOLUME_BASELINE_METRICS.csv")):
    if r["model"] == "quadratic_VV0" and r["support_role"] == "all":
        assert int(r["n"]) == 7
        mae[(r["phase"], "quadratic")] = float(r["mae"])
for k, v in TABLE.items():
    assert abs(mae[k] - v) < 5e-6, (k, mae[k], v)

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": GREY, "axes.labelcolor": INK,
                     "xtick.color": INK, "ytick.color": INK})
fig, ax = plt.subplots(figsize=(5.2, 2.2), dpi=200)
bars = [("graphite", "quadratic"), ("graphite", "DPA"), ("diamond", "quadratic"), ("diamond", "DPA")]
ys = [3.4, 2.6, 1.0, 0.2]
LAB = {"DPA": "DPA readout", "quadratic": "quadratic volume"}
for y, k in zip(ys, bars):
    col = GREY if k[1] == "quadratic" else TEAL
    ax.barh(y, mae[k], height=0.62, color=col)
    ax.text(mae[k] + 0.00015, y, f"{mae[k]:.5f}", va="center", ha="left", fontsize=8, color=INK)
ax.set_yticks(ys, [LAB[k[1]] for k in bars], fontsize=8.5)
for y, ph in ((3.0, "Graphite"), (0.6, "Diamond")):
    ax.text(-0.0003, y + 0.95, ph, ha="right", va="center", fontsize=8.5, fontweight="bold", color=INK, transform=ax.transData)
ax.set_xlim(0, 0.0112); ax.set_ylim(-0.45, 4.45)
ax.set_xlabel("Same-phase Svib(300 K) MAE, µeV/(atom·K)  ·  lower is better", fontsize=8.5)
for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
ax.tick_params(axis="y", length=0)
fig.tight_layout()
for ext in ("png", "svg", "pdf"):
    fig.savefig(OUT / f"CARBON_NEGATIVE_CONTROL.{ext}", bbox_inches="tight")
summary = {"task": "graphite/diamond fixed-volume harmonic Svib(300 K), within-phase leave-one-volume-out, n = 7 per phase",
           "mae_ueV_atom_K": {f"{ph}_{mdl}": v for (ph, mdl), v in mae.items()},
           "note": "Development-only negative control; does not test or confirm the SiO2 decision mechanism."}
(OUT / "CARBON_NEGATIVE_CONTROL.json").write_text(json.dumps(summary, indent=1) + "\n")
print(json.dumps(summary, indent=1))
