"""H1b external cross-tier stress test (ds62 assessed tier), page version of FIG_H1B_DG_CURVES.

Same data as the archived post-hoc figure: posthoc/REGENERATED_PRIMARY_PAIR_CURVES.npz, whose
regeneration check (byte-identical FROZEN_SELECTION.csv, gap MAE difference 0) is required to pass.
Restyled for the page without an in-image title. Before drawing, the script asserts the statements
used in the page caption: at 3 GPa no fit changes sign (24/24); at 10 GPa the three-source B readout
crosses at 1354 K (both scalings) and C (coordinate) at 1368 K, 377-391 K before the reference T0.
Post-hoc and descriptive; the ds62 tier is not DFT and this is not an independent confirmation.
Usage: plot_h1b_curves.py <h1b_run_dir> <out_dir>
"""
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

H, OUT = Path(sys.argv[1]), Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)
INK, TEAL, SLATE, GREY = "#203546", "#087f80", "#7d98ab", "#cfd8df"
NAME = {"q": "quartz", "trd": "tridymite", "crst": "cristobalite", "coe": "coesite", "stv": "stishovite"}
# From the frozen scripts/score_h1b.py (TESTS, SIGMA_MEV).
TESTS = {"primary_3GPa": dict(iso="3GPa", target="coe", anchor="q", two="two_q+trd", three="three_q+trd+crst"),
         "secondary_10GPa": dict(iso="10GPa", target="stv", anchor="coe", two="two_q+coe", three="three_q+trd+coe")}
SIGMA_MEV = {"3GPa": 0.12456301893223591, "10GPa": 1.319156164123266}

assert json.load(open(H / "posthoc/REGENERATION_CHECK.json"))["pass"]
z = np.load(H / "posthoc/REGENERATED_PRIMARY_PAIR_CURVES.npz"); T = z["T_K"]
ref = np.load(H / "reference/REFERENCE_DS62.npz"); assert np.array_equal(ref["T_K"], T)
gates = json.load(open(H / "reference/GATES_R.json"))["R2"]
curves = {tuple(k.split("__")): z[k] * 1000 for k in z.files if k != "T_K"}


def crossings(y):
    i = np.flatnonzero(np.sign(y[:-1]) != np.sign(y[1:]))
    return [float(T[j] - y[j] * (T[j + 1] - T[j]) / (y[j + 1] - y[j])) for j in i]


facts = {}
for tname, tc in TESTS.items():
    ks = [k for k in curves if k[0] == tname]
    facts[tname] = {"n_curves": len(ks), "T0_K": gates[tc["iso"]]["T0_K"],
                    "crossing": {"__".join(k[1:]): crossings(curves[k]) for k in ks if crossings(curves[k])}}
f3, f10 = facts["primary_3GPa"], facts["secondary_10GPa"]
assert f3["n_curves"] == 24 and not f3["crossing"], f3
assert f10["n_curves"] == 48
for key, want in (("three_q+trd+coe__coordinate__B", 1354), ("three_q+trd+coe__equal_block__B", 1354), ("three_q+trd+coe__coordinate__C", 1368)):
    got = f10["crossing"][key]; assert len(got) == 1 and abs(got[0] - want) < 1.0, (key, got)
early = [round(f10["T0_K"]) - round(v[0]) for v in f10["crossing"].values()]
assert sorted(set(early)) == [377, 391], early  # RESULT.md convention: difference of rounded temperatures (unrounded 377.1-391.8 K)

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9.5, "axes.edgecolor": "#8aa0ae", "axes.labelcolor": INK,
                     "xtick.color": INK, "ytick.color": INK})
fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.9), dpi=200)
for ax, (tname, tc) in zip(axes, TESTS.items()):
    iso, tgt, anc = tc["iso"], tc["target"], tc["anchor"]
    r = (ref[f"{iso}_{tgt}_G_eV_atom"] - ref[f"{iso}_{anc}_G_eV_atom"]) * 1000
    sig = SIGMA_MEV[iso]; T0 = gates[iso]["T0_K"]; uT = gates[iso]["reference_Tc_uncertainty_K_2sigma"]
    k2, k3 = (tname, tc["two"], "coordinate", "B"), (tname, tc["three"], "coordinate", "B")
    n_other = facts[tname]["n_curves"] - 2
    first = True
    for k in curves:
        if k[0] == tname and k not in (k2, k3):
            ax.plot(T, curves[k], color=GREY, lw=0.5, alpha=0.7, zorder=1, label=f"other fits ({n_other})" if first else None); first = False
    ax.axvspan(T0 - uT, T0 + uT, color=INK, alpha=0.07, lw=0, zorder=0, label="reference T₀ ± 2σ")
    ax.fill_between(T, r - 2 * sig, r + 2 * sig, color=INK, alpha=0.15, lw=0, zorder=2)
    ax.plot(T, r, color=INK, lw=2.6, zorder=4, label="ds62 reference (±2σ)")
    src = lambda c: ", ".join(NAME[p] for p in c.split("_", 1)[1].split("+"))
    # At 10 GPa the two-source fit lies ~600-770 meV/atom above the reference; keep the axis on the crossing region.
    shown = [r, curves[k3]] + ([curves[k2]] if iso == "3GPa" else [])
    off = "" if iso == "3GPa" else "; off scale above"
    ax.plot(T, curves[k2], color=SLATE, lw=1.8, zorder=5, label=f"DPA, two sources ({src(tc['two'])}{off})")
    ax.plot(T, curves[k3], color=TEAL, lw=2.2, ls="--", zorder=5, label=f"DPA, three sources ({src(tc['three'])})")
    ax.axhline(0, color="#8aa0ae", lw=0.7); ax.axvline(T0, color=INK, lw=0.8, ls=":")
    lo, hi = np.min(shown), np.max(shown)
    pad = 0.15 * (hi - lo); ax.set_ylim(lo - pad, hi + 1.6 * pad); ax.set_xlim(T[0], T[-1])
    ax.text(T0 + uT + 12, hi + 1.45 * pad, f"reference T₀ = {T0:.0f} ± {uT:.0f} K", fontsize=8, color=INK, va="top")
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel(f"G({NAME[tgt]}) − G({NAME[anc]})  (meV/atom)", fontsize=9)
    ax.set_title(f"{iso.replace('GPa', ' GPa')} · {NAME[tgt]} held out", fontsize=9.5, loc="left", color=INK, fontweight="bold")
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    ax.legend(fontsize=7.6, frameon=False, loc="upper left", bbox_to_anchor=(0.0, -0.14), ncol=2)
f10 = facts["secondary_10GPa"]; T10 = round(f10["T0_K"])
tb, tcv = f10["crossing"]["three_q+trd+coe__coordinate__B"][0], f10["crossing"]["three_q+trd+coe__coordinate__C"][0]
axes[1].annotate(f"three-source crossings\nDPA {tb:.0f} K · DPA + volume {tcv:.0f} K\n({T10 - round(tcv)}–{T10 - round(tb)} K early)", (tb, 0),
                 xytext=(tb - 460, 33), fontsize=8, color=TEAL, ha="center", va="center", arrowprops=dict(arrowstyle="->", color=TEAL, lw=0.7))
fig.tight_layout(w_pad=2.5)
for ext in ("png", "svg", "pdf"):
    fig.savefig(OUT / f"H1B_DG_CURVES_PAGE.{ext}", bbox_inches="tight")
(OUT / "H1B_DG_CURVES_PAGE_FACTS.json").write_text(json.dumps(facts, indent=1) + "\n")
print(json.dumps({t: {k: v for k, v in f.items()} for t, f in facts.items()}, indent=1))
