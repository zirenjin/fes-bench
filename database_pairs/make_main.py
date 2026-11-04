"""FES-Bench v0.1: the authoritative main table (CSV + LaTeX), the four main figures, one SI figure and the
key numbers quoted in the text, from the stored predictions only.

Subsets (score.py): accuracy columns use each prediction's model-valid pairs; failure/success columns use all
reference-valid pairs; paired comparisons between predictions (Fig. 2, Fig. 4, key_numbers) use the common-valid
pairs (valid for DFT and both models). Truncated-harmonic curves appear only in the 'all' diagnostic table.
Usage: python make_main.py  -> tables/, figures/
"""
import csv, json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from score import Bench, GATE

HERE = Path(__file__).parent; (HERE / "tables").mkdir(exist_ok=True); (HERE / "figures").mkdir(exist_ok=True)
B = Bench(); T = B.T
MAIN = [("DPA-3.1 (F0)", "F0_DPA-3.1-3M"), ("MatterSim-v1", "MatterSim-v1.0.0-5M"), ("Equal-weight average", "AVG_F0+MatterSim")]
COL = {"DPA-3.1 (F0)": "#0072B2", "MatterSim-v1": "#E69F00", "Equal-weight average": "#009E73"}
MK = {"DPA-3.1 (F0)": "o", "MatterSim-v1": "s", "Equal-weight average": "D"}
INK, MUTED, GRID = "#1a1a1a", "#555555", "#dddddd"
plt.rcParams.update({"font.size": 8, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False, "pdf.fonttype": 42})
SUBS = ("model_valid", "common_valid", "all")
S = {(lab, s, sub): B.score(HERE / f"predictions/{f}_{s}.npz", sub) for lab, f in MAIN for s in ("S0", "S1") for sub in SUBS}
ent = list(csv.DictReader(open(HERE / "manifest/ENTRIES.csv")))
hours = {"DPA-3.1 (F0)": sum(float(r["seconds_F0"]) for r in ent) / 3600, "MatterSim-v1": sum(float(r["seconds_MatterSim"]) for r in ent) / 3600}
hours["Equal-weight average"] = hours["DPA-3.1 (F0)"] + hours["MatterSim-v1"]
pct = lambda x: 100 * x

# ---------------- tables ----------------
COLS = [("failure_model_invalid", "Failure (\\%)", 0, pct), ("single_dFth_rms_1000K", "Single $\\delta F_{\\mathrm{th}}$", 1, None),
        ("pair_dFth_rms_1000K", "Pair $\\delta\\Delta F_{\\mathrm{th}}$", 1, None), ("cancellation_c_1000K", "$c$", 2, None),
        ("pair_dS_rms_1000K", "$\\delta\\Delta S$", 1, None), ("wrong_ordering_mean", "Wrong order (\\%)", 1, pct),
        ("ordering_success_mean", "Order success (\\%)", 1, pct), ("near_Tc_rms_e", "$e(T_{\\times}^{\\mathrm{ref}})$", 2, None),
        ("near_Tc_within_energy_budget50", "In $\\varepsilon_{50}$ (\\%)", 0, pct), ("crossing_hit50", "Hit $\\pm50$\\,K (\\%)", 0, pct),
        ("crossing_success50", "Cross.\\ success (\\%)", 0, pct), ("tail_p95_abs_e_1000K", "p95 $|e|$", 1, None),
        ("decided_fraction_lofo_quantile_1000K", "Decided", 2, None)]
SUCCESS_COLS = {"failure_model_invalid", "ordering_success_mean", "crossing_success50"}
THERMAL_ONLY = {"single_dFth_rms_1000K", "pair_dFth_rms_1000K", "cancellation_c_1000K"}


def val(m, k, f=None):
    x = m[k]
    if k.startswith("decided"): return x["0.05"]["coverage"], None
    sc = (lambda v: f(v)) if f else (lambda v: v)
    return sc(x["value"]), (None if "ci95" not in x else [sc(c) for c in x["ci95"]])


rows = []
for s in ("S0", "S1"):
    for sub in SUBS:
        for lab, _ in MAIN:
            m = S[(lab, s, sub)]["metrics"]; r = dict(setting=s, subset=sub, prediction=lab)
            for k, _, nd, f in COLS:
                v, ci = val(m, k, f); r[k] = round(v, 4); r[k + "_ci95"] = "" if ci is None else f"{ci[0]:.4f};{ci[1]:.4f}"
            r.update(pairs_in_subset=m["validity"]["pairs_in_subset"], crossing_pairs_in_subset=m["validity"]["crossing_pairs_in_subset"],
                     inference_hours_521=round(hours[lab], 2)); rows.append(r)
with open(HERE / "tables/main_table.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


BLOCKS = [("(a) Accuracy on valid pairs and ordering", ["single_dFth_rms_1000K", "pair_dFth_rms_1000K", "cancellation_c_1000K", "pair_dS_rms_1000K",
                                                      "wrong_ordering_mean", "ordering_success_mean"], False),
          ("(b) Crossings, failure and cost", ["failure_model_invalid", "near_Tc_rms_e", "near_Tc_within_energy_budget50", "crossing_hit50",
                                              "crossing_success50", "tail_p95_abs_e_1000K", "decided_fraction_lofo_quantile_1000K"], True)]


def tex_table(sub, path, caption, label):
    """Two stacked blocks so that every column stays legible at text width."""
    L = ["\\begin{table*}[!tp]", "\\centering\\footnotesize", f"\\caption{{{caption}}}", f"\\label{{{label}}}", "\\setlength{\\tabcolsep}{4pt}"]
    for title, keys, extra in BLOCKS:
        cols = [c for c in COLS if c[0] in keys]
        L += [f"\\textit{{{title}}}\\\\[2pt]", "\\resizebox{\\linewidth}{!}{%", "\\begin{tabular}{ll" + "c" * len(cols) + ("cc" if extra else "") + "}", "\\toprule",
              " & Prediction & " + " & ".join(h for _, h, _, _ in cols) + (" & Pairs & Cost (h)" if extra else "") + " \\\\", "\\midrule"]
        for s, stext in (("S0", "S0"), ("S1", "S1")):
            for n, (lab, _) in enumerate(MAIN):
                m = S[(lab, s, sub)]["metrics"]; cells = []
                for k, _, nd, f in cols:
                    if k in THERMAL_ONLY and s == "S1": cells.append("$\\uparrow$"); continue
                    v, _ = val(m, k, f); cells.append(f"{v:.{nd}f}")
                L.append((stext if n == 0 else "") + f" & {lab} & " + " & ".join(cells)
                         + (f" & {m['validity']['pairs_in_subset']} & {hours[lab]:.1f}" if extra else "") + " \\\\")
            L.append("\\midrule" if s == "S0" else "\\bottomrule")
        L += ["\\end{tabular}}", "\\\\[6pt]"]
    L += ["\\end{table*}"]; Path(path).write_text("\n".join(L) + "\n", encoding="utf-8")


CAP_BASE = ("S0: no target thermodynamic labels; S1: reference static energy given (both conditional on the reference geometry). "
            "Energies in meV\\,atom$^{-1}$ at 1000\\,K unless noted; $\\delta\\Delta S$ in $\\mu$eV\\,atom$^{-1}$\\,K$^{-1}$. "
            "Pairs touching a reference entry below the harmonic gate ($-0.05$\\,THz) have no trustworthy target and are not scored (31 of 351). "
            "Failure: share of the remaining 320 pairs for which the prediction has no valid harmonic curve (for the average: either model). "
            "Order success is computed over the reference-valid pairs that are non-degenerate ($|\\Delta G^{\\mathrm{ref}}|\\ge1$\\,meV\\,atom$^{-1}$) at each of 600, 1000 and 1500\\,K and averaged over the three temperatures; crossing success is computed over the 93 reference-valid crossing pairs; both count a model-invalid pair as a failure. "
            "In $\\varepsilon_{50}$: $|e(T_\\times^{\\mathrm{ref}})|\\le|\\Delta S^{\\mathrm{ref}}|\\times50$\\,K (an energy budget); "
            "Hit $\\pm50$\\,K: a model crossing within 50\\,K of the reference pairwise crossing (a temperature outcome). "
            "Single-phase and pair RMS are pooled over entries or pairs; $c$, ordering, crossing and budget fractions are formula-weighted. "
            "Decided: fraction of non-degenerate pairs decided by leave-one-formula-out empirical-quantile abstention at $\\alpha=0.05$ (not a coverage guarantee). "
            "Cost: summed stored per-entry wall time (forces plus phonopy; approximate). $\\uparrow$: unchanged by the static control. "
            "95\\% formula-cluster intervals are in the released table. Development evidence on previously analysed data.")
tex_table("model_valid", HERE / "tables/main_table.tex",
          "FES-Bench main table: fixed-geometry PBE harmonic task (239 formulas; 351 pairs, 320 reference-valid; 93 reference-valid pairwise crossings). "
          "Accuracy columns use each prediction's model-valid pairs (``Pairs''). " + CAP_BASE, "tab:main")
tex_table("common_valid", HERE / "tables/main_table_common_valid.tex",
          "As Table~\\ref{tab:main}, with accuracy columns restricted to the 219 pairs valid for DFT and both models (used for all paired comparisons). " + CAP_BASE,
          "tab:main-valid")
tex_table("all", HERE / "tables/main_table_all_diagnostic.tex",
          "Diagnostic only: accuracy columns over all 351 pairs, including truncated-harmonic curves of entries below the gate and reference-invalid pairs. "
          "These curves are not valid harmonic free energies. " + CAP_BASE, "tab:main-all")

# ---------------- key numbers for the text ----------------
K = {}
rv_c = B.ref_valid & B.cross; eps = B.eps50[rv_c]
K["budget"] = dict(n_crossing_pairs_reference_valid=int(rv_c.sum()), n_formulas=int(len(set(B.pf[rv_c]))),
                   eps50_median=float(np.median(eps)), eps50_iqr=[float(np.percentile(eps, 25)), float(np.percentile(eps, 75))],
                   abs_dS_median=float(np.median(eps) / 0.05), K_per_meV_at_median=float(50 / np.median(eps)),
                   on_envelope=int((rv_c & B.envelope).sum()))
for lab, f in MAIN:
    z = np.load(HERE / f"predictions/{f}_S0.npz"); G = z["E0"][:, None] + z["F"]; d = G[B.pi] - G[B.pj]
    mv = B.ref_valid & (z["fmin"][B.pi] >= GATE) & (z["fmin"][B.pj] >= GATE); c = np.flatnonzero(mv & B.cross)
    e = np.array([np.interp(B.tc[k], T, d[k]) for k in c]); P = []
    for k in c:
        pc = B.crossings(d[k])
        if pc:
            tm = min(pc, key=lambda q: abs(q - B.tc[k])); Sm = -np.interp(B.tc[k], T, np.gradient(d[k], T)); Sr = -np.interp(B.tc[k], T, np.gradient(B.dref[k], T))
            P.append((tm - B.tc[k], np.interp(B.tc[k], T, d[k]) / Sm, np.interp(B.tc[k], T, d[k]) / Sr))
    P = np.array(P)
    K[lab] = dict(n_valid_crossing=len(c), median_abs_e=float(np.median(np.abs(e))), median_e_over_eps50=float(np.median(np.abs(e) / B.eps50[c])),
                  frac_with_model_crossing=float(len(P) / len(c)), median_abs_dTc=float(np.median(np.abs(P[:, 0]))),
                  frac_model_slope_within50=float(np.mean(np.abs(P[:, 0] - P[:, 1]) <= 50)), frac_ref_slope_within50=float(np.mean(np.abs(P[:, 0] - P[:, 2]) <= 50)))
for lab, _ in MAIN:
    for s in ("S0", "S1"):
        m = S[(lab, s, "model_valid")]["metrics"]
        K[f"{lab}|{s}|model_valid"] = {k: (val(m, k, f)[0], val(m, k, f)[1]) for k, _, _, f in COLS}


def paired(a, b, k, scale=1.0, ratio=False):
    ma, mb = a["metrics"][k], b["metrics"][k]
    if ratio: v, dr = ma["value"] / mb["value"], ma["_draws"] / mb["_draws"]
    else: v, dr = (ma["value"] - mb["value"]) * scale, (ma["_draws"] - mb["_draws"]) * scale
    lo, hi = np.nanpercentile(dr, [2.5, 97.5]); return dict(value=float(v), ci95=[float(lo), float(hi)])


cF, cM, cA = (S[(lab, "S0", "common_valid")] for lab, _ in MAIN)
K["MS_vs_F0_common_valid"] = dict(single_ratio=paired(cM, cF, "single_dFth_rms_1000K", ratio=True), pair_ratio=paired(cM, cF, "pair_dFth_rms_1000K", ratio=True),
                                  pair_ratio_fw=paired(cM, cF, "pair_dFth_rms_1000K_fw", ratio=True), dc=paired(cM, cF, "cancellation_c_1000K"),
                                  dS_diff=paired(cM, cF, "pair_dS_rms_1000K"), wrong_mean_diff_points=paired(cM, cF, "wrong_ordering_mean", 100),
                                  near_Tc_diff=paired(cM, cF, "near_Tc_rms_e"))
K["AVG_vs_F0_common_valid"] = {k: paired(cA, cF, k, 100 if k == "wrong_ordering_mean" else 1) for k in ("near_Tc_rms_e", "wrong_ordering_mean", "pair_dS_rms_1000K", "pair_dG_rms_1000K")}
K["AVG_vs_MS_common_valid"] = {k: paired(cA, cM, k, 100 if k == "wrong_ordering_mean" else 1) for k in ("near_Tc_rms_e", "wrong_ordering_mean", "pair_dS_rms_1000K", "pair_dG_rms_1000K")}
sF, sM, sA = (S[(lab, "S0", "model_valid")] for lab, _ in MAIN)
K["ordering_success_diff_points"] = dict(AVG_vs_F0=paired(sA, sF, "ordering_success_mean", 100), AVG_vs_MS=paired(sA, sM, "ordering_success_mean", 100))

# ---------------- Figure 1: precision requirement ----------------
fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.7), gridspec_kw=dict(wspace=0.35))
x = np.sort(eps); ax[0].step(x, np.arange(1, len(x) + 1) / len(x), where="post", color=INK, lw=1.6, ls="--", label=f"Budget $\\varepsilon_{{50}}$ (n={len(x)})")
for lab, f in MAIN:
    z = np.load(HERE / f"predictions/{f}_S0.npz"); G = z["E0"][:, None] + z["F"]; d = G[B.pi] - G[B.pj]
    mv = B.ref_valid & (z["fmin"][B.pi] >= GATE) & (z["fmin"][B.pj] >= GATE); c = np.flatnonzero(mv & B.cross)
    e = np.sort(np.abs([np.interp(B.tc[k], T, d[k]) for k in c]))
    ax[0].step(e, np.arange(1, len(e) + 1) / len(e), where="post", color=COL[lab], lw=1.4,
               label=f"$|e|$, {lab.replace(' (F0)', '').replace('Equal-weight average', 'average')} (n={len(c)})")
ax[0].set_xscale("log"); ax[0].set_xlabel("meV atom$^{-1}$ (reference-valid crossing pairs)"); ax[0].set_ylabel("Cumulative fraction of pairs")
ax[0].axvline(np.median(eps), color=MUTED, lw=0.6, ls=":"); ax[0].text(np.median(eps) * 1.08, 0.04, f"median {np.median(eps):.2f}", color=MUTED, fontsize=7)
ax[0].legend(fontsize=6.3, loc="upper left"); ax[0].grid(axis="x", color=GRID, lw=0.4); ax[0].set_title("a  Local accuracy needed vs achieved", loc="left", fontsize=8.5)
for lab, f in MAIN[:2]:
    z = np.load(HERE / f"predictions/{f}_S0.npz"); G = z["E0"][:, None] + z["F"]; d = G[B.pi] - G[B.pj]; P = []
    mv = B.ref_valid & (z["fmin"][B.pi] >= GATE) & (z["fmin"][B.pj] >= GATE)
    for k in np.flatnonzero(mv & B.cross):
        pc = B.crossings(d[k])
        if not pc: continue
        tm = min(pc, key=lambda q: abs(q - B.tc[k])); Sm = -np.interp(B.tc[k], T, np.gradient(d[k], T))
        if Sm != 0: P.append((np.interp(B.tc[k], T, d[k]) / Sm, tm - B.tc[k]))
    P = np.array(P); ax[1].scatter(P[:, 0], P[:, 1], s=12, marker=MK[lab], facecolor="none", edgecolor=COL[lab], lw=0.9, label=f"{lab} (n={len(P)})")
lim = 1200; ax[1].plot([-lim, lim], [-lim, lim], color=MUTED, lw=0.7); ax[1].fill_between([-lim, lim], [-lim - 50, lim - 50], [-lim + 50, lim + 50], color=GRID, alpha=0.6, lw=0)
ax[1].set_xlim(-lim, lim); ax[1].set_ylim(-lim, lim); ax[1].set_xlabel("$e(T_\\times^{\\mathrm{ref}})/\\Delta S^{\\mathrm{model}}$ (K)"); ax[1].set_ylabel("Observed $\\delta T_\\times$ (K)")
ax[1].legend(fontsize=6.5, loc="upper left"); ax[1].set_title("b  Consistency with the linearised root shift", loc="left", fontsize=8.5)
fig.savefig(HERE / "figures/fig1_precision_requirement.pdf", bbox_inches="tight"); fig.savefig(HERE / "figures/fig1_precision_requirement.png", dpi=200, bbox_inches="tight"); plt.close(fig)

# ---------------- Figure 2: capability profile (common-valid pairs) ----------------
PAN = [("single_dFth_rms_1000K", "Single-phase $\\delta F_{\\mathrm{th}}$ RMS\n(meV atom$^{-1}$)", 1), ("pair_dFth_rms_1000K", "Pair $\\delta\\Delta F_{\\mathrm{th}}$ RMS\n(meV atom$^{-1}$)", 1),
       ("cancellation_c_1000K", "Common-mode\nfraction $c$", 1), ("pair_dS_rms_1000K", "Pair $\\delta\\Delta S$ RMS\n($\\mu$eV atom$^{-1}$ K$^{-1}$)", 1),
       ("wrong_ordering_mean", "Wrong ordering,\nvalid pairs (%)", 100), ("ordering_success_mean", "Ordering success,\nall pairs (%)", 100)]
fig, axs = plt.subplots(1, len(PAN), figsize=(7.6, 2.0), gridspec_kw=dict(wspace=1.0))
for a, (k, t, sc) in zip(axs, PAN):
    for n, (lab, _) in enumerate(MAIN):
        m = S[(lab, "S0", "common_valid")]["metrics"][k]; v = m["value"] * sc; lo, hi = [q * sc for q in m["ci95"]]
        a.errorbar(v, -n, xerr=[[v - lo], [hi - v]], fmt=MK[lab], color=COL[lab], ms=4.5, capsize=2, lw=1)
    a.set_yticks([0, -1, -2]); a.set_yticklabels([lab for lab, _ in MAIN] if a is axs[0] else [], fontsize=7); a.set_ylim(-2.6, 0.6)
    a.set_title(t, fontsize=7); a.grid(axis="x", color=GRID, lw=0.4); a.tick_params(axis="x", labelsize=6.5); a.xaxis.set_major_locator(plt.MaxNLocator(3))
fig.savefig(HERE / "figures/fig2_capability_profile.pdf", bbox_inches="tight"); fig.savefig(HERE / "figures/fig2_capability_profile.png", dpi=200, bbox_inches="tight"); plt.close(fig)

# ---------------- Figure 3: workflow validity and decided fraction ----------------
fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.6), gridspec_kw=dict(wspace=0.4, width_ratios=[1, 1.2]))
nE = len(ent); cnt = [("DFT reference", sum(float(r["fmin_dft_THz"]) < GATE for r in ent), INK), ("DPA-3.1 (F0)", sum(float(r["fmin_F0_THz"]) < GATE for r in ent), COL["DPA-3.1 (F0)"]),
                      ("MatterSim-v1", sum(float(r["fmin_MatterSim_THz"]) < GATE for r in ent), COL["MatterSim-v1"])]
rel = {k: sum(r["relax_pipeline_F0"] == k for r in ent) for k in ("unstable_at_mlip_geometry", "spacegroup_changed", "missing_or_failed")}
labels = [c[0] + "\n(DFT geometry)" for c in cnt] + ["DPA-3.1 (F0)\nrelaxed pipeline"]
for n, c in enumerate(cnt):
    v = c[1] / nE * 100; ax[0].barh(-n, v, color=c[2], height=0.6); ax[0].text(v + 0.3, -n, f"{c[1]}", va="center", fontsize=7, color=INK)
left = 0
for k, col in (("unstable_at_mlip_geometry", COL["DPA-3.1 (F0)"]), ("spacegroup_changed", "#7fb5d8"), ("missing_or_failed", "#bcd8ea")):
    w = rel[k] / nE * 100; ax[0].barh(-3, w, left=left, color=col, height=0.6, hatch="////" if k == "unstable_at_mlip_geometry" else None, edgecolor="white", lw=0.3); left += w
ax[0].text(left + 0.3, -3, f"{sum(rel.values())}: 44 unstable,\n6 symmetry, 1 failed", va="center", fontsize=6.3, color=INK)
ax[0].set_yticks([0, -1, -2, -3]); ax[0].set_yticklabels(labels, fontsize=6.5); ax[0].set_xlabel("Entries without a valid harmonic $F$ (% of 521)")
ax[0].set_xlim(0, 25); ax[0].set_title("a  Harmonic validity and pipeline failures", loc="left", fontsize=8.5)
alphas = [0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.5]; it = B.it[1000]
for lab, f in MAIN:
    z = np.load(HERE / f"predictions/{f}_S0.npz"); G = z["E0"][:, None] + z["F"]; d = (G[B.pi] - G[B.pj])[:, it]
    r = B.lofo(d, B.dref[:, it], B.common_valid, alphas)
    cv = [r[str(a)]["coverage"] for a in alphas]; rk = [100 * r[str(a)]["selective_risk"] for a in alphas]
    ax[1].plot(cv, rk, marker=MK[lab], ms=3.5, color=COL[lab], lw=1.2, label=lab)
    k5 = alphas.index(0.05); ax[1].plot(cv[k5], rk[k5], marker=MK[lab], ms=7, mfc="none", mec=COL[lab], mew=1.2)
ax[1].set_xlabel("Decided fraction of non-degenerate pairs"); ax[1].set_ylabel("Wrong among decided (%)")
ax[1].set_xlim(0, 1); ax[1].set_ylim(bottom=0); ax[1].grid(color=GRID, lw=0.4); ax[1].legend(fontsize=6.5, loc="upper left")
ax[1].set_title("b  Abstention at 1000 K (LOFO empirical quantile; ring: $\\alpha=0.05$)", loc="left", fontsize=8.5)
fig.savefig(HERE / "figures/fig3_pipeline_coverage.pdf", bbox_inches="tight"); fig.savefig(HERE / "figures/fig3_pipeline_coverage.png", dpi=200, bbox_inches="tight"); plt.close(fig)

# ---------------- Figure 4 (main: A and average) and Figure S1 (all routes), common-valid pairs ----------------
base = {M: S[(lab, "S0", "common_valid")] for M, lab in (("F0", "DPA-3.1 (F0)"), ("MS", "MatterSim-v1"))}
MET = [("near_Tc_rms_e", "$\\Delta$ $e(T_\\times^{\\mathrm{ref}})$ RMS\n(meV atom$^{-1}$)", 1), ("wrong_ordering_mean", "$\\Delta$ wrong ordering\n(points)", 100),
       ("pair_dS_rms_1000K", "$\\Delta$ pair $\\delta\\Delta S$ RMS\n($\\mu$eV atom$^{-1}$ K$^{-1}$)", 1), ("pair_dG_rms_1000K", "$\\Delta$ pair $\\Delta G$ RMS, 1000 K\n(meV atom$^{-1}$)", 1)]
ALL = [("A: location", "routeA_F0_S0", "F0"), ("B: + uniform scale", "routeB_F0_S0", "F0"), ("C: spectral", "routeC_F0_S0", "F0"),
       ("A: location", "routeA_MS_S0", "MS"), ("B: + uniform scale", "routeB_MS_S0", "MS"), ("C: spectral", "routeC_MS_S0", "MS"),
       ("Equal-weight average", "AVG_F0+MatterSim_S0", "F0"), ("Equal-weight average", "AVG_F0+MatterSim_S0", "MS")]
scored = {f: B.score(HERE / f"predictions/{f}.npz", "common_valid") for _, f, _ in ALL}
diffs = {}


def forest(items, path):
    fig, axs = plt.subplots(1, len(MET), figsize=(7.2, 0.5 + 0.3 * len(items)), gridspec_kw=dict(wspace=0.25))
    for a, (k, t, s) in zip(axs, MET):
        for n, (lab, f, M) in enumerate(items):
            R = scored[f]; dr = (R["metrics"][k]["_draws"] - base[M]["metrics"][k]["_draws"]) * s; v = (R["metrics"][k]["value"] - base[M]["metrics"][k]["value"]) * s
            lo, hi = np.nanpercentile(dr, [2.5, 97.5]); diffs[f"{lab} vs {M}|{k}"] = dict(value=float(v), ci95=[float(lo), float(hi)])
            col = COL["Equal-weight average"] if lab.startswith("Equal") else (COL["DPA-3.1 (F0)"] if M == "F0" else COL["MatterSim-v1"])
            sig = hi < 0 or lo > 0
            a.errorbar(v, -n, xerr=[[v - lo], [hi - v]], fmt="o", color=col, mfc=col if sig else "white", ms=4, capsize=2, lw=1)
        a.axvline(0, color=MUTED, lw=0.7); a.set_yticks(-np.arange(len(items)))
        a.set_yticklabels([f"{lab} (vs {'F0' if M == 'F0' else 'MatterSim'})" if lab.startswith("Equal") else f"{lab} [{'F0' if M == 'F0' else 'MatterSim'}]" for lab, _, M in items] if a is axs[0] else [], fontsize=6.5)
        a.set_title(t, fontsize=7); a.grid(axis="x", color=GRID, lw=0.4); a.tick_params(axis="x", labelsize=6.5); a.set_ylim(-len(items) + 0.4, 0.6)
    fig.text(0.5, -0.06 if len(items) < 6 else -0.03, "Corrected minus uncorrected model, S0, 219 common-valid pairs; outer out-of-fold predictions; 95% paired formula-cluster "
             "bootstrap, unadjusted for multiple comparisons. Filled: interval excludes 0. Negative is better.", ha="center", fontsize=6.3, color=MUTED)
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight"); fig.savefig(path.with_suffix(".png"), dpi=200, bbox_inches="tight"); plt.close(fig)


forest([x for x in ALL if x[0].startswith(("A:", "Equal"))], HERE / "figures/fig4_intervention_boundary")
forest(ALL, HERE / "figures/figS1_all_routes")
(HERE / "tables/intervention_differences_S0_common_valid.json").write_text(json.dumps(diffs, indent=1) + "\n")
K["interventions_common_valid"] = diffs
K["intervention_events_common_valid"] = {f: scored[f]["metrics"]["events"] for _, f, _ in ALL} | {lab: base[M]["metrics"]["events"] for M, lab in (("F0", "F0 none"), ("MS", "MS none"))}
K["intervention_tail_common_valid"] = {f: (scored[f]["metrics"]["tail_p95_abs_e_1000K"]["value"], scored[f]["metrics"]["tail_n_pairs_max_abs_e_gt20"]["value"]) for _, f, _ in ALL}
(HERE / "tables/key_numbers.json").write_text(json.dumps(K, indent=1, default=float) + "\n")
print(json.dumps({k: v for k, v in K.items() if not k.endswith("|model_valid")}, indent=1, default=float))
