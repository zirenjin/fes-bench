"""Stable-phase band comparison for the SiO2 held-out-cristobalite transfer experiment.

One protocol only: frozen A/B/C readouts, held-out beta-cristobalite (no target thermal labels),
candidate pool {beta-quartz, P6_3/mmc tridymite, beta-cristobalite}, 851-2499 K, PBE-D3(BJ)
reference. Curves come from the verified post-hoc regeneration of the frozen run
(fes-decision-metric-audit-20261030/curves/SIO2_ABC_TASKCKPT.npz). Every accuracy, gap MAE and
predicted event list is asserted against the frozen PRIMARY_COMPARISON.csv, ENVELOPE_METRICS.csv
and EVENTS.csv before drawing.

Figure A: three-source models (sources beta-quartz, P6_3/mmc, C222_1), Volume-only / DPA /
DPA+volume under the primary coordinate scaling, plus DPA+volume under the pre-declared
equal-block sensitivity scaling. (Volume-only and DPA are identical under both scalings.)
Figure B: two-source {beta-quartz, P6_3/mmc} vs three-source, A/B/C, coordinate scaling.
Both are scored on the same candidate pool, as archived.
Usage: plot_stable_phase_comparison.py <curves.npz> <abc_run_dir> <sio2_reference_phase_dir> <out_dir>
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
NAMES = {"quartz_beta": "β-quartz", "tridymite_p63mmc": "β-tridymite (P6₃/mmc)", "cristobalite_beta": "β-cristobalite"}
COLOR = {"quartz_beta": "#c7d5df", "tridymite_p63mmc": "#9a7bb6", "cristobalite_beta": "#5aa469"}
INK, ERR = "#203546", "#c2410c"
TWO, THREE = "low_quartz_beta+tridymite_p63mmc", "high_all_three"
ARMNAME = {"A": "Volume-only", "B": "DPA", "C": "DPA + volume"}

c = np.load(NPZ); T = c["T_K"]
ref = {}
for p in POOL:
    rows = list(csv.DictReader(open(REF / p / "reference_G.csv")))
    assert np.array_equal(np.array([float(r["T_K"]) for r in rows]), T)
    ref[p] = np.array([float(r["G_eV_per_atom"]) for r in rows])


def envelope(curves):
    lab = np.array(POOL)[np.argmin(np.column_stack([curves[p] for p in POOL]), axis=1)]
    ev = []
    for i in np.flatnonzero(lab[:-1] != lab[1:]):  # frozen scorer's interpolation rule
        a, b = str(lab[i]), str(lab[i + 1]); d0 = curves[a][i] - curves[b][i]; d1 = curves[a][i + 1] - curves[b][i + 1]
        f = abs(d0) / (abs(d0) + abs(d1)) if d0 * d1 < 0 else .5
        ev.append((a, b, float(T[i] + f * (T[i + 1] - T[i]))))
    return lab, ev


ref_lab, ref_ev = envelope(ref)
T0 = ref_ev[0][2]
prim = {(r["n_sources"], r["scaling"], r["arm"]): r for r in csv.DictReader(open(RUN / "PRIMARY_COMPARISON.csv"))}
env = {(r["configuration"], r["scaling"], r["arm"]): r for r in csv.DictReader(open(RUN / "ENVELOPE_METRICS.csv"))
       if r["target"] == "cristobalite_beta" and r["candidate_pool"] == ";".join(POOL)}
evs = {(r["configuration"], r["scaling"], r["arm"]): json.loads(r["predicted"]) for r in csv.DictReader(open(RUN / "EVENTS.csv"))
       if r["target"] == "cristobalite_beta" and r["candidate_pool"] == ";".join(POOL)}


def model(cfg, scaling, arm):
    curves = {p: c[f"cristobalite_beta__{cfg}__{scaling}__{arm}__{p}"] for p in POOL}
    lab, ev = envelope(curves)
    e = env[(cfg, scaling, arm)]
    acc = 100 * np.mean(lab == ref_lab)
    assert abs(acc - 100 * float(e["stable_accuracy"])) < 1e-9, (cfg, scaling, arm, acc)
    arch = evs[(cfg, scaling, arm)]
    assert [(a, b) for a, b, _ in ev] == [(x["from"], x["to"]) for x in arch]
    assert all(abs(t - x["T_K"]) < 1e-6 for (_, _, t), x in zip(ev, arch))
    nsrc = "3" if cfg == THREE else "2"
    gaps = [np.mean(abs((curves["cristobalite_beta"] - curves[a]) - (ref["cristobalite_beta"] - ref[a]))) * 1000 for a in ("quartz_beta", "tridymite_p63mmc")]
    p = prim[(nsrc, scaling, arm)]
    assert abs(np.mean(gaps) - float(p["gap_MAE_meV_atom"])) < 1e-6
    return {"lab": lab, "ev": ev, "acc": acc, "gap": float(np.mean(gaps)), "hits": int(e["hits"]), "missed": int(e["missed"]), "false": int(e["false"])}


def intervals(lab):
    s = np.r_[0, np.flatnonzero(lab[1:] != lab[:-1]) + 1]; e = np.r_[s[1:], len(T)]
    return [(T[a], T[b - 1], lab[a]) for a, b in zip(s, e)]


def draw(rows, fname, title, groups=None, keys=None):
    n = len(rows)
    fig = plt.figure(figsize=(10.4, 0.62 * n + 1.9), dpi=200)
    gs = fig.add_gridspec(1, 2, width_ratios=[4.0, 1.9], wspace=0.03)
    ax, tx = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
    err_rows = []
    for k, (label, m) in enumerate(rows):
        y = n - 1 - k
        lab = ref_lab if m is None else m["lab"]
        for a, b, ph in intervals(lab):
            ax.add_patch(plt.Rectangle((a, y - 0.26), b - a + 1, 0.52, color=COLOR[ph], lw=0))
        if m is not None:
            wrong = lab != ref_lab
            for a, b, ph in intervals(np.where(wrong, lab, "ok")):
                if ph == "ok":
                    continue
                ax.add_patch(plt.Rectangle((a, y - 0.37), b - a + 1, 0.07, color=ERR, lw=0))
                err_rows.append({"figure": fname, "row": (keys[k] if keys else label), "T_start_K": a, "T_end_K": b, "predicted_phase": ph,
                                 "reference_phase": ";".join(sorted(set(ref_lab[(T >= a) & (T <= b)])))})
            for _, to, t in m["ev"]:
                ax.plot([t, t], [y + 0.26, y + 0.33], color=INK, lw=0.8)
                ax.text(t, y + 0.34, f"{t:.0f} K", ha="center", va="bottom", fontsize=6.6, color=INK)
        ax.text(T[0] - 25, y, label, ha="right", va="center", fontsize=8.6, color=INK, fontweight="bold" if m is None else "normal")
        if m is None:
            for x in (0.04, 0.40, 0.76):
                tx.text(x, y, "—", va="center", fontsize=8.6, color=INK)
        else:
            tx.text(0.04, y, f"{m['gap']:.3f}", va="center", fontsize=8.6, color=INK)
            tx.text(0.40, y, f"{m['acc']:.1f}%", va="center", fontsize=8.6, color=INK, fontweight="bold")
            tx.text(0.76, y, f"{m['hits']} / {m['missed']} / {m['false']}", va="center", fontsize=8.6, color=INK)
    if groups:
        for ytop, text in groups:
            ax.text(T[0] - 25, ytop, text, ha="right", va="bottom", fontsize=8, color="#536b7b", style="italic")
    ax.axvline(T0, color=INK, lw=0.9, ls="--")
    ax.text(T0 + 12, n - 0.42, f"reference Tc {T0:.0f} K", ha="left", va="bottom", fontsize=7.6, color=INK)
    ax.set_xlim(T[0], T[-1]); ax.set_ylim(-0.6, n - 0.12)
    ax.set_yticks([]); ax.set_xlabel("Temperature (K)")
    for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
    tx.set_xlim(0, 1.15); tx.set_ylim(ax.get_ylim()); tx.axis("off")
    hy = n - 0.42
    for x, h in ((0.04, "ΔG MAE\n(meV/atom)"), (0.40, "Stable-phase\naccuracy"), (0.76, "Transitions\nhit / miss / false")):
        tx.text(x, hy, h, va="bottom", fontsize=7.6, color="#536b7b")
    handles = [plt.Rectangle((0, 0), 1, 1, color=COLOR[p]) for p in POOL] + [plt.Rectangle((0, 0), 1, 0.2, color=ERR)]
    fig.legend(handles, [NAMES[p] for p in POOL] + ["wrong stable phase"], fontsize=7.8, frameon=False, ncol=4,
               loc="lower center", bbox_to_anchor=(0.45, -0.02))
    if title:  # page panels carry letters in the legend instead of in-image titles
        fig.suptitle(title, fontsize=9.6, color=INK, x=0.02, ha="left", y=0.995)
    for ext in ("png", "svg", "pdf"):
        fig.savefig(OUT / f"{fname}.{ext}", bbox_inches="tight")
    plt.close(fig)
    return err_rows


plt.rcParams.update({"font.family": "DejaVu Sans", "axes.edgecolor": "#8aa0ae", "xtick.color": INK, "axes.labelcolor": INK})
A_rows = [("DFT reference", None)] + [(ARMNAME[a], model(THREE, "coordinate", a)) for a in "ABC"] + [("DPA + volume (equal-block)", model(THREE, "equal_block", "C"))]
B_rows = [("DFT reference", None)] + [(f"{ARMNAME[a]}", model(TWO, "coordinate", a)) for a in "ABC"] + [(f"{ARMNAME[a]} ", model(THREE, "coordinate", a)) for a in "ABC"]
errs = draw(A_rows, "FIG_A_STABLE_PHASE_MAIN",
            None)
errs += draw(B_rows, "FIG_B_SOURCE_SUPERVISION",
             None,
             groups=[(5.45, "two sources: β-quartz, P6₃/mmc"), (2.45, "three sources: + C222₁ tridymite")],
             keys=["DFT reference"] + [f"{ARMNAME[a]} (2 sources)" for a in "ABC"] + [f"{ARMNAME[a]} (3 sources)" for a in "ABC"])
with open(OUT / "STABLE_PHASE_ERROR_INTERVALS.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(errs[0])); w.writeheader(); w.writerows(errs)
with open(OUT / "STABLE_PHASE_SOURCE_DATA.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    cols = [("A", "reference", None)] + [("A", l, m) for l, m in A_rows[1:]] + [("B", l.strip() + (" (2 src)" if i < 3 else " (3 src)"), m) for i, (l, m) in enumerate(B_rows[1:])]
    w.writerow(["T_K"] + [f"{f}:{l}" for f, l, _ in cols])
    for i in range(len(T)):
        w.writerow([T[i]] + [ref_lab[i] if m is None else m["lab"][i] for _, _, m in cols])
summary = {"protocol": "held-out beta-cristobalite; no target thermal labels; pool quartz/P6_3mmc/cristobalite; 851-2499 K; PBE-D3(BJ) reference; Si/O task checkpoint",
           "reference_Tc_K": T0,
           "rows": {f"{f}:{l.strip()}" + (" (2 src)" if f == "B" and i < 3 else " (3 src)" if f == "B" else ""): {"acc_percent": m["acc"], "gap_MAE_meV_atom": m["gap"], "events": m["ev"], "hits_missed_false": [m["hits"], m["missed"], m["false"]]}
                    for f, rows_ in (("A", A_rows[1:]), ("B", B_rows[1:])) for i, (l, m) in enumerate(rows_)}}
(OUT / "STABLE_PHASE_SUMMARY.json").write_text(json.dumps(summary, indent=1, default=float) + "\n")
print(json.dumps(summary, indent=1, default=float))


# ---- Combined page figure: one reference, two- and three-source groups, both DPA+volume scalings ----
# "per-dimension" is the frozen run's "coordinate" scaling (source-only standardization of each coordinate);
# Volume-only and DPA have a single feature group, so both scalings give the same model.
LBL = [("Volume only", "coordinate", "A"), ("DPA", "coordinate", "B"),
       ("DPA + volume (per-dimension)", "coordinate", "C"), ("DPA + volume (equal-block)", "equal_block", "C")]
GROUPS = [("Two source phases · β-quartz, P6₃/mmc tridymite", TWO), ("Three source phases · + C222₁ tridymite", THREE)]
comb = {(cfg, lab): model(cfg, sc, arm) for _, cfg in GROUPS for lab, sc, arm in LBL}
for cfg in (TWO, THREE):  # identical-model claim behind the shared labels
    for arm in "AB":
        assert all(np.array_equal(c[f"cristobalite_beta__{cfg}__coordinate__{arm}__{p}"], c[f"cristobalite_beta__{cfg}__equal_block__{arm}__{p}"]) or
                   np.allclose(c[f"cristobalite_beta__{cfg}__coordinate__{arm}__{p}"], c[f"cristobalite_beta__{cfg}__equal_block__{arm}__{p}"], atol=1e-12) for p in POOL)
d2, d3 = comb[(TWO, "DPA")], comb[(THREE, "DPA")]
assert (round(d2["gap"], 3), round(d3["gap"], 3), round(d2["acc"], 2), round(d3["acc"], 2)) == (0.853, 0.315, 97.57, 82.96)
assert [(a, b) for a, b, _ in d2["ev"]] == [("quartz_beta", "cristobalite_beta")]
assert [(a, b) for a, b, _ in d3["ev"]] == [("quartz_beta", "tridymite_p63mmc"), ("tridymite_p63mmc", "cristobalite_beta")]

ypos, rows_c, heads = {}, [], []
y = 0.0
rows_c.append((y, "DFT reference (PBE-D3(BJ))", None, None)); y -= 1.55
for head, cfg in GROUPS:
    heads.append((y + 0.62, head))
    for lab, _, _ in LBL:
        rows_c.append((y, lab, comb[(cfg, lab)], cfg)); y -= 1.0
    y -= 0.55
ylo = rows_c[-1][0] - 0.62
fig = plt.figure(figsize=(11.2, 8.9), dpi=200)
gs = fig.add_gridspec(2, 2, width_ratios=[4.6, 1.45], height_ratios=[3.0, 5.4], wspace=0.015, hspace=0.12)
ax, tx = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])
mx, lx = fig.add_subplot(gs[0, 0], sharex=ax), fig.add_subplot(gs[0, 1])
HL = "#fff4d6"

# Top panel: cristobalite margin to its nearest competitor, three sources (the hero view).
TEAL, AMBER, GREENZ = "#087f80", "#a96613", "#5aa469"
crv = lambda cfg, sc, arm: {p: c[f"cristobalite_beta__{cfg}__{sc}__{arm}__{p}"] for p in POOL}
margin = lambda d: (d["cristobalite_beta"] - np.minimum(d["quartz_beta"], d["tridymite_p63mmc"])) * 1000
m_ref, m_dpa, m_dpv = margin(ref), margin(crv(THREE, "coordinate", "B")), margin(crv(THREE, "equal_block", "C"))
TB = [t for a, b, t in comb[(THREE, "DPA")]["ev"] if b == "cristobalite_beta"][0]
assert not (comb[(THREE, "DPA + volume (equal-block)")]["lab"] == "cristobalite_beta").any() and (m_dpv > 0).all()
assert abs(np.interp(TB, T, m_dpa)) < 0.02 and abs(np.interp(T0, T, m_ref)) < 0.01
mx.axhspan(-3, 0, color=GREENZ, alpha=0.07, lw=0)
mx.axhline(0, color="#8aa0ae", lw=0.8)
mx.text(T[-1] - 15, -0.45, "below 0: cristobalite is the stable phase", ha="right", va="top", fontsize=7.6, color="#3d7a4a")
mx.plot(T, m_ref, color=INK, lw=2.4, label="DFT reference (PBE-D3(BJ))", zorder=4)
mx.plot(T, m_dpa, color=TEAL, lw=1.9, ls="--", label="DPA", zorder=5)
mx.plot(T, m_dpv, color=AMBER, lw=2.1, ls=(0, (1.0, 1.3)), label="DPA + volume (equal-block)", zorder=5)
mx.axvline(T0, color=INK, lw=0.9, ls="--", zorder=1)
mx.vlines(TB, -0.8, 0.0, color=TEAL, lw=0.8, ls=":")
mx.annotate(f"reference {T0:.0f} K", xy=(T0, 0), xytext=(T0 - 200, -0.5), fontsize=7.8, color=INK, ha="center", va="center",
            arrowprops=dict(arrowstyle="->", color=INK, lw=0.7))
mx.annotate(f"DPA: {TB:.0f} K, after a\ntridymite interval", xy=(TB, 0), xytext=(TB + 110, 1.45), fontsize=7.8, color=TEAL, ha="center",
            va="center", arrowprops=dict(arrowstyle="->", color=TEAL, lw=0.7))
mx.annotate("DPA + volume (equal-block):\nnever stable", xy=(2400, float(np.interp(2400, T, m_dpv))), xytext=(2330, 0.85), fontsize=7.8,
            color=AMBER, ha="center", va="center", arrowprops=dict(arrowstyle="->", color=AMBER, lw=0.7))
mx.set_ylim(-0.8, 2.4)
mx.set_ylabel("Margin (meV/atom)", fontsize=8.6)
mx.set_title("Cristobalite margin, G(cristobalite) − min G(quartz, tridymite)", fontsize=8.6, loc="left", color=INK, pad=4)
mx.tick_params(labelbottom=False, labelsize=8)
for s in ("top", "right"): mx.spines[s].set_visible(False)
lx.axis("off")
lx.legend(*mx.get_legend_handles_labels(), title="Three source phases", title_fontsize=7.8, fontsize=7.8, frameon=False, loc="center left")
for yy, lab, m, cfg in rows_c:
    if lab == "DPA":
        for a_ in (ax, tx):
            a_.axhspan(yy - 0.47, yy + 0.47, color=HL, lw=0, zorder=0)
    labs = ref_lab if m is None else m["lab"]
    for a, b, ph in intervals(labs):
        ax.add_patch(plt.Rectangle((a, yy - 0.26), b - a + 1, 0.52, color=COLOR[ph], lw=0, zorder=2))
    if m is not None:
        wrong = labs != ref_lab
        for a, b, ph in intervals(np.where(wrong, labs, "ok")):
            if ph != "ok":
                ax.add_patch(plt.Rectangle((a, yy - 0.40), b - a + 1, 0.08, color=ERR, lw=0, zorder=2))
        for _, _, t in m["ev"]:
            ax.plot([t, t], [yy + 0.26, yy + 0.34], color=INK, lw=0.8, zorder=3)
            ax.text(t, yy + 0.35, f"{t:.0f} K", ha="center", va="bottom", fontsize=6.6, color=INK, zorder=3)
    else:
        for _, _, t in ref_ev:
            ax.text(t - 10, yy + 0.33, f"reference {t:.0f} K", ha="right", va="bottom", fontsize=6.8, color=INK, fontweight="bold")
    bold = m is None or lab == "DPA"
    ax.text(T[0] - 22, yy, lab, ha="right", va="center", fontsize=8.6, color=INK, fontweight="bold" if bold else "normal")
    if m is not None:
        fw = "bold" if lab == "DPA" else "normal"
        tx.text(0.03, yy, f"{m['gap']:.3f}", va="center", fontsize=8.6, color=INK, fontweight=fw)
        tx.text(0.36, yy, f"{m['acc']:.2f}%", va="center", fontsize=8.6, color=INK, fontweight=fw)
        tx.text(0.73, yy, f"{m['hits']} / {m['missed']} / {m['false']}", va="center", fontsize=8.6, color=INK, fontweight=fw)
for yy, text in heads:
    ax.text(T[0] - 22, yy, text, ha="right", va="center", fontsize=8, color="#536b7b", style="italic")
    ax.plot([T[0], T[-1]], [yy - 0.12] * 2, color="#dbe4ea", lw=0.7, zorder=1)
ax.axvline(T0, color=INK, lw=0.9, ls="--", zorder=1)
ax.set_xlim(T[0], T[-1]); ax.set_ylim(ylo, 0.95)
ax.set_yticks([]); ax.set_xlabel("Temperature (K)", fontsize=9)
for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
tx.set_xlim(0, 1.05); tx.set_ylim(ax.get_ylim()); tx.axis("off")
for x, h in ((0.03, "ΔG MAE\n(meV/atom)"), (0.36, "Stable-phase\naccuracy"), (0.73, "Transitions\nhit / miss / false")):
    tx.text(x, 0.55, h, va="bottom", fontsize=7.6, color="#536b7b")
handles = [plt.Rectangle((0, 0), 1, 1, color=COLOR[p]) for p in POOL] + [plt.Rectangle((0, 0), 1, 0.2, color=ERR),
           plt.Line2D([], [], color=INK, ls="--", lw=0.9), plt.Rectangle((0, 0), 1, 1, color=HL)]
fig.legend(handles, [NAMES[p] for p in POOL] + ["wrong stable phase", f"reference transition ({T0:.0f} K)", "DPA, two vs three sources"],
           fontsize=7.6, frameon=False, ncol=6, loc="lower center", bbox_to_anchor=(0.5, -0.01))
for ext in ("png", "svg", "pdf"):
    fig.savefig(OUT / f"FIG_PHASE_DECISIONS.{ext}", bbox_inches="tight")
plt.close(fig)
with open(OUT / "FIG_PHASE_DECISIONS_SOURCE_DATA.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    cols = [(f"{'2src' if cfg == TWO else '3src'}:{lab}", m) for _, lab, m, cfg in rows_c[1:]]
    w.writerow(["T_K", "reference"] + [k for k, _ in cols])
    for i in range(len(T)):
        w.writerow([T[i], ref_lab[i]] + [m["lab"][i] for _, m in cols])
comb_summary = {"reference_events": ref_ev, "scaling_note": "per-dimension = frozen-run 'coordinate' scaling",
                "rows": {f"{'2src' if cfg == TWO else '3src'}:{lab}": {"gap_MAE_meV_atom": m["gap"], "acc_percent": m["acc"], "events": m["ev"],
                                                                       "hits_missed_false": [m["hits"], m["missed"], m["false"]]} for _, lab, m, cfg in rows_c[1:]}}
(OUT / "FIG_PHASE_DECISIONS_SUMMARY.json").write_text(json.dumps(comb_summary, indent=1, default=float) + "\n")
