"""Figure 1 — why relative free energies are hard to learn.

(a) Reference ΔG = G(high-T phase) − G(low-T phase) against T − T_c: every
    crossing is aligned at the origin, so the size of ΔG near a transition is
    directly visible against a ±1 meV/atom band.
(b) T_c shift caused by a 1 meV/atom ΔG error, 1/|dΔG/dT| from a linear fit of
    the reference ΔG(T) over T_c ± FIT_HALF_WIDTH_K.
(c) Imaginary-mode fraction of the low-T (open) and high-T (filled) phase of
    each pair on a symlog axis, with the 1% QH-reliability threshold.

Rows of (b) and (c) are the same phase pairs in the same order; the coloured
row labels of (b) also serve as the legend for (a).

Reads only data/processed/ and result/tables/phase_inventory.csv.
Usage: python src/figures/figure1_difficulty.py --repo-root . --output result/figures/figure1
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

FIT_HALF_WIDTH_K = 25.0
X_WINDOW_K = 600.0
Y_LIMIT_MEV = 10.0
QH_THRESHOLD_PCT = 1.0

SYSTEM_LABEL = {"hf": "Hf", "ti": "Ti", "zr": "Zr", "sio2": "SiO$_2$"}
PHASE_SHORT = {"hcp": "hcp", "bcc": "bcc", "quartz_beta": "β-qz",
               "cristobalite_beta": "β-crs", "tridymite_p63mmc": "β-trd"}
# Okabe–Ito, keyed by (system, frozenset of phases)
PAIR_COLOR = {
    ("sio2", frozenset({"quartz_beta", "cristobalite_beta"})): "#E69F00",
    ("sio2", frozenset({"quartz_beta", "tridymite_p63mmc"})): "#D55E00",
    ("hf", frozenset({"hcp", "bcc"})): "#0072B2",
    ("ti", frozenset({"hcp", "bcc"})): "#009E73",
    ("zr", frozenset({"hcp", "bcc"})): "#CC79A7",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_reference(root: Path, system: str, phase: str) -> pd.DataFrame:
    df = pd.read_csv(root / "data/processed" / system / phase / "reference_G.csv")
    return df[["T_K", "G_eV_per_atom"]]


def find_pairs(root: Path) -> list[dict]:
    """Every phase pair with exactly one reference crossing, oriented high-T minus low-T."""
    pairs = []
    for sysdir in sorted((root / "data/processed").iterdir()):
        sj = sysdir / "system.json"
        if not sj.exists():
            continue
        phases = json.loads(sj.read_text())["phases"]
        for a, b in itertools.combinations(phases, 2):
            m = load_reference(root, sysdir.name, a).merge(
                load_reference(root, sysdir.name, b), on="T_K", suffixes=("_a", "_b"))
            T = m.T_K.to_numpy(float)
            d = (m.G_eV_per_atom_a - m.G_eV_per_atom_b).to_numpy(float) * 1000.0  # meV, a − b
            nz = np.flatnonzero(d != 0)
            roots = []
            for i, j in zip(nz[:-1], nz[1:]):
                if np.sign(d[i]) != np.sign(d[j]):
                    roots.append(T[i] - d[i] * (T[j] - T[i]) / (d[j] - d[i]))
            if len(roots) != 1:
                continue  # no crossing (near-degenerate) or multiple crossings
            tc = roots[0]
            # low-T phase = lower G below the crossing
            below = d[T < tc]
            low, high = (a, b) if np.mean(below) < 0 else (b, a)
            dG = d if high == a else -d                       # G(high) − G(low), meV
            w = (T >= tc - FIT_HALF_WIDTH_K) & (T <= tc + FIT_HALF_WIDTH_K)
            coef, cov = np.polyfit(T[w], dG[w], 1, cov=True)
            slope, slope_err = coef[0], float(np.sqrt(cov[0, 0]))
            resid = dG[w] - np.polyval(coef, T[w])
            r2 = 1 - np.sum(resid**2) / np.sum((dG[w] - dG[w].mean())**2)
            key = (sysdir.name, frozenset({a, b}))
            pairs.append(dict(
                system=sysdir.name, low=low, high=high, Tc=tc, T=T, dG=dG,
                slope=slope, slope_err=slope_err, r2=r2, n_fit=int(w.sum()),
                K_per_meV=1.0 / abs(slope), K_per_meV_err=slope_err / slope**2,
                color=PAIR_COLOR.get(key, "#555555"),
                label=f"{SYSTEM_LABEL[sysdir.name]}  {PHASE_SHORT[low]} / {PHASE_SHORT[high]}",
            ))
    return sorted(pairs, key=lambda p: -p["K_per_meV"])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo-root", default=".")
    ap.add_argument("--output", default="result/figures/figure1")
    args = ap.parse_args()
    root, out = Path(args.repo_root).resolve(), Path(args.output)
    out = out if out.is_absolute() else root / out
    out.mkdir(parents=True, exist_ok=True)

    pairs = find_pairs(root)
    inv = pd.read_csv(root / "result/tables/phase_inventory.csv", dtype=str, keep_default_na=False)
    imag = {(r.system, r.phase): 100.0 * float(r.imaginary_fraction) for r in inv.itertuples()}
    minfreq = {(r.system, r.phase): r.minimum_frequency_THz for r in inv.itertuples()}
    for p in pairs:
        if p["r2"] < 0.95:
            raise SystemExit(f"slope fit r² = {p['r2']:.3f} < 0.95 for {p['label']}")

    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8, "axes.labelsize": 8,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.linewidth": 0.8,
        "axes.spines.top": False, "axes.spines.right": False,
        "text.color": "black", "axes.labelcolor": "black",
        "xtick.color": "black", "ytick.color": "black",
    })
    fig = plt.figure(figsize=(7.0, 2.5))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.3, 1.0, 1.0], wspace=0.08,
                          left=0.085, right=0.985, bottom=0.2, top=0.9)
    ax_a = fig.add_subplot(gs[0])
    ax_b = fig.add_subplot(gs[1])
    ax_c = fig.add_subplot(gs[2], sharey=ax_b)
    # make room for (b)'s row labels between (a) and (b)
    pos_a = ax_a.get_position(); ax_a.set_position([pos_a.x0, pos_a.y0, pos_a.width * 0.86, pos_a.height])
    pos_b = ax_b.get_position(); ax_b.set_position([pos_b.x0 + 0.075, pos_b.y0, pos_b.width - 0.06, pos_b.height])
    pos_c = ax_c.get_position(); ax_c.set_position([pos_c.x0 + 0.03, pos_c.y0, pos_c.width - 0.03, pos_c.height])

    # ---------- (a) ----------
    ax_a.axhspan(-1, 1, color="0.88", zorder=0, lw=0)
    ax_a.axhline(0, color="0.6", lw=0.5, zorder=1)
    ax_a.axvline(0, color="0.6", lw=0.5, zorder=1)
    rows_a = []
    for p in pairs:
        x = p["T"] - p["Tc"]
        m = np.abs(x) <= X_WINDOW_K
        ls = (0, (3, 1.5)) if "trd" in p["label"] else "-"   # the two SiO2 pairs nearly coincide
        ax_a.plot(x[m], p["dG"][m], color=p["color"], lw=1.5, ls=ls, zorder=3, solid_capstyle="round")
        rows_a += [dict(pair=p["label"], T_minus_Tc_K=xi, delta_G_meV_per_atom=yi)
                   for xi, yi in zip(x[m], p["dG"][m])]
    ax_a.plot(0, 0, "o", color="black", ms=3, zorder=4)
    ax_a.set_xlim(-X_WINDOW_K, X_WINDOW_K)
    ax_a.set_ylim(-Y_LIMIT_MEV, Y_LIMIT_MEV)
    ax_a.set_yticks([-10, -5, 0, 5, 10])
    ax_a.set_xlabel("T − T$_c$ (K)")
    ax_a.set_ylabel("ΔG (meV/atom)")
    ax_a.text(X_WINDOW_K * 0.97, 1.3, "±1 meV/atom", ha="right", va="bottom", fontsize=6.5, color="0.35")
    ax_a.text(X_WINDOW_K * 0.97, Y_LIMIT_MEV * 0.93, "low-T phase stable", ha="right", va="top", fontsize=6.5)
    ax_a.text(-X_WINDOW_K * 0.97, -Y_LIMIT_MEV * 0.93, "high-T phase stable", ha="left", va="bottom", fontsize=6.5)

    # ---------- (b) ----------
    y = np.arange(len(pairs))[::-1]
    vals = [p["K_per_meV"] for p in pairs]
    # fit uncertainties are < 2% of the bar length; kept in panel_b.csv, not drawn
    ax_b.barh(y, vals, height=0.6, color=[p["color"] for p in pairs])
    for yi, v in zip(y, vals):
        ax_b.text(v + max(vals) * 0.03, yi, f"{v:.0f} K", va="center", fontsize=7)
    ax_b.set_yticks(y)
    ax_b.set_yticklabels([p["label"] for p in pairs])
    for tick, p in zip(ax_b.get_yticklabels(), pairs):
        tick.set_color(p["color"]); tick.set_fontweight("bold")
    ax_b.tick_params(axis="y", length=0)
    ax_b.set_xlim(0, max(vals) * 1.3)
    ax_b.set_xlabel("T$_c$ shift per 1 meV/atom (K)")
    ax_b.set_ylim(-0.6, len(pairs) - 0.4)

    # ---------- (c) ----------
    ax_c.axvspan(QH_THRESHOLD_PCT, 100, color="#fbe3e3", zorder=0, lw=0)
    ax_c.axvline(QH_THRESHOLD_PCT, color="#c0392b", lw=0.8, ls="--", zorder=1)
    ax_c.text(1.25, -0.5, "QH unreliable", color="#c0392b", fontsize=6.5, va="bottom")
    rows_c = []
    for yi, p in zip(y, pairs):
        lo, hi = imag[(p["system"], p["low"])], imag[(p["system"], p["high"])]
        ax_c.plot([lo, hi], [yi, yi], color=p["color"], lw=1.0, zorder=2)
        ax_c.plot(lo, yi, "o", mfc="white", mec=p["color"], mew=1.3, ms=5, zorder=3)
        ax_c.plot(hi, yi, "o", color=p["color"], ms=5, zorder=3)
        for xv, ph in ((lo, p["low"]), (hi, p["high"])):
            tag = PHASE_SHORT[ph] + ("†" if ph == "quartz_beta" else "")
            ax_c.annotate(tag, (xv, yi), xytext=(0, 5), textcoords="offset points",
                          ha="center", va="bottom", fontsize=6)
        rows_c += [dict(pair=p["label"], phase=ph, role=role, imaginary_pct=v,
                        qh_reliable=v <= QH_THRESHOLD_PCT, minimum_frequency_THz=minfreq[(p["system"], ph)])
                   for ph, role, v in ((p["low"], "low_T", lo), (p["high"], "high_T", hi))]
    ax_c.set_xscale("symlog", linthresh=0.1, linscale=0.6)
    ax_c.set_xlim(-0.03, 60)
    ax_c.set_xticks([0, 0.1, 1, 10])
    ax_c.set_xticklabels(["0", "0.1", "1", "10"])
    ax_c.set_xlabel("Imaginary modes (%)")
    ax_c.tick_params(axis="y", left=False, labelleft=False)
    ax_c.spines["left"].set_visible(False)

    for ax, lab in ((ax_a, "(a)"), (ax_b, "(b)"), (ax_c, "(c)")):
        ax.text(-0.02, 1.04, lab, transform=ax.transAxes, fontsize=9, fontweight="bold", ha="right", va="bottom")

    # every pair: at least one phase above threshold
    assert all(max(imag[(p["system"], p["low"])], imag[(p["system"], p["high"])]) > QH_THRESHOLD_PCT for p in pairs)

    for ext in ("pdf", "png", "svg"):
        fig.savefig(out / f"figure1.{ext}", dpi=300, bbox_inches="tight")
    pd.DataFrame(rows_a).to_csv(out / "panel_a.csv", index=False)
    pd.DataFrame([{k: p[k] for k in ("label", "system", "low", "high", "Tc", "slope", "slope_err",
                                      "r2", "n_fit", "K_per_meV", "K_per_meV_err")} for p in pairs]
                 ).rename(columns={"slope": "slope_meV_per_atom_per_K"}).to_csv(out / "panel_b.csv", index=False)
    pd.DataFrame(rows_c).to_csv(out / "panel_c.csv", index=False)
    plt.close(fig)

    # Provenance is generated here so the figure remains reproducible without
    # changing the supplied plotting/data-selection logic above.
    input_paths = [root / "result/tables/phase_inventory.csv"]
    for sysdir in sorted((root / "data/processed").iterdir()):
        system_json = sysdir / "system.json"
        if not system_json.exists():
            continue
        input_paths.append(system_json)
        phases = json.loads(system_json.read_text(encoding="utf-8"))["phases"]
        input_paths.extend(sysdir / phase / "reference_G.csv" for phase in phases)
    input_paths = sorted({path.resolve() for path in input_paths if path.exists()})
    try:
        commit = subprocess.check_output(
            ["git", "-C", str(root), "rev-parse", "HEAD"], text=True
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        commit = "unknown"
    meta = {
        "script": "src/figures/figure1_difficulty.py",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "git_commit": commit,
        "fit_half_width_K": FIT_HALF_WIDTH_K,
        "input_sha256": [
            {"path": str(path.relative_to(root)), "sha256": sha256(path)}
            for path in input_paths
        ],
        "pair_order": [p["label"] for p in pairs],
        "outputs": ["figure1.pdf", "figure1.png", "figure1.svg", "panel_a.csv", "panel_b.csv", "panel_c.csv", "caption.md"],
    }
    (out / "figure1.meta.json").write_text(json.dumps(meta, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    panel_b = pd.read_csv(out / "panel_b.csv")
    panel_c = pd.read_csv(out / "panel_c.csv")
    max_shift = panel_b.loc[panel_b["K_per_meV"].idxmax()]
    max_imag = float(panel_c["imaginary_pct"].max())
    caption = (
        "Figure 1. Why relative free energies are hard to learn. (a) Reference relative "
        "free energies are aligned at each transition temperature and shown against "
        "T − T_c; the shaded band is ±1 meV/atom. (b) A 1 meV/atom error produces a "
        f"maximum shift of {max_shift['K_per_meV']:.1f} K at {max_shift['label']}, "
        "using the inverse absolute slope from a linear fit over T_c ± 25 K. (c) "
        f"Low-T phases are open markers and high-T phases filled; the largest plotted "
        f"imaginary-mode fraction is {max_imag:.2f}%. The dashed 1% line marks the "
        "QH-reliability threshold. The near-degenerate β-cristobalite–β-tridymite "
        "pair has no reference crossing and is omitted. Abbreviations: qz = quartz, "
        "crs = cristobalite, and trd = tridymite."
    )
    (out / "caption.md").write_text(caption + "\n", encoding="utf-8")
    print(pd.read_csv(out / "panel_b.csv").to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
