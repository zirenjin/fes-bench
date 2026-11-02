"""Thermodynamic error-covariance audit (PROTOCOL.md Q1-Q4). Same-geometry, fixed-volume harmonic F = E0 + F_vib.

Inputs: FINAL_GROUPS.csv from select_entries.py and one MLIP result directory (per-entry .npz from
mlip_vs_dft_phonons.py). All errors are MLIP minus DFT in meV/atom.
Implementation notes fixed before running:
- pair orientation (i - j) is arbitrary, so Spearman correlations use the symmetrized set (both orientations) and an
  orientation-invariant uncentred correlation is also reported;
- the cancellation fraction uses raw second moments, c = 1 - <dD^2> / <d_i^2 + d_j^2>, so a shared bias counts as
  common mode (that is the quantity a phase-relative method can exploit);
- formula-weighted means (primary) average within each formula first; cluster bootstrap over formulas, 2000 draws.
Usage: analyze_covariance.py <FINAL_GROUPS.csv> <mlip_dir> <out_dir> [--exclude-large N]
"""
import csv, itertools, json, sys
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr

groups_csv, mlip_dir, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]); out.mkdir(parents=True, exist_ok=True)
EXCL = int(sys.argv[sys.argv.index("--exclude-large") + 1]) if "--exclude-large" in sys.argv else None
TS = [300, 600, 1000, 1500]; rng = np.random.default_rng(20261031)

groups = {r["formula"]: r["mp_ids"].split(";") for r in csv.DictReader(open(groups_csv))}
if EXCL:
    groups = {f: m for f, m in groups.items() if len(m) < EXCL}
entry, missing = {}, []
for f, ms in groups.items():
    for m in ms:
        p = mlip_dir / f"{m}.npz"
        if not p.exists():
            missing.append(m); continue
        z = np.load(p); T = z["T"]; idx = {t: int(np.flatnonzero(T == t)[0]) for t in TS}
        Fd, Fm = z["F_dft"], z["F_mlip"]
        entry[m] = dict(formula=f, dE0=(float(z["E_mlip_eV_atom"]) - float(z["E_dft_eV_atom"])) * 1000, dZPE=float(Fm[0] - Fd[0]),
                        dFth={t: float((Fm[idx[t]] - Fm[0]) - (Fd[idx[t]] - Fd[0])) for t in TS},
                        Fref={t: float(z["E_dft_eV_atom"]) * 1000 + float(Fd[idx[t]]) for t in TS})
pairs = []
for f, ms in groups.items():
    ok = [m for m in ms if m in entry]
    for i, j in itertools.combinations(ok, 2):
        a0 = (entry[i]["dE0"] - entry[j]["dE0"]) + (entry[i]["dZPE"] - entry[j]["dZPE"])
        pairs.append(dict(formula=f, i=i, j=j, dE0=entry[i]["dE0"] - entry[j]["dE0"], dZPE=entry[i]["dZPE"] - entry[j]["dZPE"], a=a0,
                          **{f"b{t}": entry[i]["dFth"][t] - entry[j]["dFth"][t] for t in TS},
                          **{f"ref{t}": entry[i]["Fref"][t] - entry[j]["Fref"][t] for t in TS}))
formulas = sorted({p["formula"] for p in pairs})
byf = {f: [p for p in pairs if p["formula"] == f] for f in formulas}


def fw_mean(fn, fs):
    vals = [fn(byf[f]) for f in fs]; vals = [v for v in vals if v is not None and np.isfinite(v)]
    return float(np.mean(vals)) if vals else float("nan")


def boot(fn):
    point = fw_mean(fn, formulas); draws = [fw_mean(fn, list(rng.choice(formulas, len(formulas)))) for _ in range(2000)]
    lo, hi = np.nanpercentile(draws, [2.5, 97.5]); return dict(value=point, ci95=[float(lo), float(hi)])


def pair_weighted(fn):
    return float(fn(pairs))


def cancel(comp, t=None):
    """1 - <dD^2>/<d_i^2 + d_j^2> within a set of pairs, for component comp."""
    def single(m):
        e = entry[m]; return e["dFth"][t] if comp == "Fth" else (e["dE0"] + e["dZPE"] if comp == "anchor" else e[comp])
    def fn(ps):
        num = np.mean([(single(p["i"]) - single(p["j"])) ** 2 for p in ps]); den = np.mean([single(p["i"]) ** 2 + single(p["j"]) ** 2 for p in ps])
        return 1 - num / den if den > 0 else None
    return fn


res = {"n_formulas": len(formulas), "n_pairs": len(pairs), "n_entries": len(entry), "missing_mlip": missing, "exclude_large_groups_ge": EXCL}
# Q1
q1 = {"dE0": boot(cancel("dE0")), "dZPE": boot(cancel("dZPE")), "anchor(dE0+dZPE)": boot(cancel("anchor"))}
for t in TS: q1[f"dFth_{t}K"] = boot(cancel("Fth", t))
q1["rms_single_vs_pair_meV"] = {c: dict(single=float(np.sqrt(np.mean([entry[m][c] ** 2 for m in entry]))), pair=float(np.sqrt(np.mean([p[c] ** 2 for p in pairs]))))
                                for c in ("dE0", "dZPE")}
for t in TS:
    q1["rms_single_vs_pair_meV"][f"dFth_{t}K"] = dict(single=float(np.sqrt(np.mean([entry[m]["dFth"][t] ** 2 for m in entry]))),
                                                     pair=float(np.sqrt(np.mean([p[f"b{t}"] ** 2 for p in pairs]))))
res["Q1_common_mode"] = q1
# Q2-Q4
for t in TS:
    b = f"b{t}"; r = f"ref{t}"
    opp = lambda ps: np.mean([np.sign(p["a"]) * np.sign(p[b]) < 0 for p in ps])
    cancels = lambda ps: np.mean([abs(p["a"] + p[b]) < max(abs(p["a"]), abs(p[b])) for p in ps])
    harm_anchor = lambda ps: np.mean([abs(p[b]) > abs(p["a"] + p[b]) for p in ps])      # correct anchor only -> error b
    harm_thermal = lambda ps: np.mean([abs(p["a"]) > abs(p["a"] + p[b]) for p in ps])   # correct thermal only -> error a
    nd = [p for p in pairs if abs(p[r]) >= 1.0]
    def wrong(kind):
        def fn(ps):
            ps = [p for p in ps if abs(p[r]) >= 1.0]
            if not ps: return None
            err = {"none": lambda p: p["a"] + p[b], "anchor_only_corrected": lambda p: p[b], "thermal_only_corrected": lambda p: p["a"]}[kind]
            return np.mean([np.sign(p[r] + err(p)) != np.sign(p[r]) for p in ps])
        return fn
    A = np.array([p["a"] for p in pairs]); B = np.array([p[b] for p in pairs])
    res[f"T{t}"] = dict(
        Q2_opposite_sign_fraction=boot(opp), Q2_partial_cancellation_fraction=boot(cancels),
        Q2_spearman_symmetrized=float(spearmanr(np.r_[A, -A], np.r_[B, -B]).correlation),
        Q2_uncentred_corr=float(np.sum(A * B) / np.sqrt(np.sum(A * A) * np.sum(B * B))),
        Q3_anchor_only_correction_increases_error=boot(harm_anchor), Q3_thermal_only_correction_increases_error=boot(harm_thermal),
        Q3_ordering_wrong=dict(none=boot(wrong("none")), anchor_only_corrected=boot(wrong("anchor_only_corrected")),
                               thermal_only_corrected=boot(wrong("thermal_only_corrected"))),
        n_pairs_nondegenerate=len(nd), n_pairs_near_degenerate=len(pairs) - len(nd),
        pair_weighted=dict(opposite_sign=pair_weighted(opp), partial_cancellation=pair_weighted(cancels),
                           anchor_only_harm=pair_weighted(harm_anchor), thermal_only_harm=pair_weighted(harm_thermal)),
        rms_a_meV=float(np.sqrt(np.mean(A ** 2))), rms_b_meV=float(np.sqrt(np.mean(B ** 2))), rms_total_meV=float(np.sqrt(np.mean((A + B) ** 2))))
res["Q4_sign_of_uncentred_corr_by_T"] = {t: float(np.sign(res[f"T{t}"]["Q2_uncentred_corr"])) for t in TS}
with open(out / "PAIRS.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(pairs[0])); w.writeheader(); w.writerows(pairs)
with open(out / "ENTRIES.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["mp_id", "formula", "dE0", "dZPE"] + [f"dFth_{t}K" for t in TS])
    for m, e in entry.items(): w.writerow([m, e["formula"], e["dE0"], e["dZPE"]] + [e["dFth"][t] for t in TS])
(out / "RESULTS.json").write_text(json.dumps(res, indent=1) + "\n")
print(json.dumps(res, indent=1))
