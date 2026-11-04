"""PROTOCOL.md (relaxation channel) analyses R1-R5.

Per entry (meV/atom): G_dft = E_dft + F_dft at the DFT geometry; G_fixed = MLIP E + F_vib at the DFT geometry (existing
fixed-geometry run); G_pipe = MLIP E + F_vib at the MLIP-relaxed geometry. Per pair: e_fixed = dG_fixed - dG_dft,
d_relax = dG_pipe - dG_fixed (static part d(E_relaxed - E_at_DFT_geometry), vibrational part the remainder),
e_pipe = e_fixed + d_relax. Formula-weighted means, cluster bootstrap over formulas (2000 draws).
Usage: analyze_relax.py <FINAL_GROUPS.csv> <relax_dir> <fixed_dir> <out.json>
"""
import csv, itertools, json, sys
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr

groups_csv, relax_dir, fixed_dir, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), Path(sys.argv[4])
TS = [300, 1000, 1500]; rng = np.random.default_rng(20261101)
groups = {r["formula"]: r["mp_ids"].split(";") for r in csv.DictReader(open(groups_csv))}
E, excl = {}, {"missing_or_failed": [], "not_converged": [], "spacegroup_changed": [], "unstable_at_mlip_geometry": []}
for ms in groups.values():
    for m in ms:
        pr = relax_dir / f"{m}.npz"
        if not pr.exists():
            excl["missing_or_failed"].append(m); continue
        r, z = np.load(pr), np.load(fixed_dir / f"{m}.npz")
        T = z["T"]; ix = {t: int(np.flatnonzero(T == t)[0]) for t in TS}
        bad = None if bool(r["converged"]) else "not_converged"
        bad = bad or (None if str(r["sg_before"]) == str(r["sg_after"]) else "spacegroup_changed")
        bad = bad or (None if float(r["fmin_relaxed"]) >= -0.05 else "unstable_at_mlip_geometry")
        if bad:
            excl[bad].append(m); continue
        e_fix = float(z["E_mlip_eV_atom"]) * 1000
        assert abs(e_fix - float(r["E_dftgeom_mlip"])) < 0.5, (m, e_fix, float(r["E_dftgeom_mlip"]))  # same model, same geometry
        E[m] = dict(dft={t: float(z["E_dft_eV_atom"]) * 1000 + float(z["F_dft"][ix[t]]) for t in TS},
                    fixed={t: e_fix + float(z["F_mlip"][ix[t]]) for t in TS},
                    pipe={t: float(r["E_relaxed_mlip"]) + float(r["F_vib_relaxed"][ix[t]]) for t in TS},
                    static=float(r["E_relaxed_mlip"]) - float(r["E_dftgeom_mlip"]), R0=float(r["R0_rms"]),
                    dV=float(r["vol_relaxed"]) / float(r["vol_dft"]) - 1)
P = []
for f, ms in groups.items():
    ok = [m for m in ms if m in E]
    for i, j in itertools.combinations(ok, 2):
        row = dict(formula=f, i=i, j=j, static=E[i]["static"] - E[j]["static"], R0max=max(E[i]["R0"], E[j]["R0"]))
        for t in TS:
            d = lambda k: E[i][k][t] - E[j][k][t]
            row[f"ref{t}"] = d("dft"); row[f"efix{t}"] = d("fixed") - d("dft"); row[f"relax{t}"] = d("pipe") - d("fixed")
            row[f"epipe{t}"] = d("pipe") - d("dft"); row[f"relaxvib{t}"] = row[f"relax{t}"] - row["static"]
            row[f"srel_i{t}"] = E[i]["pipe"][t] - E[i]["fixed"][t]; row[f"srel_j{t}"] = E[j]["pipe"][t] - E[j]["fixed"][t]
        P.append(row)
F = np.array([p["formula"] for p in P]); forms = sorted(set(F)); idx = {f: np.flatnonzero(F == f) for f in forms}
col = lambda k: np.array([p[k] for p in P])
draws = [rng.choice(forms, len(forms)) for _ in range(2000)]


def fw(v, fs=forms, mask=None):
    m = np.ones(len(P), bool) if mask is None else mask
    vals = [v[idx[f]][m[idx[f]]].mean() for f in fs if m[idx[f]].any()]
    return float(np.mean(vals)) if vals else float("nan")


def boot(fn):
    return dict(value=fn(forms), ci95=[float(x) for x in np.nanpercentile([fn(s) for s in draws], [2.5, 97.5])])


def coverage(err, pred, D, nd, alpha=0.05):
    thr = np.array([np.quantile(np.abs(err[F != F[k]]), 1 - alpha) for k in range(len(P))])
    acc = (np.abs(pred) > thr) & nd
    return dict(coverage=float(acc.sum() / nd.sum()), selective_risk=float((acc & (np.sign(pred) != np.sign(D))).sum() / max(acc.sum(), 1)))


res = dict(n_entries_used=len(E), n_pairs=len(P), n_formulas=len(forms), exclusions={k: len(v) for k, v in excl.items()}, excluded_ids=excl,
           entry_volume_change_pct=dict(median=float(np.median([100 * e["dV"] for e in E.values()])),
                                        p10_p90=[float(np.percentile([100 * e["dV"] for e in E.values()], q)) for q in (10, 90)]),
           entry_static_relaxation_meV=dict(median=float(np.median([e["static"] for e in E.values()]))))
st = col("static")
for t in TS:
    ef, rl, ep, rv, D = col(f"efix{t}"), col(f"relax{t}"), col(f"epipe{t}"), col(f"relaxvib{t}"), col(f"ref{t}")
    si, sj = col(f"srel_i{t}"), col(f"srel_j{t}"); nd = np.abs(D) >= 1.0
    rmsb = lambda v: boot(lambda fs: np.sqrt(fw(v ** 2, fs)))
    blk = dict(
        R1_rms=dict(e_fixed=rmsb(ef), d_relax=rmsb(rl), d_relax_static=rmsb(st), d_relax_vib=rmsb(rv), e_pipeline=rmsb(ep)),
        R2_relax_cancellation_c=boot(lambda fs: 1 - fw(rl ** 2, fs) / fw(si ** 2 + sj ** 2, fs)),
        R3=dict(uncentred_corr_relax_vs_efixed=float(np.sum(rl * ef) / np.sqrt(np.sum(rl ** 2) * np.sum(ef ** 2))),
                frac_pipeline_smaller=boot(lambda fs: fw((np.abs(ep) < np.abs(ef)).astype(float), fs))),
        R4=dict(wrong_fixed=boot(lambda fs: fw((np.sign(D + ef) != np.sign(D)).astype(float), fs, nd)),
                wrong_pipeline=boot(lambda fs: fw((np.sign(D + ep) != np.sign(D)).astype(float), fs, nd)),
                coverage_fixed=coverage(ef, D + ef, D, nd), coverage_pipeline=coverage(ep, D + ep, D, nd)),
        R5_spearman_abs_epipe=dict(abs_static_relax_diff=float(spearmanr(np.abs(st), np.abs(ep)).correlation),
                                   R0_pair_max=float(spearmanr(col("R0max"), np.abs(ep)).correlation)),
        n_nondegenerate=int(nd.sum()))
    res[f"T{t}"] = blk
out.write_text(json.dumps(res, indent=1) + "\n"); print(json.dumps(res, indent=1))
