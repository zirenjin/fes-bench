"""Exploratory (post-hoc): relaxation-channel decisions reported over ALL primary pairs, not only the successful subset.

Pairs are classed as evaluable (both entries pass the pipeline) or failed (any entry missing, not converged, space-group
change, or unstable at the MLIP geometry). Over all non-degenerate pairs (|dF_DFT| >= 1 meV/atom):
  fixed geometry on all pairs; fixed geometry on the evaluable subset; pipeline on the evaluable subset;
  pipeline over all pairs with failures as non-decisions (coverage, risk) and, as a bound, failures counted as wrong.
Also: is the failed set enriched in hard pairs (fixed-geometry |e| on failed vs evaluable)?
Held-out check of the R5 feature: LOFO conformal abstention at alpha = 0.05, plain (|err|) vs normalised
(|err| / (q0 + |d_static|), q0 = median |d_static| over the other formulas). Pair-weighted.
Usage: relax_coverage.py <FINAL_GROUPS.csv> <relax_dir> <fixed_dir> <out.json>
"""
import csv, itertools, json, sys
from pathlib import Path
import numpy as np

groups_csv, rdir, fdir, out = map(Path, sys.argv[1:5])
TS = [300, 1000, 1500]
groups = {r["formula"]: r["mp_ids"].split(";") for r in csv.DictReader(open(groups_csv))}
E = {}
for ms in groups.values():
    for m in ms:
        z = np.load(fdir / f"{m}.npz"); T = z["T"]; ix = {t: int(np.flatnonzero(T == t)[0]) for t in TS}
        e = dict(dft={t: float(z["E_dft_eV_atom"]) * 1000 + float(z["F_dft"][ix[t]]) for t in TS},
                 fixed={t: float(z["E_mlip_eV_atom"]) * 1000 + float(z["F_mlip"][ix[t]]) for t in TS}, ok=False)
        p = rdir / f"{m}.npz"
        if p.exists():
            r = np.load(p)
            e["ok"] = bool(r["converged"]) and str(r["sg_before"]) == str(r["sg_after"]) and float(r["fmin_relaxed"]) >= -0.05
            e["pipe"] = {t: float(r["E_relaxed_mlip"]) + float(r["F_vib_relaxed"][ix[t]]) for t in TS}
            e["static"] = float(r["E_relaxed_mlip"]) - float(r["E_dftgeom_mlip"])
        E[m] = e
P = [(f, i, j) for f, ms in groups.items() for i, j in itertools.combinations(ms, 2)]
F = np.array([f for f, _, _ in P]); ok = np.array([E[i]["ok"] and E[j]["ok"] for _, i, j in P])


def lofo_cov(err, pred, D, nd, unit, alpha=0.05):
    acc = np.zeros(len(err), bool)
    for f in set(F[nd]):
        tr, te = (F != f) & nd, (F == f) & nd
        thr = np.quantile(np.abs(err[tr]) / unit[tr], 1 - alpha); acc[te] = np.abs(pred[te]) > thr * unit[te]
    wrong = np.sign(pred) != np.sign(D)
    return dict(coverage=float(acc[nd].mean()), selective_risk=float((acc & wrong)[nd].sum() / max(acc[nd].sum(), 1)))


res = dict(n_pairs_all=len(P), n_pairs_evaluable=int(ok.sum()), n_formulas_all=len(set(F)), n_formulas_with_evaluable=len(set(F[ok])))
for t in TS:
    D = np.array([E[i]["dft"][t] - E[j]["dft"][t] for _, i, j in P]); nd = np.abs(D) >= 1.0
    ef = np.array([(E[i]["fixed"][t] - E[j]["fixed"][t]) - d for (_, i, j), d in zip(P, D)])
    ep = np.array([((E[i]["pipe"][t] - E[j]["pipe"][t]) - d) if o else np.nan for (_, i, j), d, o in zip(P, D, ok)])
    wf = (np.sign(D + ef) != np.sign(D)); wp = (np.sign(D + np.nan_to_num(ep)) != np.sign(D))
    blk = dict(n_nondegenerate_all=int(nd.sum()), n_nondegenerate_evaluable=int((nd & ok).sum()),
               wrong_fixed_all=float(wf[nd].mean()), wrong_fixed_evaluable=float(wf[nd & ok].mean()), wrong_fixed_failed=float(wf[nd & ~ok].mean()),
               wrong_pipeline_evaluable=float(wp[nd & ok].mean()),
               pipeline_over_all_failures_as_nondecision=dict(decided=float((nd & ok).sum() / nd.sum()), wrong_among_all=float((wp & ok)[nd].mean())),
               pipeline_over_all_failures_as_wrong=float(((wp & ok) | ~ok)[nd].mean()),
               fixed_abs_err_median=dict(evaluable=float(np.median(np.abs(ef[ok]))), failed=float(np.median(np.abs(ef[~ok])))),
               fixed_rms=dict(all=float(np.sqrt(np.mean(ef ** 2))), evaluable=float(np.sqrt(np.mean(ef[ok] ** 2))), failed=float(np.sqrt(np.mean(ef[~ok] ** 2)))))
    # R5 held-out check, evaluable pairs only
    st = np.array([abs(E[i]["static"] - E[j]["static"]) if o else np.nan for (_, i, j), o in zip(P, ok)])
    m = ok & nd; epp = np.where(ok, ep, 0.0); pred = D + epp
    q0 = np.zeros(len(P))
    for f in set(F[m]):
        q0[F == f] = np.median(st[m & (F != f)])
    blk["R5_heldout_alpha0.05"] = dict(plain=lofo_cov(epp, pred, D, m, np.ones(len(P))), normalised=lofo_cov(epp, pred, D, m, np.where(ok, q0 + np.nan_to_num(st), 1.0)))
    res[f"T{t}"] = blk
out.write_text(json.dumps(res, indent=1) + "\n"); print(json.dumps(res, indent=1))
