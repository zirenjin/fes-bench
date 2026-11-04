"""PROTOCOL.md dataset 1: Alexandria fixed-geometry harmonic pairs, F0, arms none / static / Tc / T0ref(300 K).
Usage: anchored_alexandria.py <FINAL_GROUPS.csv> <F0_mlip_dir> <out.json>
"""
import csv, itertools, json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).parent))
from anchor_common import crossings, match_events, tc_anchor, side_correct  # noqa: E402

groups_csv, mdir, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
TE = [600, 1000, 1500]; LO, HI = 50.0, 1500.0; rng = np.random.default_rng(20261102)
groups = {r["formula"]: r["mp_ids"].split(";") for r in csv.DictReader(open(groups_csv))}
E = {}
for ms in groups.values():
    for m in ms:
        z = np.load(mdir / f"{m}.npz"); E[m] = dict(T=z["T"], ed=float(z["E_dft_eV_atom"]) * 1000, em=float(z["E_mlip_eV_atom"]) * 1000,
                                                    Fd=z["F_dft"], Fm=z["F_mlip"])
T = E[next(iter(E))]["T"]; it = {t: int(np.flatnonzero(T == t)[0]) for t in TE + [300]}
rows = []
for f, ms in groups.items():
    for i, j in itertools.combinations(ms, 2):
        a, b = E[i], E[j]
        dref = (a["ed"] + a["Fd"]) - (b["ed"] + b["Fd"]); d0 = (a["em"] + a["Fm"]) - (b["em"] + b["Fm"])
        arms = {"none": d0, "static": (a["ed"] + a["Fm"]) - (b["ed"] + b["Fm"]), "T0ref": d0 - d0[it[300]] + dref[it[300]]}
        tc, dtc = tc_anchor(T, dref, d0, LO)
        if tc is not None and tc <= HI:
            arms["Tc"] = dtc
        rows.append(dict(formula=f, dref=dref, arms=arms, tc=tc if (tc is not None and tc <= HI) else None,
                         dS=-np.gradient(d0 - dref, T) * 1000))
F = np.array([r["formula"] for r in rows]); forms = sorted(set(F))


def cluster_ci(v, mask):
    fs = [f for f in forms if ((F == f) & mask).any()]
    num = {f: v[(F == f) & mask].sum() for f in fs}; den = {f: ((F == f) & mask).sum() for f in fs}
    d = [sum(num[f] for f in s) / max(sum(den[f] for f in s), 1) for s in (rng.choice(fs, len(fs)) for _ in range(2000))]
    return [float(x) for x in np.percentile(d, [2.5, 97.5])]


res = dict(n_pairs=len(rows), n_formulas=len(forms),
           applicability={a: float(np.mean([a in r["arms"] for r in rows])) for a in ("none", "static", "Tc", "T0ref")},
           response_dS_rms_ueV={str(t): float(np.sqrt(np.mean([r["dS"][it[t]] ** 2 for r in rows]))) for t in TE})
for arm in ("none", "static", "Tc", "T0ref"):
    app = np.array([arm in r["arms"] for r in rows]); blk = dict(n_applicable=int(app.sum()))
    for t in TE:
        D = np.array([r["dref"][it[t]] for r in rows]); P = np.array([r["arms"][arm][it[t]] if arm in r["arms"] else np.nan for r in rows])
        nd = app & (np.abs(D) >= 1.0)
        if arm == "Tc":
            nd &= np.array([r["tc"] is None or abs(t - r["tc"]) >= 100 for r in rows])
        err = P - D; wrong = (np.sign(P) != np.sign(D)).astype(float)
        acc = np.zeros(len(rows), bool)
        for f in set(F[nd]):
            tr, te = (F != f) & nd, (F == f) & nd
            if tr.sum() > 20:
                thr = np.quantile(np.abs(err[tr]), 0.95); acc[te] = np.abs(P[te]) > thr
        blk[f"T{t}"] = dict(n_nondegenerate=int(nd.sum()), rms_e=float(np.sqrt(np.nanmean(err[app] ** 2))),
                            wrong=float(wrong[nd].mean()), wrong_ci95=cluster_ci(wrong, nd),
                            coverage_a005=float(acc[nd].mean()), selective_risk=float((acc & (wrong > 0))[nd].sum() / max(acc[nd].sum(), 1)))
    ev = [match_events(crossings(T[(T >= LO)], r["dref"][T >= LO]), crossings(T[T >= LO], r["arms"][arm][T >= LO]),
                       provided=r["tc"] if arm == "Tc" else None) for r in rows if arm in r["arms"]]
    blk["events"] = {k: int(sum(e[k] for e in ev)) for k in ("reference", "predicted", "hits", "missed", "false")}
    if arm == "Tc":
        sc = [side_correct(T, r["dref"], r["arms"]["Tc"], r["tc"], LO, HI) for r in rows if "Tc" in r["arms"]]
        sc = [s for s in sc if s is not None]
        blk["side_correct_mean"] = float(np.mean(sc)); blk["side_correct_all_points_frac"] = float(np.mean([s == 1.0 for s in sc]))
        blk["anchored_pair_extra_crossings"] = int(sum(len([c for c in crossings(T[T >= LO], r["arms"]["Tc"][T >= LO]) if abs(c - r["tc"]) > 100])
                                                       for r in rows if "Tc" in r["arms"]))
    res[arm] = blk
out.write_text(json.dumps(res, indent=1) + "\n"); print(json.dumps(res, indent=1))
