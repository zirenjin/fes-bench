"""PROTOCOL.md parts 1-2: root-finding sensitivity and error-type triage on reference-crossing pairs.
Usage: tc_sensitivity.py <FINAL_GROUPS.csv> <curves_F0.npz> <curves_MS.npz> <out.json>
"""
import csv, itertools, json, sys
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr

groups_csv, out = Path(sys.argv[1]), Path(sys.argv[-1])
LO, HI, GATE, NB = 50.0, 1500.0, -0.05, 2000
groups = {r["formula"]: r["mp_ids"].split(";") for r in csv.DictReader(open(groups_csv))}


def crossings(T, d):
    s = np.sign(d); k = np.flatnonzero(s[:-1] * s[1:] < 0)
    return [float(T[i] - d[i] * (T[i + 1] - T[i]) / (d[i + 1] - d[i])) for i in k]


def load(p):
    z = np.load(p); ix = {m: k for k, m in enumerate(z["ids"])}
    return z["T"], {m: dict(Gd=z["E_dft"][k] + z["F_dft"][k], Gm=z["E_mlip"][k] + z["F_mlip"][k], Gs=z["E_dft"][k] + z["F_mlip"][k],
                            fd=float(z["fmin_dft"][k]), fm=float(z["fmin_mlip"][k])) for m, k in ix.items()}


def nearest(cs, tc):
    return min(cs, key=lambda c: abs(c - tc)) if cs else None


def triage(T, E):
    rows = []; sel = (T >= LO) & (T <= HI); Ts = T[sel]
    for f, ms in groups.items():
        for i, j in itertools.combinations(ms, 2):
            a, b = E[i], E[j]; dref = a["Gd"] - b["Gd"]; cr = crossings(Ts, dref[sel])
            if not cr:
                continue
            tc = cr[0]; dmod = a["Gm"] - b["Gm"]; dst = a["Gs"] - b["Gs"]
            S_ref = -np.interp(tc, T, np.gradient(dref, T)) * 1000; S_mod = -np.interp(tc, T, np.gradient(dmod, T)) * 1000  # ueV/(atom K)
            e = float(np.interp(tc, T, dmod)); e_st = float(np.interp(tc, T, dst))
            tm = nearest(crossings(Ts, dmod[sel]), tc); ts = nearest(crossings(Ts, dst[sel]), tc)
            anch = dmod - e; far = sel & (np.abs(T - tc) >= 100)
            resp_bad = bool(np.any(np.sign(anch[far]) != np.sign(dref[far])))
            applic = min(a["fd"], b["fd"], a["fm"], b["fm"]) < GATE
            dT = None if tm is None else tm - tc
            cls = "A_applicability" if applic else ("B_correct" if dT is not None and abs(dT) <= 50 else ("C_response" if resp_bad else "D_location"))
            rows.append(dict(formula=f, i=i, j=j, tc_ref=tc, S_ref=float(S_ref), S_mod=float(S_mod), e=e, e_static=e_st,
                             eps50=abs(S_ref) * 50 / 1000, eps100=abs(S_ref) * 100 / 1000, dT=dT,
                             dT_lin=e / S_ref * 1000 if S_ref else None, dT_resp=e / S_mod * 1000 if S_mod else None,
                             dT_static=None if ts is None else ts - tc, cls=cls, S_sign_wrong=bool(np.sign(S_mod) != np.sign(S_ref)),
                             static_fixes=(ts is not None and abs(ts - tc) <= 50)))
    return rows


def fw_boot(rows, fn, seed=20261102):
    """Formula-weighted fraction with a formula-cluster bootstrap; fn(row) -> 0/1 or None (excluded)."""
    byf = {}
    for r in rows:
        v = fn(r)
        if v is not None: byf.setdefault(r["formula"], []).append(float(v))
    fs = sorted(byf); m = {f: np.mean(byf[f]) for f in fs}; rng = np.random.default_rng(seed)
    d = [np.mean([m[f] for f in rng.choice(fs, len(fs))]) for _ in range(NB)]
    return dict(value=float(np.mean([m[f] for f in fs])), ci95=[float(x) for x in np.percentile(d, [2.5, 97.5])], n_pairs=sum(len(v) for v in byf.values()), n_formulas=len(fs))


res = dict(protocol_sha256=open(Path(__file__).parents[1] / "PROTOCOL.sha256").read().split()[0])
all_rows = {}
for name, p in (("F0", sys.argv[2]), ("MatterSim", sys.argv[3])):
    T, E = load(p); rows = triage(T, E); all_rows[name] = rows
    eps = np.array([r["eps50"] for r in rows]); S = np.abs([r["S_ref"] for r in rows]); e = np.abs([r["e"] for r in rows])
    hit = [r for r in rows if r["dT"] is not None]
    blk = dict(n_pairs=len(rows), n_formulas=len({r["formula"] for r in rows}),
               S_ref_abs_ueV=dict(median=float(np.median(S)), p25=float(np.percentile(S, 25)), p75=float(np.percentile(S, 75))),
               eps50_meV=dict(median=float(np.median(eps)), p25=float(np.percentile(eps, 25)), p75=float(np.percentile(eps, 75))),
               e_abs_meV=dict(median=float(np.median(e)), p25=float(np.percentile(e, 25)), p75=float(np.percentile(e, 75))),
               within_budget_50K=fw_boot(rows, lambda r: abs(r["e"]) <= r["eps50"]),
               within_budget_100K=fw_boot(rows, lambda r: abs(r["e"]) <= r["eps100"]),
               within_budget_50K_static=fw_boot(rows, lambda r: abs(r["e_static"]) <= r["eps50"]),
               model_crossing_found=fw_boot(rows, lambda r: r["dT"] is not None),
               linear_check=dict(n=len(hit),
                                 spearman_dT_vs_lin=float(spearmanr([r["dT"] for r in hit], [r["dT_lin"] for r in hit]).correlation) if len(hit) > 2 else None,
                                 frac_within_50K_lin=float(np.mean([abs(r["dT"] - r["dT_lin"]) <= 50 for r in hit])) if hit else None,
                                 frac_within_50K_resp=float(np.mean([abs(r["dT"] - r["dT_resp"]) <= 50 for r in hit])) if hit else None,
                                 median_abs_dT=float(np.median([abs(r["dT"]) for r in hit])) if hit else None),
               classes={c: fw_boot(rows, lambda r, c=c: r["cls"] == c) for c in ("A_applicability", "B_correct", "C_response", "D_location")},
               class_counts={c: sum(r["cls"] == c for r in rows) for c in ("A_applicability", "B_correct", "C_response", "D_location")},
               C_with_wrong_S_sign=sum(r["cls"] == "C_response" and r["S_sign_wrong"] for r in rows),
               D_fixed_by_static=sum(r["cls"] == "D_location" and r["static_fixes"] for r in rows),
               missed_by_class={c: sum(r["cls"] == c and r["dT"] is None for r in rows) for c in ("A_applicability", "C_response", "D_location")})
    res[name] = blk
# agreement of the two models' class labels on the same pairs (excluding A, whose gate is model-specific)
key = lambda r: (r["i"], r["j"])
f0 = {key(r): r["cls"] for r in all_rows["F0"]}; msd = {key(r): r["cls"] for r in all_rows["MatterSim"]}
res["class_crosstab_F0_rows_MS_cols"] = {a: {b: sum(f0[k] == a and msd[k] == b for k in f0) for b in sorted(set(msd.values()))} for a in sorted(set(f0.values()))}
out.write_text(json.dumps(res, indent=1) + "\n")
with open(out.with_suffix(".pairs.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["model"] + list(all_rows["F0"][0])); w.writeheader()
    for n, rs in all_rows.items():
        for r in rs: w.writerow(dict(model=n, **r))
print(json.dumps({k: v for k, v in res.items()}, indent=1))
