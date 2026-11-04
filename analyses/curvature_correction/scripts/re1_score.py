"""PROTOCOL_RE1.md scoring. Hyperparameters (cutoff, lambda) are chosen per (method, model, budget) inside the outer
training formulas (v0.1 SPLITS repeat 1, 5 folds) by the frozen objective J, averaged over the two label repeats; the
selected configuration is applied to the test fold. Each repeat is scored separately with the v0.1 scorer (S1) and the
paired formula-cluster bootstrap draws are averaged over repeats.
Usage: re1_score.py <data_re1_dir> <out.json>
"""
import csv, json, sys, tempfile
from pathlib import Path
import numpy as np
V01 = Path("C:/users/ziren/documents/dpa-adapt/outputs/fes-bench-v0.1-20261102"); sys.path.insert(0, str(V01))
from score import Bench, GATE  # noqa: E402

data, out = Path(sys.argv[1]), Path(sys.argv[2]); B = Bench(); T = B.T; NE = len(B.ids); i1000 = B.it[1000]; it4 = [B.it[t] for t in (300, 600, 1000, 1500)]
ref = np.load(V01 / "predictions/reference_PBE.npz"); E0 = ref["E0"]; Z = {m: dict(np.load(data / f"{m}.npz")) for m in B.ids}
BUD = (0.1, 0.25, 0.5, 1.0); REP = (1, 2); LAMS = ("1e-06", "0.001", "0.1", "10"); CUTS = ("rc", "rc4"); PRIMARY = 0.25
fold_of = {r["formula"]: int(r["outer_fold_repeat1"]) for r in csv.DictReader(open(V01 / "manifest/SPLITS.csv"))}
efold = np.array([fold_of[B.E[n]["formula"]] for n in range(NE)]); pfold = efold[B.pi]
tmp = Path(tempfile.mkdtemp())


def arrays(key_F, key_f):
    F = np.array([Z[m].get(key_F, np.full(len(T), np.nan)) for m in B.ids], float); fm = np.array([float(Z[m].get(key_f, -np.inf)) for m in B.ids])
    fm = np.where(np.isfinite(fm), fm, -np.inf); return np.where(np.isfinite(F), F, 0.0), fm


def J(F, fm, fold_mask):
    """Frozen selection objective on training formulas: formula-weighted mean of min(|e(1000 K)|, 20), invalid -> 20."""
    G = E0[:, None] + F; e = np.abs((G[B.pi, i1000] - G[B.pj, i1000]) - B.dref[:, i1000]); ok = (fm[B.pi] >= GATE) & (fm[B.pj] >= GATE)
    v = np.where(ok, np.minimum(e, 20.0), 20.0); m = B.ref_valid & fold_mask
    s, n = B.per_formula(np.where(m, v, np.nan), B.pf); pf = np.where(n > 0, s / np.where(n > 0, n, 1), np.nan); return float(np.nanmean(pf))


def assemble(method, M, b):
    """Out-of-fold assembly: per test fold, the (cut, lam) minimising J on the training folds (averaged over repeats)."""
    tag = (lambda c, l, r: f"M1_{c}_{l}_{b}_{r}") if method == "M1" else (lambda c, l, r: f"{method}_{M}_{c}_{l}_{b}_{r}")
    cache = {(c, l, r): arrays("F_" + tag(c, l, r), "fmin_" + tag(c, l, r)) for c in CUTS for l in LAMS for r in REP}
    outF = {r: np.zeros((NE, len(T))) for r in REP}; outf = {r: np.full(NE, -np.inf) for r in REP}; chosen = {}
    for k in range(5):
        train = pfold != k; best = min(((np.mean([J(*cache[(c, l, r)], train) for r in REP]), c, l) for c in CUTS for l in LAMS))
        chosen[k] = best[1:]; te = efold == k
        for r in REP:
            F, fm = cache[(best[1], best[2], r)]; outF[r][te] = F[te]; outf[r][te] = fm[te]
    return outF, outf, chosen


def save(name, F, fm):
    p = tmp / f"{name}.npz"; np.savez(p, ids=np.array(B.ids), T=T, E0=E0, F=F, fmin=fm, info=json.dumps(dict(name=name, setting="S1"))); return p


res = dict(protocol_sha256=open(Path(__file__).parents[1] / "PROTOCOL_RE1.sha256").read().split()[0],
           proj_equiv_check=dict(median=float(np.median([Z[m]["proj_equiv_check"] for m in B.ids if "proj_equiv_check" in Z[m]])),
                                 max=float(np.max([Z[m]["proj_equiv_check"] for m in B.ids if "proj_equiv_check" in Z[m]]))))
nd = np.array([int(Z[m]["n_disp"]) for m in B.ids]); nat = np.array([int(Z[m]["natom_sc"]) for m in B.ids])
res["reference_cost"] = {str(b): {str(r): dict(configs=int(sum(int(Z[m][f"nconf_{b}_{r}"]) for m in B.ids)), fraction=float(sum(int(Z[m][f"nconf_{b}_{r}"]) for m in B.ids) / nd.sum()),
                                               supercell_atom_configs=int(sum(int(Z[m][f"nconf_{b}_{r}"]) * nat[n] for n, m in enumerate(B.ids)))) for r in REP} for b in BUD}
res["reference_cost"]["full_configs"] = int(nd.sum())
S = {}   # (method, M, b) -> list of per-repeat score dicts (model_valid) and their prediction paths
for M in ("F0", "MS"):
    F0_, f0_ = arrays(f"F_M0_{M}", f"fmin_M0_{M}"); p = save(f"M0_{M}", F0_, f0_); s = B.score(p, "model_valid"); S[("M0", M, None)] = [(p, s), (p, s)]
    for b in BUD:
        for method in ("M1", "M2", "M3"):
            Fd, fd, ch = assemble(method, M, b); S[(method, M, b)] = [(save(f"{method}_{M}_{b}_{r}", Fd[r], fd[r]), None) for r in REP]
            S[(method, M, b)] = [(pth, B.score(pth, "model_valid")) for pth, _ in S[(method, M, b)]]
            res.setdefault("chosen", {})[f"{method}_{M}_{b}"] = {str(k): list(v) for k, v in ch.items()}


def avg_metric(lst, k):
    v = np.mean([s["metrics"][k]["value"] for _, s in lst]); dr = np.mean([s["metrics"][k]["_draws"] for _, s in lst], axis=0); return float(v), dr


def paired(a, b, k="crossing_success50"):
    va, da = avg_metric(a, k); vb, db = avg_metric(b, k); lo, hi = np.nanpercentile(da - db, [2.5, 97.5]); return dict(value=va - vb, ci95=[float(lo), float(hi)])


def errs(a, b):
    """Error metrics on pairs valid for both a and b (per repeat, averaged), plus pairs worsened by > 10 meV/atom vs b."""
    outd = {}; worse = []
    for (pa, _), (pb, _) in zip(a, b):
        za, zb = np.load(pa), np.load(pb); va = B.ref_valid & (za["fmin"][B.pi] >= GATE) & (za["fmin"][B.pj] >= GATE); vb = B.ref_valid & (zb["fmin"][B.pi] >= GATE) & (zb["fmin"][B.pj] >= GATE)
        mask = va & vb; sa = B.score(pa, mask)["metrics"]; sb = B.score(pb, mask)["metrics"]
        for k in ("pair_dFth_rms_1000K", "pair_dS_rms_1000K", "near_Tc_rms_e", "wrong_ordering_mean"):
            outd.setdefault(k, []).append((sa[k]["value"], sb[k]["value"]))
        Ga = E0[:, None] + za["F"]; Gb = E0[:, None] + zb["F"]; ea = (Ga[B.pi] - Ga[B.pj]) - B.dref; eb = (Gb[B.pi] - Gb[B.pj]) - B.dref
        worse.append(int((mask & ((np.abs(ea[:, it4]) - np.abs(eb[:, it4])).max(1) > 10)).sum())); outd.setdefault("n_compared", []).append(int(mask.sum()))
    o = {k: [float(np.mean([x[0] for x in v])), float(np.mean([x[1] for x in v]))] for k, v in outd.items() if k != "n_compared"}
    o["n_compared"] = float(np.mean(outd["n_compared"])); o["n_worsened_gt10"] = float(np.mean(worse)); return o


for M in ("F0", "MS"):
    blk = {}; base = S[("M0", M, None)]
    blk["M0"] = dict(success=avg_metric(base, "crossing_success50")[0], failure=avg_metric(base, "failure_model_invalid")[0])
    for b in BUD:
        bb = {}
        for method in ("M1", "M2", "M3"):
            L = S[(method, M, b)]; bb[method] = dict(success=avg_metric(L, "crossing_success50")[0], failure=avg_metric(L, "failure_model_invalid")[0],
                                                     ordering_success=avg_metric(L, "ordering_success_mean")[0], hit50_valid=avg_metric(L, "crossing_hit50")[0])
        bb["M2_minus_M1"] = paired(S[("M2", M, b)], S[("M1", M, b)]); bb["M2_minus_M0"] = paired(S[("M2", M, b)], base)
        bb["M3_minus_M1"] = paired(S[("M3", M, b)], S[("M1", M, b)]); bb["M2_minus_M3"] = paired(S[("M2", M, b)], S[("M3", M, b)])
        bb["errors_M2_vs_M0"] = errs(S[("M2", M, b)], base); bb["errors_M2_vs_M1"] = errs(S[("M2", M, b)], S[("M1", M, b)])
        blk[str(b)] = bb
    p = blk[str(PRIMARY)]; nbase = int(B.ref_valid.sum())
    tier1 = (p["M2_minus_M1"]["ci95"][0] > 0) and (p["M2_minus_M0"]["ci95"][0] > 0) and (p["errors_M2_vs_M0"]["n_worsened_gt10"] <= 0.02 * nbase)
    tier2 = tier1 and p["M2"]["success"] >= 0.5 and p["errors_M2_vs_M0"]["wrong_ordering_mean"][0] <= p["errors_M2_vs_M0"]["wrong_ordering_mean"][1]
    ni = paired(S[("M2", M, PRIMARY)], S[("M2", M, 1.0)])
    near_same = all((blk[str(b)]["M2_minus_M1"]["ci95"][0] <= 0 <= blk[str(b)]["M2_minus_M1"]["ci95"][1]) for b in BUD) and abs(p["M2_minus_M1"]["value"]) < 0.03
    blk["verdict"] = dict(tier1_prior_effective=bool(tier1), tier2_practical=bool(tier2), noninferiority_25_vs_100=ni, noninferior=bool(ni["ci95"][0] > -0.05),
                          near_identical_M2_M1=bool(near_same))
    res[M] = blk


def strip(o):
    if isinstance(o, dict): return {str(k): strip(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)): return [strip(v) for v in o]
    if isinstance(o, (np.floating, np.integer, np.bool_)): return o.item()
    return o


out.write_text(json.dumps(strip(res), indent=1) + "\n")
f = lambda d: f"{d['value']:+.3f} [{d['ci95'][0]:+.3f},{d['ci95'][1]:+.3f}]"
print("proj_equiv_check", res["proj_equiv_check"]); print("reference cost", json.dumps(res["reference_cost"]))
for M in ("F0", "MS"):
    blk = res[M]; print(f"== {M}: M0 success {blk['M0']['success']:.3f} failure {blk['M0']['failure']:.3f}")
    for b in BUD:
        bb = blk[str(b)]
        print(f"  budget {b:4.2f}: success M1 {bb['M1']['success']:.3f} M2 {bb['M2']['success']:.3f} M3 {bb['M3']['success']:.3f} | M2-M1 {f(bb['M2_minus_M1'])} M2-M0 {f(bb['M2_minus_M0'])} "
              f"M3-M1 {f(bb['M3_minus_M1'])} M2-M3 {f(bb['M2_minus_M3'])} | fail M1/M2 {bb['M1']['failure']:.2f}/{bb['M2']['failure']:.2f} | M2 vs M0 errs {json.dumps({k: [round(x, 3) for x in v] if isinstance(v, list) else v for k, v in bb['errors_M2_vs_M0'].items()})}")
    print("  VERDICT", json.dumps(blk["verdict"]))
