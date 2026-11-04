"""Exploratory (post-hoc, not in PROTOCOL.md): residual structure and held-out-formula abstention for phase ordering.

Part A, residual structure (eps = MLIP - DFT, pair level):
  thermal  eps_th = dFth_MLIP - dFth_DFT at the DFT geometry; total eps = dE0 + dZPE + eps_th.
  Reports corr(eps, D), Spearman(|eps|, |D|), RMS eps by |D| tertile, between-formula variance share, corr of eps_th across T.
Part B, abstention (no Gaussian assumption): for each held-out formula, the threshold is the (1 - alpha) quantile of |eps|
  over all pairs of the other formulas, at the same T ("T-aware") or pooled over the four T ("pooled").
  A pair is decided only if |dF_MLIP| > threshold; a wrong decision requires |eps| > |dF_MLIP|, so P(decided and wrong)
  <= alpha is expected under formula exchangeability. Arms: M0 (raw MLIP) and M1s (static DFT dE0 supplied).
  Ordering metrics use non-degenerate pairs (|dF_DFT| >= 1 meV/atom), as in Q3.
Usage: abstention_check.py <FINAL_GROUPS.csv> <mlip_dir> <out.json>
"""
import csv, itertools, json, sys
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr

groups_csv, mlip_dir, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
TS = [300, 600, 1000, 1500]; ALPHAS = [0.10, 0.05, 0.02]; rng = np.random.default_rng(20261101)
groups = {r["formula"]: r["mp_ids"].split(";") for r in csv.DictReader(open(groups_csv))}
E = {}
for f, ms in groups.items():
    for m in ms:
        p = mlip_dir / f"{m}.npz"
        if not p.exists():
            continue
        z = np.load(p); T = z["T"]; ix = {t: int(np.flatnonzero(T == t)[0]) for t in TS}
        Fd, Fm = z["F_dft"], z["F_mlip"]; ed, em = float(z["E_dft_eV_atom"]) * 1000, float(z["E_mlip_eV_atom"]) * 1000
        E[m] = dict(ed=ed, em=em, zd=float(Fd[0]), zm=float(Fm[0]),
                    thd={t: float(Fd[ix[t]] - Fd[0]) for t in TS}, thm={t: float(Fm[ix[t]] - Fm[0]) for t in TS})
P = []
for f, ms in groups.items():
    for i, j in itertools.combinations([m for m in ms if m in E], 2):
        d = lambda k: E[i][k] - E[j][k]
        row = dict(formula=f, dE0=d("em") - d("ed"), dZPE=d("zm") - d("zd"))
        for t in TS:
            Dth = E[i]["thd"][t] - E[j]["thd"][t]; Mth = E[i]["thm"][t] - E[j]["thm"][t]
            row[f"Dth{t}"] = Dth; row[f"eth{t}"] = Mth - Dth
            row[f"D{t}"] = d("ed") + d("zd") + Dth                     # DFT total pair gap
        P.append(row)
F = np.array([p["formula"] for p in P]); forms = sorted(set(F)); col = lambda k: np.array([p[k] for p in P])
res = {"note": "exploratory post-hoc check; not a protocol analysis", "n_pairs": len(P), "n_formulas": len(forms)}

# Part A
A = {}
for t in TS:
    Dth, eth, D = col(f"Dth{t}"), col(f"eth{t}"), col(f"D{t}")
    eps = col("dE0") + col("dZPE") + eth
    blk = {}
    for name, x, e in (("thermal", Dth, eth), ("total", D, eps)):
        q = np.quantile(np.abs(x), [1 / 3, 2 / 3]); bins = np.digitize(np.abs(x), q)
        fm = np.array([e[F == f].mean() for f in F])
        blk[name] = dict(pearson_eps_D=float(np.corrcoef(np.r_[e, -e], np.r_[x, -x])[0, 1]),   # symmetrized: orientation-invariant
                         spearman_abs_eps_abs_D=float(spearmanr(np.abs(e), np.abs(x)).correlation),
                         rms_eps_by_absD_tertile=[float(np.sqrt(np.mean(e[bins == k] ** 2))) for k in range(3)],
                         absD_tertile_edges=[float(v) for v in q],
                         between_formula_share_of_eps2=float(np.sum(fm ** 2) / np.sum(e ** 2)))
    A[f"T{t}"] = blk
eth_T = np.array([col(f"eth{t}") for t in TS])
A["thermal_eps_corr_across_T"] = {f"{a}-{b}": float(np.corrcoef(eth_T[i], eth_T[j])[0, 1]) for (i, a), (j, b) in itertools.combinations(enumerate(TS), 2)}
res["A_residual_structure"] = A

# Part B
err_by_arm_T = {}
for t in TS:
    eth = col(f"eth{t}")
    err_by_arm_T[("M0", t)] = col("dE0") + col("dZPE") + eth
    err_by_arm_T[("M1s", t)] = col("dZPE") + eth


def cluster_ci(v, mask):
    fs = [f for f in forms if ((F == f) & mask).any()]
    num = {f: v[(F == f) & mask].sum() for f in fs}; den = {f: ((F == f) & mask).sum() for f in fs}
    draws = []
    for _ in range(2000):
        s = rng.choice(fs, len(fs)); d = sum(den[f] for f in s); draws.append(sum(num[f] for f in s) / d if d else np.nan)
    return [float(x) for x in np.nanpercentile(draws, [2.5, 97.5])]


B = {}
for arm in ("M0", "M1s"):
    for t in TS:
        D = col(f"D{t}"); e = err_by_arm_T[(arm, t)]; pred = D + e; nd = np.abs(D) >= 1.0
        wrong = (np.sign(pred) != np.sign(D)) & nd
        blk = dict(no_abstention=dict(wrong_rate=float(wrong[nd].mean()), n=int(nd.sum())))
        for mode in ("T-aware", "pooled"):
            for a in ALPHAS:
                thr = np.zeros(len(P))
                for f in forms:
                    tr = F != f
                    pool = np.abs(e[tr]) if mode == "T-aware" else np.abs(np.concatenate([err_by_arm_T[(arm, s)][tr] for s in TS]))
                    thr[F == f] = np.quantile(pool, 1 - a)
                dec = np.abs(pred) > thr
                acc = dec & nd
                blk[f"{mode}_alpha{a}"] = dict(
                    coverage=float(acc.sum() / nd.sum()), joint_wrong_and_decided=float((acc & wrong).sum() / nd.sum()),
                    selective_risk=float((acc & wrong).sum() / max(acc.sum(), 1)),
                    selective_risk_ci95=cluster_ci((acc & wrong).astype(float), acc) if acc.any() else None,
                    mean_threshold_meV=float(thr.mean()))
        # label-free ranking curve: sort by |pred|; oracle sorts by |D|
        order = np.argsort(-np.abs(pred[nd])); w = wrong[nd][order]
        risk = np.cumsum(w) / np.arange(1, len(w) + 1); cov = np.arange(1, len(w) + 1) / len(w)
        ow = wrong[nd][np.argsort(-np.abs(D[nd]))]; orisk = np.cumsum(ow) / np.arange(1, len(ow) + 1)
        blk["AURC"] = float(risk.mean()); blk["AURC_oracle_true_gap_order"] = float(orisk.mean())
        blk["risk_at_coverage"] = {str(c): float(risk[int(c * len(w)) - 1]) for c in (0.5, 0.7, 0.9, 1.0)}
        B[f"{arm}_T{t}"] = blk
res["B_abstention"] = B
out.write_text(json.dumps(res, indent=1) + "\n")
for k, v in A.items():
    if k.startswith("T"):
        print(k, {n: {kk: (round(x, 3) if isinstance(x, float) else [round(y, 2) for y in x]) for kk, x in b.items()} for n, b in v.items()})
print("eps_th corr across T", {k: round(v, 3) for k, v in A["thermal_eps_corr_across_T"].items()})
for k, b in B.items():
    print(k, "none", round(b["no_abstention"]["wrong_rate"] * 100, 1), "AURC", round(b["AURC"] * 100, 2), "oracle", round(b["AURC_oracle_true_gap_order"] * 100, 2))
    for m in ("T-aware", "pooled"):
        print("   ", m, {a: (round(b[f"{m}_alpha{a}"]["coverage"], 3), round(b[f"{m}_alpha{a}"]["joint_wrong_and_decided"], 4), round(b[f"{m}_alpha{a}"]["selective_risk"], 4)) for a in ALPHAS})
