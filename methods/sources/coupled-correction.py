"""Exploratory (post-hoc, not in PROTOCOL.md): can a coupled anchor/thermal correction beat partial corrections?

Same pairs, same errors and same non-degenerate ordering metric as analyze_covariance.py Q3. Arms per temperature:
  M0  no correction                                  error = dE0 + dZPE + b
  M1s static DFT dE0 on both phases (cheap budget)   error = dZPE + b
  M1  dE0 + dZPE corrected (needs DFT phonons)       error = b
  M2  thermal corrected only                         error = dE0 + dZPE
  M3s M1s budget + coupling learned on other formulas: remainder (dZPE + b) predicted from the known dE0 error
  M3z zero target labels: total error predicted from MLIP-only pair features (dZPE_mlip, dFth_mlip(T)), learned on other formulas
  M3sz M1s budget + MLIP-only features for the remainder
Couplings are leave-one-formula-out least squares through the origin (pair orientation-invariant).
Usage: coupled_correction_check.py <FINAL_GROUPS.csv> <mlip_dir> <out.json>
"""
import csv, itertools, json, sys
from pathlib import Path
import numpy as np

groups_csv, mlip_dir, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
TS = [300, 600, 1000, 1500]; rng = np.random.default_rng(20261101)
groups = {r["formula"]: r["mp_ids"].split(";") for r in csv.DictReader(open(groups_csv))}
E = {}
for f, ms in groups.items():
    for m in ms:
        p = mlip_dir / f"{m}.npz"
        if not p.exists():
            continue
        z = np.load(p); T = z["T"]; ix = {t: int(np.flatnonzero(T == t)[0]) for t in TS}
        Fd, Fm = z["F_dft"], z["F_mlip"]; ed, em = float(z["E_dft_eV_atom"]) * 1000, float(z["E_mlip_eV_atom"]) * 1000
        E[m] = dict(dE0=em - ed, dZPE=float(Fm[0] - Fd[0]), zpe_m=float(Fm[0]),
                    fth_m={t: float(Fm[ix[t]] - Fm[0]) for t in TS}, dFth={t: float((Fm[ix[t]] - Fm[0]) - (Fd[ix[t]] - Fd[0])) for t in TS},
                    Fref={t: ed + float(Fd[ix[t]]) for t in TS})
P = []
for f, ms in groups.items():
    for i, j in itertools.combinations([m for m in ms if m in E], 2):
        d = lambda k: E[i][k] - E[j][k]
        P.append(dict(formula=f, dE0=d("dE0"), dZPE=d("dZPE"), zpe_m=d("zpe_m"),
                      **{f"b{t}": E[i]["dFth"][t] - E[j]["dFth"][t] for t in TS},
                      **{f"fthm{t}": E[i]["fth_m"][t] - E[j]["fth_m"][t] for t in TS},
                      **{f"ref{t}": E[i]["Fref"][t] - E[j]["Fref"][t] for t in TS}))
F = np.array([p["formula"] for p in P]); forms = sorted(set(F))
col = lambda k: np.array([p[k] for p in P])


def lofo(X, y):
    """Leave-one-formula-out prediction of y from X (no intercept)."""
    yhat = np.zeros_like(y)
    for f in forms:
        tr, te = F != f, F == f
        k, *_ = np.linalg.lstsq(X[tr], y[tr], rcond=None); yhat[te] = X[te] @ k
    k_all, *_ = np.linalg.lstsq(X, y, rcond=None)
    return yhat, k_all.tolist()


def fw(v, mask):
    vals = [v[(F == f) & mask].mean() for f in forms if ((F == f) & mask).any()]
    return float(np.mean(vals))


def fw_boot_diff(v1, v2, mask):
    fs = [f for f in forms if ((F == f) & mask).any()]
    per = {f: (v1[(F == f) & mask].mean(), v2[(F == f) & mask].mean()) for f in fs}
    draws = []
    for _ in range(2000):
        s = rng.choice(fs, len(fs)); draws.append(np.mean([per[f][0] - per[f][1] for f in s]))
    return [float(x) for x in np.percentile(draws, [2.5, 97.5])]


res = {"note": "exploratory post-hoc check; not a protocol analysis", "n_pairs": len(P), "n_formulas": len(forms)}
dE0, dZPE, zpe_m = col("dE0"), col("dZPE"), col("zpe_m")
for t in TS:
    b, ref, fthm = col(f"b{t}"), col(f"ref{t}"), col(f"fthm{t}")
    a = dE0 + dZPE; tot = a + b; nd = np.abs(ref) >= 1.0
    rem = dZPE + b
    rem_hat_s, k_s = lofo(dE0[:, None], rem)
    tot_hat_z, k_z = lofo(np.c_[zpe_m, fthm], tot)
    rem_hat_sz, k_sz = lofo(np.c_[zpe_m, fthm], rem)
    err = dict(M0=tot, M1s=rem, M1=b, M2=a, M3s=rem - rem_hat_s, M3z=tot - tot_hat_z, M3sz=rem - rem_hat_sz)
    wrong = {k: (np.sign(ref + e) != np.sign(ref)).astype(float) for k, e in err.items()}
    res[f"T{t}"] = dict(
        cross_term_share=float(2 * np.sum(a * b) / np.sum(tot ** 2)),           # 2<ab>/<(a+b)^2>; <0 means cancellation
        rms_meV={k: float(np.sqrt(np.mean(e ** 2))) for k, e in err.items()},
        wrong_ordering_fw={k: fw(w, nd) for k, w in wrong.items()},
        diff_ci95={"M3s-M1s": fw_boot_diff(wrong["M3s"], wrong["M1s"], nd), "M3z-M0": fw_boot_diff(wrong["M3z"], wrong["M0"], nd),
                   "M3sz-M1s": fw_boot_diff(wrong["M3sz"], wrong["M1s"], nd), "M1s-M1": fw_boot_diff(wrong["M1s"], wrong["M1"], nd)},
        coupling_full_fit=dict(M3s_k_dE0=k_s, M3z_k_zpe_fth=k_z, M3sz_k_zpe_fth=k_sz), n_nondegenerate=int(nd.sum()))
out.write_text(json.dumps(res, indent=1) + "\n"); print(json.dumps(res, indent=1))
