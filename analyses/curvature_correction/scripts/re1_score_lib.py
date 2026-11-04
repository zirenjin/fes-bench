"""Library extract of re1_score.py (selection and out-of-fold assembly only; identical code)."""
import csv, json, sys, tempfile
from pathlib import Path
import numpy as np
V01 = Path("C:/users/ziren/documents/dpa-adapt/outputs/fes-bench-v0.1-20261102"); sys.path.insert(0, str(V01))
from score import Bench, GATE  # noqa: E402

data, out = Path(sys.argv[1]), Path(sys.argv[2]) if len(sys.argv) > 2 else None; B = Bench(); T = B.T; NE = len(B.ids); i1000 = B.it[1000]; it4 = [B.it[t] for t in (300, 600, 1000, 1500)]
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


