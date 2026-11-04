"""PROTOCOL.md: routes A (location), B (location + uniform frequency scale), C (spectral correction + re-integration),
D (model combination), against 'none', under settings S0 (no reference information) and S1 (static dE0 provided).
Nested, formula-grouped validation: outer 5 folds x 3 repeats, inner 4 folds. Runs on GenSi (sklearn).
Usage: routes.py <FINAL_GROUPS.csv> <curves_F0.npz> <curves_MS.npz> <spectra_dir> <out_dir>
"""
import csv, itertools, json, sys, warnings
from pathlib import Path
import numpy as np
from joblib import Parallel, delayed
from sklearn.linear_model import Ridge
from sklearn.ensemble import GradientBoostingRegressor
warnings.filterwarnings("ignore")

groups_csv, cF0, cMS, spec_dir, out = Path(sys.argv[1]), sys.argv[2], sys.argv[3], Path(sys.argv[4]), Path(sys.argv[5]); out.mkdir(parents=True, exist_ok=True)
MODELS, SETTINGS, ROUTES = ("F0", "MS"), ("S0", "S1"), ("none", "A", "B", "C")
REPEATS, KO, KI, NB = ((1,) if __import__("os").environ.get("ROUTES_SMOKE") else (1, 2, 3)), 5, 4, 2000
CFGS = [("ridge", a) for a in (0.1, 1.0, 10.0, 100.0)] + [("gbm", n) for n in (100, 300)]
KB, HMEV = 0.08617333262, 4.135667696  # meV/K, meV/THz
KNOTS = np.log([0.5, 2.0, 6.0, 20.0])

# ---------------- data ----------------
groups = {r["formula"]: r["mp_ids"].split(";") for r in csv.DictReader(open(groups_csv))}
ids = [m for ms in groups.values() for m in ms]; idx = {m: k for k, m in enumerate(ids)}; NE = len(ids)
form = np.array([f for f, ms in groups.items() for _ in ms]); forms = sorted(groups)
pairs = np.array([(idx[i], idx[j]) for ms in groups.values() for i, j in itertools.combinations(ms, 2)]); pform = form[pairs[:, 0]]
cur = {}
for name, p in (("F0", cF0), ("MS", cMS)):
    z = np.load(p); k = {m: n for n, m in enumerate(z["ids"])}; o = [k[m] for m in ids]
    cur[name] = dict(F=z["F_mlip"][o], E=z["E_mlip"][o]); cur["DFT"] = dict(F=z["F_dft"][o], E=z["E_dft"][o]); T = z["T"]
TS = T >= 50; NT = len(T)
Gref = cur["DFT"]["E"][:, None] + cur["DFT"]["F"]
spec = [dict(np.load(spec_dir / f"{e}.npz")) for e in ids]
BW = {M: np.array([s[f"{M}_w"] for s in spec]) for M in MODELS}; BM = {M: np.array([s[f"{M}_m"] for s in spec]) for M in MODELS}
BD = {M: np.array([s[f"{M}_d"] for s in spec]) for M in MODELS}; BWD = {M: np.array([s[f"{M}_wd"] for s in spec]) for M in MODELS}


def f_mode(lnnu, T_):
    """Per-mode harmonic F (meV), incl. ZPE; lnnu (..., nb), T_ (nt,) -> (..., nt)."""
    hv = HMEV * np.exp(lnnu)[..., None]; kt = KB * np.maximum(T_, 1e-9)
    return hv / 2 + np.where(T_ > 0, kt * np.log1p(-np.exp(-hv / kt)), 0.0)


def F_bins(M, m):
    return np.einsum("eb,ebt->et", BW[M], f_mode(m, T))


# gate (DEVIATIONS.md D1: bins hold every model mode with nu > 0, as phonopy integrates): re-integration reproduces the stored model F at 1000 K within 0.05 meV/atom for >= 99% of entries
# DEVIATIONS.md D2: the gate compares with the same-run phonopy F (binning fidelity); recompute reproducibility against the
# stored F is reported separately, and route C is applied as a difference on top of the stored F.
i1000 = int(np.flatnonzero(T == 1000)[0]); gate = {}; FB0 = {M: F_bins(M, BM[M]) for M in MODELS}
for M in MODELS:
    dev = np.abs(FB0[M][:, i1000] - np.array([s[f"{M}_Fphon"][i1000] for s in spec]))
    rep_dev = np.abs(np.array([s[f"{M}_Fphon"][i1000] for s in spec]) - cur[M]["F"][:, i1000])
    gate[M] = dict(frac_within_0p05=float(np.mean(dev <= 0.05)), max_dev=float(dev.max()), p99_dev=float(np.quantile(dev, 0.99)),
                   recompute_vs_stored=dict(median=float(np.median(rep_dev)), p99=float(np.quantile(rep_dev, 0.99)), max=float(rep_dev.max())))
C_OK = {M: gate[M]["frac_within_0p05"] >= 0.99 for M in MODELS}


def hats(m):
    """Piecewise-linear hat basis on ln(nu), constant beyond end knots; m (...,) -> (..., 4)."""
    x = np.clip(m, KNOTS[0], KNOTS[-1]); B = np.zeros(m.shape + (4,))
    for k in range(4):
        y = np.zeros_like(x)
        if k > 0:
            s = (x >= KNOTS[k - 1]) & (x <= KNOTS[k]); y[s] = (x[s] - KNOTS[k - 1]) / (KNOTS[k] - KNOTS[k - 1])
        if k < 3:
            s = (x >= KNOTS[k]) & (x <= KNOTS[k + 1]); y[s] = (KNOTS[k + 1] - x[s]) / (KNOTS[k + 1] - KNOTS[k])
        if k == 0: y[x <= KNOTS[0]] = 1.0
        if k == 3: y[x >= KNOTS[3]] = 1.0
        B[..., k] = y
    return B


HB = {M: hats(BM[M]) for M in MODELS}


def centred(v, fset=None):
    v = np.asarray(v, float).copy(); f = form
    for g in np.unique(f):
        s = f == g; v[s] -= v[s].mean(axis=0)
    return v


def setting_G(M, S):
    E = cur[M]["E"] if S == "S0" else cur["DFT"]["E"]
    return E[:, None] + cur[M]["F"], E


def features(M, S):
    G, E = setting_G(M, S); F = cur[M]["F"]; it = {t: int(np.flatnonzero(T == t)[0]) for t in (0, 300, 600, 1000, 1500)}
    Sx = -np.gradient(F, T, axis=1)
    mo = np.column_stack([E] + [F[:, it[t]] for t in (0, 300, 600, 1000, 1500)] + [Sx[:, it[300]], Sx[:, it[1000]]] + [np.array([s[f"{M}_spec"] for s in spec])])
    st = np.array([s["struct"] for s in spec]); co = np.array([s["comp"] for s in spec])
    X = np.column_stack([mo, st, co, centred(np.column_stack([mo, st]))])
    return np.nan_to_num(X)


U = {M: cur[M]["F"] - T * np.gradient(cur[M]["F"], T, axis=1) for M in MODELS}
dU = {M: np.gradient(U[M], T, axis=1) for M in MODELS}
XF = {(M, S): features(M, S) for M in MODELS for S in SETTINGS}
GS = {(M, S): setting_G(M, S) for M in MODELS for S in SETTINGS}

# reference crossings (first in 50-1500 K) for metric 1
def crossings(d):
    s = np.sign(d[TS]); Tt = T[TS]; dd = d[TS]; k = np.flatnonzero(s[:-1] * s[1:] < 0)
    return [float(Tt[i] - dd[i] * (Tt[i + 1] - Tt[i]) / (dd[i + 1] - dd[i])) for i in k]


dref = Gref[pairs[:, 0]] - Gref[pairs[:, 1]]
tcref = np.array([(crossings(d) or [np.nan])[0] for d in dref]); cross = np.isfinite(tcref)
Sref_tc = np.array([-np.interp(t, T, np.gradient(d, T)) if np.isfinite(t) else np.nan for t, d in zip(tcref, dref)])
eps50 = np.abs(Sref_tc) * 50


# ---------------- learners ----------------
def fit_predict(cfg, Xtr, ytr, Xte_list):
    mu, sd = Xtr.mean(0), Xtr.std(0); sd[sd == 0] = 1
    Z = lambda X: (X - mu) / sd
    y2 = ytr.reshape(len(ytr), -1)
    if cfg[0] == "ridge":
        r = Ridge(alpha=cfg[1]).fit(Z(Xtr), y2); return [r.predict(Z(X)).reshape((len(X),) + ytr.shape[1:]) for X in Xte_list]
    outs = [np.zeros((len(X), y2.shape[1])) for X in Xte_list]
    for c in range(y2.shape[1]):
        g = GradientBoostingRegressor(max_depth=2, learning_rate=0.05, subsample=0.8, n_estimators=cfg[1], random_state=0).fit(Z(Xtr), y2[:, c])
        for o, X in zip(outs, Xte_list): o[:, c] = g.predict(Z(X))
    return [o.reshape((len(o),) + ytr.shape[1:]) for o in outs]


def folds(fs, k, seed):
    rng = np.random.default_rng(seed); p = list(rng.permutation(fs)); return [set(x) for x in np.array_split(p, k)]


def ents_of(fset):
    return np.flatnonzero(np.isin(form, list(fset)))


def pairs_of(fset):
    return np.flatnonzero(np.isin(pform, list(fset)))


def pair_loss(Gc, ps):
    """Mean over pairs and T in 50-1500 of the squared pair error; Gc is a full (NE, NT) array (only used rows matter)."""
    e = (Gc[pairs[ps, 0]] - Gc[pairs[ps, 1]]) - dref[ps]; return float(np.mean(e[:, TS] ** 2))


def centre_rows(v, rows):
    v = v.copy()
    for g in np.unique(form[rows]):
        s = rows[form[rows] == g]; v[s] -= v[s].mean(axis=0)
    return v


# ---------------- routes ----------------
def stage1(route, M, S, rows):
    """Stage-1 targets (absolute, per phase) for the given rows."""
    G, E = GS[(M, S)]; r = Gref - G
    if route == "B":
        s = -np.gradient(Gref, T, axis=1) + np.gradient(G, T, axis=1)   # S_ref - S_used
        y = np.zeros((len(rows), 1))
        for n, e in enumerate(rows):
            A = np.vstack([np.column_stack([np.ones(TS.sum()), U[M][e, TS]]), np.column_stack([np.zeros(TS.sum()), -100 * dU[M][e, TS]])])
            b = np.r_[r[e, TS], 100 * s[e, TS]]; y[n, 0] = np.linalg.lstsq(A, b, rcond=None)[0][1]
        return y
    if route == "C":
        y = np.zeros((len(rows), 4))
        for n, e in enumerate(rows):
            w = BWD[M][e]; B = HB[M][e]; A = B.T @ (w[:, None] * B) + 1e-3 * np.eye(4); y[n] = np.linalg.solve(A, B.T @ (w * BD[M][e]))
        return y
    return None


def apply1(route, M, S, rows, pred, Gfull):
    G, E = GS[(M, S)]; Gn = Gfull.copy()
    if route == "B":
        Gn[rows] = G[rows] + pred[:, :1] * U[M][rows]
    elif route == "C":
        mp = BM[M][rows] + np.einsum("ebk,ek->eb", HB[M][rows], pred); Gn[rows] = G[rows] + F_bins_rows(M, rows, mp) - FB0[M][rows]
    else:
        Gn[rows] = G[rows]
    return Gn


def F_bins_rows(M, rows, m):
    return np.einsum("eb,ebt->et", BW[M][rows], f_mode(m, T))


def r2(y, p):
    y = y.reshape(len(y), -1); p = p.reshape(len(p), -1)
    return float(np.mean([1 - np.sum((y[:, c] - p[:, c]) ** 2) / max(np.sum((y[:, c] - y[:, c].mean()) ** 2), 1e-12) for c in range(y.shape[1])]))


def run_route(route, M, S, fset_tr, fset_te, seed):
    """Returns corrected G for test rows, inner-OOF corrected G for train rows, diagnostics."""
    G, _ = GS[(M, S)]; X = XF[(M, S)]; tr, te = ents_of(fset_tr), ents_of(fset_te)
    Gout = G.copy(); Goof = G.copy(); diag = {}
    if route == "none":
        return Gout, Goof, diag
    inner = folds(sorted(fset_tr), KI, seed * 100 + 7)
    # stage 1 (B, C)
    G1 = G.copy()
    if route in ("B", "C"):
        y = np.zeros((NE,) + ((1,) if route == "B" else (4,))); y[tr] = stage1(route, M, S, tr)
        best, oof_best = None, None
        for cfg in CFGS:
            oof = np.zeros_like(y)
            for fk in inner:
                vi = ents_of(fk); ti = np.setdiff1d(tr, vi); oof[vi] = fit_predict(cfg, X[ti], y[ti], [X[vi]])[0]
            l = float(np.mean((oof[tr] - y[tr]) ** 2))
            if best is None or l < best[0]: best, oof_best = (l, cfg), oof
        cfg1 = best[1]; ptr_in, pte = fit_predict(cfg1, X[tr], y[tr], [X[tr], X[te]])
        y[te] = stage1(route, M, S, te)
        G1 = apply1(route, M, S, tr, oof_best[tr], G1)      # cross-fitted on training formulas
        G1 = apply1(route, M, S, te, pte, G1)
        diag.update(stage1_cfg=list(cfg1), stage1_r2_train=r2(y[tr], ptr_in), stage1_r2_test=r2(y[te], pte), _s1_test=(te, y[te], pte))
    # stage 2 (location constant, formula-centred)
    a = np.zeros(NE); a[tr] = centre_rows(np.mean((Gref - G1)[:, TS], axis=1), tr)[tr]
    best, oof_best = None, None
    for cfg in CFGS:
        oof = np.zeros(NE)
        for fk in inner:
            vi = ents_of(fk); ti = np.setdiff1d(tr, vi); oof[vi] = fit_predict(cfg, X[ti], a[ti], [X[vi]])[0]
        oof = centre_rows(oof, tr); Gc = G1.copy(); Gc[tr] += oof[tr, None]
        l = pair_loss(Gc, pairs_of(fset_tr))
        if best is None or l < best[0]: best, oof_best = (l, cfg), oof
    cfg2 = best[1]; atr_in, ate = fit_predict(cfg2, X[tr], a[tr], [X[tr], X[te]])
    full = np.zeros(NE); full[te] = ate; full = centre_rows(full, te); full_tr = np.zeros(NE); full_tr[tr] = atr_in; full_tr = centre_rows(full_tr, tr)
    Gout = G1.copy(); Gout[te] += full[te, None]
    Goof = G1.copy(); Goof[tr] += oof_best[tr, None]
    Gin = G1.copy(); Gin[tr] += full_tr[tr, None]
    a_te = np.zeros(NE); a_te[te] = centre_rows(np.mean((Gref - G1)[:, TS], axis=1), te)[te]
    diag.update(stage2_cfg=list(cfg2), stage2_r2_train=r2(a[tr], full_tr[tr]), stage2_r2_test=r2(a_te[te], full[te]),
                pairloss_train=pair_loss(Gin, pairs_of(fset_tr)), pairloss_test=pair_loss(Gout, pairs_of(fset_te)), _s2_test=(te, a_te[te], full[te]))
    return Gout, Goof, diag


def zform(S):
    """Formula-level features for D3."""
    GF, _ = GS[("F0", S)]; GM, _ = GS[("MS", S)]; it = [int(np.flatnonzero(T == t)[0]) for t in (300, 1000, 1500)]
    dis = np.abs((GF[pairs[:, 0]] - GF[pairs[:, 1]]) - (GM[pairs[:, 0]] - GM[pairs[:, 1]]))[:, it]
    ln = np.abs(np.array([s["F0_spec"][0] - s["MS_spec"][0] for s in spec])); co = np.array([s["comp"] for s in spec])
    Z = {}
    for f in forms:
        ps = pform == f; es = form == f; Z[f] = np.r_[dis[ps].mean(0), ln[es].mean(), co[es][0]]
    return Z


ZF = {S: zform(S) for S in SETTINGS}


def d_combine(wf, S):
    GF, _ = GS[("F0", S)]; GM, _ = GS[("MS", S)]; w = np.array([wf[f] for f in form])[:, None]
    return w * GF + (1 - w) * GM   # one weight per formula, so pair differences are w*dF0 + (1-w)*dMS


def fit_d3(S, ftr, lam):
    from scipy.optimize import minimize
    Zt = np.array([ZF[S][f] for f in ftr]); mu, sd = Zt.mean(0), Zt.std(0); sd[sd == 0] = 1
    GF, _ = GS[("F0", S)]; GM, _ = GS[("MS", S)]
    ps = pairs_of(set(ftr)); dF = (GF[pairs[ps, 0]] - GF[pairs[ps, 1]])[:, TS]; dM = (GM[pairs[ps, 0]] - GM[pairs[ps, 1]])[:, TS]; dr = dref[ps][:, TS]
    fi = {f: k for k, f in enumerate(ftr)}; pk = np.array([fi[f] for f in pform[ps]])
    def obj(v):
        w = 1 / (1 + np.exp(-(v[0] + ((Zt - mu) / sd) @ v[1:]))); wp = w[pk][:, None]
        e = wp * dF + (1 - wp) * dM - dr; return float(np.mean(e ** 2) + lam * np.sum(v[1:] ** 2))
    v = minimize(obj, np.zeros(Zt.shape[1] + 1), method="L-BFGS-B").x
    return lambda fs: {f: float(1 / (1 + np.exp(-(v[0] + ((ZF[S][f] - mu) / sd) @ v[1:])))) for f in fs}


def run_d(variant, S, fset_tr, fset_te, seed):
    ftr, fte = sorted(fset_tr), sorted(fset_te); inner = folds(ftr, KI, seed * 100 + 7)
    if variant == "D1":
        wfn = lambda par, fit_f: (lambda fs: {f: 0.5 for f in fs}); grid = [None]
    elif variant == "D2":
        wfn = lambda par, fit_f: (lambda fs: {f: par for f in fs}); grid = list(np.round(np.linspace(0, 1, 11), 2))
    else:
        wfn = lambda par, fit_f: fit_d3(S, fit_f, par); grid = [0.1, 1.0, 10.0]
    best = None
    for par in grid:
        Goof = np.zeros((NE, NT)); wo = {}
        for fk in inner:
            w = wfn(par, [f for f in ftr if f not in fk])(sorted(fk)); wo.update(w)
        Goof = d_combine({**{f: 0.5 for f in forms}, **wo}, S); l = pair_loss(Goof, pairs_of(fset_tr))
        if best is None or l < best[0]: best = (l, par, Goof)
    w = wfn(best[1], ftr)(forms); Gout = d_combine(w, S)
    return Gout, best[2], dict(param=best[1], w_test={f: w[f] for f in fte})


# ---------------- one repeat ----------------
def one_repeat(rep, M, S, route):
    fo = folds(forms, KO, 1000 + rep); Gtest = np.zeros((NE, NT)); thr = {}; diags = []
    for k, fte in enumerate(fo):
        ftr = set(forms) - fte
        if route.startswith("D"):
            Gout, Goof, dg = run_d(route, S, ftr, fte, 1000 + rep * 10 + k)
        else:
            Gout, Goof, dg = run_route(route, M, S, ftr, fte, 1000 + rep * 10 + k)
        te = ents_of(fte); Gtest[te] = Gout[te]
        ptr = pairs_of(ftr); e = (Goof[pairs[ptr, 0]] - Goof[pairs[ptr, 1]]) - dref[ptr]
        for t in (600, 1000, 1500):
            i = int(np.flatnonzero(T == t)[0]); nd = np.abs(dref[ptr, i]) >= 1.0
            for f in fte: thr[(f, t)] = float(np.quantile(np.abs(e[nd, i]), 0.95))
        diags.append(dg)
    return dict(rep=rep, M=M, S=S, route=route, G=Gtest, thr=thr, diags=diags)


jobs = [(rep, M, S, r) for rep in REPEATS for S in SETTINGS for M in MODELS for r in ROUTES if (r != "C" or C_OK[M])]
jobs += [(rep, "D", S, d) for rep in REPEATS for S in SETTINGS for d in ("D1", "D2", "D3")]
results = Parallel(n_jobs=min(len(jobs), 48), verbose=0)(delayed(one_repeat)(*j) for j in jobs)
np.savez_compressed(out / "OOF_G.npz", **{f"{r['rep']}_{r['M']}_{r['S']}_{r['route']}": r["G"] for r in results}, T=T, ids=np.array(ids))

# ---------------- metrics ----------------
it = {t: int(np.flatnonzero(T == t)[0]) for t in (300, 600, 1000, 1500)}
rng_b = np.random.default_rng(20261102); BOOTI = rng_b.integers(0, len(forms), (NB, len(forms)))
FI = np.array([forms.index(f) for f in pform])


def per_pair(G):
    d = G[pairs[:, 0]] - G[pairs[:, 1]]; e = d - dref
    etc = np.array([np.interp(t, T, x) if np.isfinite(t) else np.nan for t, x in zip(tcref, d)])
    wrong = {t: np.where(np.abs(dref[:, it[t]]) >= 1.0, (np.sign(d[:, it[t]]) != np.sign(dref[:, it[t]])).astype(float), np.nan) for t in (600, 1000, 1500)}
    dS = -np.gradient(e, T, axis=1) * 1000
    return d, e, etc, wrong, dS


def perf(v):
    """Per-formula mean of the non-nan pair values (nan where a formula has none)."""
    ok = np.isfinite(v); num = np.bincount(FI[ok], v[ok], len(forms)); den = np.bincount(FI[ok], None, len(forms))
    return np.where(den > 0, num / np.where(den > 0, den, 1), np.nan)


def fw(v, draw=None):
    """Formula-weighted mean; draw = None (all formulas) or an (NB, n_formulas) index array -> (NB,)."""
    pf = perf(v)
    return float(np.nanmean(pf)) if draw is None else np.nanmean(pf[draw], axis=1)


def events(d):
    hits = miss = false = 0; dT = []
    for x, r in zip(d, dref):
        rc, pc = crossings(r), crossings(x); used = set()
        for c in rc:
            j = next((j for j, p in enumerate(pc) if j not in used and abs(p - c) <= 100), None)
            if j is None: miss += 1
            else: hits += 1; used.add(j)
        false += len(pc) - len(used)
        if rc and pc: dT.append(min(abs(p - rc[0]) for p in pc))
    return dict(hits=hits, missed=miss, false=false, n_ref=hits + miss, median_abs_dTc=float(np.median(dT)) if dT else None)


def metrics(R, base=None):
    d, e, etc, wrong, dS = per_pair(R["G"])
    m = dict(M1_rms_e_tc=float(np.sqrt(fw(np.where(cross, etc ** 2, np.nan)))), M1_within_eps50=fw(np.where(cross, (np.abs(etc) <= eps50).astype(float), np.nan)))
    m["M2_wrong"] = {t: fw(wrong[t]) for t in (600, 1000, 1500)}; m["M2_wrong_mean"] = float(np.mean(list(m["M2_wrong"].values())))
    m["M3_pair_rms"] = {t: float(np.sqrt(np.mean(e[:, it[t]] ** 2))) for t in (300, 600, 1000, 1500)}
    m["M4_dS_rms"] = {t: float(np.sqrt(np.mean(dS[:, it[t]] ** 2))) for t in (600, 1000, 1500)}
    m["M5_events"] = events(d)
    m["M6_p95_e1000"] = float(np.quantile(np.abs(e[:, it[1000]]), 0.95)); m["M6_n_max_gt20"] = int(np.sum(np.abs(e[:, TS]).max(1) > 20))
    cov = {}
    for t in (600, 1000, 1500):
        nd = np.abs(dref[:, it[t]]) >= 1.0; th = np.array([R["thr"][(f, t)] for f in pform]); dec = nd & (np.abs(d[:, it[t]]) > th)
        cov[t] = dict(coverage=float(dec.sum() / nd.sum()), selective_risk=float(np.sum(dec & (wrong[t] == 1)) / max(dec.sum(), 1)))
    m["M7_conformal"] = cov
    if base is not None:
        _, eb, etcb, wrongb, _ = per_pair(base["G"])
        worse = (np.abs(e[:, [it[t] for t in (300, 600, 1000, 1500)]]) - np.abs(eb[:, [it[t] for t in (300, 600, 1000, 1500)]])).max(1) > 10
        m["M6_n_worsened_gt10"] = int(worse.sum())
        def bootdiff(fn):
            v = fn(BOOTI); return dict(value=float(fn(None)), ci95=[float(x) for x in np.nanpercentile(v, [2.5, 97.5])])
        m["vs_none"] = dict(
            M1_rms_diff=bootdiff(lambda dr: np.sqrt(fw(np.where(cross, etc ** 2, np.nan), dr)) - np.sqrt(fw(np.where(cross, etcb ** 2, np.nan), dr))),
            M2_wrong_mean_diff_points=bootdiff(lambda dr: 100 * np.mean([fw(wrong[t], dr) - fw(wrongb[t], dr) for t in (600, 1000, 1500)], axis=0)),
            M3_rms1000_diff=bootdiff(lambda dr: np.sqrt(fw(e[:, it[1000]] ** 2, dr)) - np.sqrt(fw(eb[:, it[1000]] ** 2, dr))))
    return m


res = dict(protocol_sha256=open(Path(__file__).parents[1] / "PROTOCOL.sha256").read().split()[0], gate=gate, C_runnable=C_OK,
           n=dict(entries=NE, pairs=len(pairs), formulas=len(forms), crossing_pairs=int(cross.sum())), runs={})
byk = {(r["rep"], r["M"], r["S"], r["route"]): r for r in results}
for (rep, M, S, route), R in byk.items():
    bases = [byk[(rep, b, S, "none")] for b in (MODELS if M == "D" else (M,))]
    key = f"rep{rep}/{M}/{S}/{route}"
    if route == "none":
        res["runs"][key] = metrics(R)
    else:
        res["runs"][key] = {f"vs_{b['M']}": metrics(R, b) for b in bases}
    dg = R["diags"]
    if route in ("A", "B", "C"):
        def pooled(name):
            te = np.concatenate([x[name][0] for x in dg]); y = np.concatenate([x[name][1] for x in dg]); p = np.concatenate([x[name][2] for x in dg]); return r2(y, p)
        res["runs"][key]["M8"] = dict(stage2_r2_train=float(np.mean([x["stage2_r2_train"] for x in dg])), stage2_r2_test_pooled=pooled("_s2_test"),
                                      pairloss_train=float(np.mean([x["pairloss_train"] for x in dg])), pairloss_test=float(np.mean([x["pairloss_test"] for x in dg])),
                                      stage2_cfgs=[x["stage2_cfg"] for x in dg])
        if route in ("B", "C"):
            res["runs"][key]["M8"].update(stage1_r2_train=float(np.mean([x["stage1_r2_train"] for x in dg])), stage1_r2_test_pooled=pooled("_s1_test"),
                                          stage1_cfgs=[x["stage1_cfg"] for x in dg])
    if route.startswith("D"):
        res["runs"][key]["D_params"] = [x["param"] for x in dg]


def verdict(m, n_pairs):
    v = m["vs_none"]; c1 = v["M1_rms_diff"]["ci95"][1] < 0; c2 = v["M2_wrong_mean_diff_points"]["value"] <= 0
    c3 = m["M6_n_worsened_gt10"] <= 0.02 * n_pairs
    return dict(M1_improved=c1, M2_not_worse=c2, tail_ok_worsened=c3, label="success" if (c1 and c2 and c3) else ("partial" if (c1 or c2 or c3) else "failure"))


for key, blk in res["runs"].items():
    rep, M, S, route = key.split("/")
    if route == "none": continue
    for b, m in blk.items():
        if not b.startswith("vs_"): continue
        base = res["runs"][f"{rep}/{b[3:]}/{S}/none"]
        vd = verdict(m, len(pairs)); vd["tail_ok_p95"] = m["M6_p95_e1000"] <= base["M6_p95_e1000"]
        if not vd["tail_ok_p95"] and vd["label"] == "success": vd["label"] = "partial"
        m["verdict"] = vd
(out / "ROUTES.json").write_text(json.dumps(res, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o)) + "\n")
print(json.dumps(dict(gate=gate, C_runnable=C_OK), indent=1))
for key, blk in res["runs"].items():
    if not key.startswith("rep1/"): continue
    if key.endswith("/none"):
        print(f"{key:18s} M1 {blk['M1_rms_e_tc']:.2f} eps50 {blk['M1_within_eps50']:.2f} wrong {100*blk['M2_wrong_mean']:.1f} rms1000 {blk['M3_pair_rms'][1000]:.2f} p95 {blk['M6_p95_e1000']:.1f} n>20 {blk['M6_n_max_gt20']} ev {blk['M5_events']}")
        continue
    for b, m in blk.items():
        if not b.startswith("vs_"): continue
        v = m["vs_none"]
        print(f"{key:18s} {b:6s} M1 {m['M1_rms_e_tc']:.2f} (d {v['M1_rms_diff']['value']:+.2f} [{v['M1_rms_diff']['ci95'][0]:+.2f},{v['M1_rms_diff']['ci95'][1]:+.2f}]) eps50 {m['M1_within_eps50']:.2f} "
              f"wrong {100*m['M2_wrong_mean']:.1f} (d {v['M2_wrong_mean_diff_points']['value']:+.1f}) rms1000 {m['M3_pair_rms'][1000]:.2f} dS1000 {m['M4_dS_rms'][1000]:.1f} "
              f"p95 {m['M6_p95_e1000']:.1f} worse>10 {m['M6_n_worsened_gt10']} n>20 {m['M6_n_max_gt20']} -> {m['verdict']['label']}")
