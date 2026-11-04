"""PROTOCOL_RE1_MECH.md scoring: pre-registered directional prediction rho(G, P1) - rho(G, P2) > 0, plus descriptive D3,
identifiability and representation-bias summaries.
G  = |e_M1(1000 K)| - |e_M2(1000 K)| at 25% (RE1 out-of-fold predictions, S1), averaged over repeats;
P1 = |D3_0,pair(1000 K)| - |D3_MLIP,pair(1000 K)|, averaged over repeats;  P2 = -(mean forceRMS of the two phases).
Usage: mech_score.py <data_mech_dir> <data_re1_dir> <RE1.json-unused> <out.json>
"""
import json, sys
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr
V01 = Path("C:/users/ziren/documents/dpa-adapt/outputs/fes-bench-v0.1-20261102"); sys.path.insert(0, str(V01)); sys.path.insert(0, str(Path(__file__).parent))
from score import Bench, GATE  # noqa: E402

dm, dr, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[4]); B = Bench(); T = B.T; i1000 = B.it[1000]; NE = len(B.ids)
ref = np.load(V01 / "predictions/reference_PBE.npz"); E0 = ref["E0"]; Fref = ref["F"]
Zm = {m: dict(np.load(dm / f"{m}.npz")) for m in B.ids}
# reuse the RE1 out-of-fold assembly (same selection rule) for M1/M2 at the primary budget
import importlib.util
spec = importlib.util.spec_from_file_location("re1s", Path(__file__).parent / "re1_score_lib.py")
REP = (1, 2); PRIMARY = 0.25; rng = np.random.default_rng(20261104)
fth = lambda F: F[..., i1000] - F[..., 0]


_argv = sys.argv; sys.argv = ["x", str(dr)]; LIB = importlib.util.module_from_spec(spec); spec.loader.exec_module(LIB); sys.argv = _argv


def oof(method, M):
    """The RE1 selection for (method, M, 25%) (identical code), per-repeat (F, fmin)."""
    Fd, fd, _ = LIB.assemble(method, M, PRIMARY); return Fd, fd


res = dict(protocol_sha256=open(Path(__file__).parents[1] / "PROTOCOL_RE1_MECH.sha256").read().split()[0])
app = np.array([bool(Zm[m]["applicable"]) for m in B.ids]); res["n_entries_applicable"] = int(app.sum()); res["n_entries_npar_gt2000_or_none"] = int((~app).sum())
for M in ("F0", "MS"):
    F1, f1 = oof("M1", M); F2, f2 = oof("M2", M)
    Gs, okG = [], np.ones(len(B.P), bool)
    for r in REP:
        G1 = E0[:, None] + F1[r]; G2 = E0[:, None] + F2[r]
        e1 = np.abs((G1[B.pi, i1000] - G1[B.pj, i1000]) - B.dref[:, i1000]); e2 = np.abs((G2[B.pi, i1000] - G2[B.pj, i1000]) - B.dref[:, i1000])
        okG &= (f1[r][B.pi] >= GATE) & (f1[r][B.pj] >= GATE) & (f2[r][B.pi] >= GATE) & (f2[r][B.pj] >= GATE); Gs.append(e1 - e2)
    Gv = np.mean(Gs, axis=0)
    P1s, okP = [], np.ones(len(B.P), bool)
    for r in REP:
        def pair_d3(tag):
            v = np.array([fth(Zm[m][f"F_D3_{tag}_{PRIMARY}_{r}"]) - fth(Zm[m]["F_ref_A1"]) if app[n] else np.nan for n, m in enumerate(B.ids)])
            fm = np.array([float(Zm[m][f"fmin_D3_{tag}_{PRIMARY}_{r}"]) if app[n] else -np.inf for n, m in enumerate(B.ids)]); return v[B.pi] - v[B.pj], fm
        d0, fm0 = pair_d3("zero"); dM, fmM = pair_d3(M)
        okP &= np.isfinite(d0) & np.isfinite(dM) & (fm0[B.pi] >= GATE) & (fm0[B.pj] >= GATE) & (fmM[B.pi] >= GATE) & (fmM[B.pj] >= GATE); P1s.append(np.abs(d0) - np.abs(dM))
    P1 = np.mean(P1s, axis=0); frms = np.array([float(Zm[m][f"forceRMS_{M}"]) for m in B.ids]); P2 = -(frms[B.pi] + frms[B.pj]) / 2
    unit = B.ref_valid & app[B.pi] & app[B.pj] & okG & okP
    # also a looser unit dropping the D3-validity requirement on the zero prior (descriptive)
    idx = np.flatnonzero(unit); fsel = B.pf[idx]
    def rho(x, y): return spearmanr(x, y).correlation if len(x) > 3 else np.nan
    def test(Gx, P1x, P2x, idx_):
        fs_ = B.pf[idx_]; forms_ = np.unique(fs_)
        if len(idx_) < 10: return dict(n_pairs=int(len(idx_)), n_formulas=int(len(forms_)), verdict=f"not evaluable (n = {len(idx_)})")
        r1, r2 = rho(Gx[idx_], P1x[idx_]), rho(Gx[idx_], P2x[idx_]); draws = []
        for _ in range(2000):
            pick = rng.choice(forms_, len(forms_)); ii = np.concatenate([idx_[fs_ == f] for f in pick]); draws.append(rho(Gx[ii], P1x[ii]) - rho(Gx[ii], P2x[ii]))
        lo, hi = np.nanpercentile(draws, [2.5, 97.5])
        return dict(n_pairs=int(len(idx_)), n_formulas=int(len(forms_)), rho_G_P1=float(r1), rho_G_P2=float(r2), diff=float(r1 - r2), diff_ci95=[float(lo), float(hi)],
                    verdict="mechanism prediction supported" if lo > 0 else "not supported", G_median=float(np.median(Gx[idx_])), G_frac_positive=float(np.mean(Gx[idx_] > 0)))
    blk = dict(preregistered=test(Gv, P1, P2, idx), unit_counts=dict(reference_valid=int(B.ref_valid.sum()), A1_applicable_both=int((B.ref_valid & app[B.pi] & app[B.pj]).sum()),
                                                                        plus_M1_M2_valid=int((B.ref_valid & app[B.pi] & app[B.pj] & okG).sum()), plus_D3_valid=int(len(idx))))
    # POST HOC (labelled; not the pre-registered test): failures kept in the denominator at the RE1 selection cap of 20 meV/atom
    Gc, P1c = [], []
    for r in REP:
        G1 = E0[:, None] + F1[r]; G2 = E0[:, None] + F2[r]
        e1 = np.where((f1[r][B.pi] >= GATE) & (f1[r][B.pj] >= GATE), np.minimum(np.abs((G1[B.pi, i1000] - G1[B.pj, i1000]) - B.dref[:, i1000]), 20), 20)
        e2 = np.where((f2[r][B.pi] >= GATE) & (f2[r][B.pj] >= GATE), np.minimum(np.abs((G2[B.pi, i1000] - G2[B.pj, i1000]) - B.dref[:, i1000]), 20), 20)
        Gc.append(e1 - e2)
        def capd3(tag):
            v = np.array([fth(Zm[m][f"F_D3_{tag}_{PRIMARY}_{r}"]) - fth(Zm[m]["F_ref_A1"]) if app[n] else np.nan for n, m in enumerate(B.ids)])
            fm = np.array([float(Zm[m][f"fmin_D3_{tag}_{PRIMARY}_{r}"]) if app[n] else -np.inf for n, m in enumerate(B.ids)])
            pv = np.abs(v[B.pi] - v[B.pj]); ok = (fm[B.pi] >= GATE) & (fm[B.pj] >= GATE); return np.where(ok, np.minimum(pv, 20), 20)
        P1c.append(capd3("zero") - capd3(M))
    unit_ph = B.ref_valid & app[B.pi] & app[B.pj]
    blk["post_hoc_capped"] = test(np.mean(Gc, axis=0), np.mean(P1c, axis=0), P2, np.flatnonzero(unit_ph))
    blk["post_hoc_capped"]["note"] = "post hoc, descriptive: invalid fits and invalid D3 spectra enter at the 20 meV/atom cap; not the pre-registered test"
    # descriptive D3 by budget, identifiability, representation bias
    d3 = {}
    for b in (0.1, 0.25, 0.5, 1.0):
        for tag in ("zero", M):
            vals = [abs(float(fth(Zm[m][f"F_D3_{tag}_{b}_1"]) - fth(Zm[m]["F_ref_A1"]))) for n, m in enumerate(B.ids) if app[n] and float(Zm[m][f"fmin_D3_{tag}_{b}_1"]) >= GATE]
            inval = sum(1 for n, m in enumerate(B.ids) if app[n] and float(Zm[m][f"fmin_D3_{tag}_{b}_1"]) < GATE)
            d3[f"{tag}_{b}"] = dict(median_abs_single_meV=float(np.median(vals)) if vals else None, n_valid=len(vals), n_invalid=inval)
        rk = [float(Zm[m][f"rank_frac_{b}_1"]) for n, m in enumerate(B.ids) if app[n]]; d3[f"rank_frac_{b}"] = dict(median=float(np.median(rk)), p10=float(np.percentile(rk, 10)))
    rep_bias = [abs(float(fth(Zm[m]["F_ref_A1"]) - fth(Fref[n]))) for n, m in enumerate(B.ids) if app[n] and float(Zm[m]["fmin_ref_A1"]) >= GATE]
    blk["descriptive"] = dict(D3=d3, representation_bias_single_abs_meV_median=float(np.median(rep_bias)), representation_bias_invalid=int(sum(1 for n, m in enumerate(B.ids) if app[n] and float(Zm[m]["fmin_ref_A1"]) < GATE)),
                              theta_err_norm_median=float(np.median([float(Zm[m][f"theta_err_norm_{M}"]) for n, m in enumerate(B.ids) if app[n]])))
    res[M] = blk
out.write_text(json.dumps(res, indent=1, default=float) + "\n"); print(json.dumps(res, indent=1, default=float))
