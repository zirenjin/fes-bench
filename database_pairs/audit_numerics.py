"""Claim-evidence audit, numerical part: sensitivity of the +-50 K budget to the 50 K grid and the derivative rule,
and how many transition pairs involve harmonic-invalid entries.
Usage: python audit_numerics.py  -> audit/NUMERICS.json
"""
import json
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicSpline
from score import Bench, GATE

HERE = Path(__file__).parent; (HERE / "audit").mkdir(exist_ok=True)
B = Bench(); T = B.T; out = {}
# reference transition pairs: linear-interpolated root + np.gradient slope (as used) vs cubic-spline root + spline slope
rows = []
for k in np.flatnonzero(B.cross & B.ref_valid):   # reference-valid crossing pairs only (review M2)
    d = B.dref[k]; tc = B.tc[k]; S_lin = -np.interp(tc, T, np.gradient(d, T)) * 1000
    cs = CubicSpline(T, d); roots = [r for r in cs.roots(extrapolate=False) if 50 <= r <= 1500]
    tcs = min(roots, key=lambda r: abs(r - tc)) if roots else np.nan; S_sp = -float(cs(tcs, 1)) * 1000 if roots else np.nan
    # local secant slope over the bracketing grid interval
    i = int(np.searchsorted(T, tc)) - 1; S_sec = -(d[i + 1] - d[i]) / (T[i + 1] - T[i]) * 1000
    rows.append(dict(tc=tc, tc_spline=tcs, S_grad=S_lin, S_spline=S_sp, S_secant=S_sec, edge=bool(tc < 100 or tc > 1450)))
S_g = np.array([r["S_grad"] for r in rows]); S_s = np.array([r["S_spline"] for r in rows]); S_c = np.array([r["S_secant"] for r in rows])
dtc = np.array([abs(r["tc"] - r["tc_spline"]) for r in rows])
rel = np.abs(S_s - S_g) / np.abs(S_g)
out["budget_derivative"] = dict(n_pairs=len(rows), tc_linear_vs_spline_K=dict(median=float(np.nanmedian(dtc)), p90=float(np.nanpercentile(dtc, 90)), max=float(np.nanmax(dtc))),
                                abs_dS_rel_diff_spline_vs_gradient=dict(median=float(np.nanmedian(rel)), p90=float(np.nanpercentile(rel, 90)), max=float(np.nanmax(rel))),
                                abs_dS_rel_diff_secant_vs_gradient=dict(median=float(np.nanmedian(np.abs(S_c - S_g) / np.abs(S_g))), p90=float(np.nanpercentile(np.abs(S_c - S_g) / np.abs(S_g), 90))),
                                eps50_median_meV=dict(gradient=float(np.median(np.abs(S_g)) * 0.05), spline=float(np.nanmedian(np.abs(S_s)) * 0.05), secant=float(np.median(np.abs(S_c)) * 0.05)),
                                n_sign_disagreement_spline=int(np.nansum(np.sign(S_s) != np.sign(S_g))), n_edge_pairs=int(sum(r["edge"] for r in rows)))
# within-budget fractions under each derivative rule, for the three main predictions (S0)
cidx = np.flatnonzero(B.cross & B.ref_valid); frac = {}
for name in ("F0_DPA-3.1-3M_S0", "MatterSim-v1.0.0-5M_S0", "AVG_F0+MatterSim_S0"):
    z = np.load(HERE / f"predictions/{name}.npz"); G = z["E0"][:, None] + z["F"]; d = G[B.pi] - G[B.pj]
    e = np.array([np.interp(B.tc[k], T, d[k]) for k in cidx])
    frac[name] = {rule: float(np.mean(np.abs(e) <= np.abs(S) * 0.05)) for rule, S in (("gradient", S_g), ("spline", S_s), ("secant", S_c))}
out["within_eps50_pair_weighted_by_rule"] = frac
# harmonic validity among transition pairs
fmin = {n: np.load(HERE / f"predictions/{f}.npz")["fmin"] for n, f in (("DFT", "reference_PBE"), ("F0", "F0_DPA-3.1-3M_S0"), ("MS", "MatterSim-v1.0.0-5M_S0"))}
inv = {n: (v[B.pi] < GATE) | (v[B.pj] < GATE) for n, v in fmin.items()}
out["transition_pairs_with_invalid_entry"] = {n: int((v & B.cross).sum()) for n, v in inv.items()} | dict(any=int(((inv["DFT"] | inv["F0"] | inv["MS"]) & B.cross).sum()), of=int(B.cross.sum()))
out["all_pairs_with_invalid_entry"] = {n: int(v.sum()) for n, v in inv.items()} | dict(any=int((inv["DFT"] | inv["F0"] | inv["MS"]).sum()), of=len(B.pi))
(HERE / "audit/NUMERICS.json").write_text(json.dumps(out, indent=1) + "\n"); print(json.dumps(out, indent=1))
