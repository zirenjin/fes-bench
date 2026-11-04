"""PROTOCOL_1c_A1.md scoring: A1_proj vs base_proj (gate criteria), A1_proj vs R3_proj (head-to-head), oracle6_proj
(descriptive). Same scorer and S1 setting as gates 1 and 1b.
Usage: gate1c_score.py <data1c_dir> <out.json>
"""
import json, sys, tempfile
from pathlib import Path
import numpy as np
V01 = Path("C:/users/ziren/documents/dpa-adapt/outputs/fes-bench-v0.1-20261102"); sys.path.insert(0, str(V01))
from score import Bench, GATE  # noqa: E402

data, out = Path(sys.argv[1]), Path(sys.argv[2]); B = Bench(); T = B.T; it4 = [B.it[t] for t in (300, 600, 1000, 1500)]
ref = np.load(V01 / "predictions/reference_PBE.npz"); E0 = ref["E0"]; Z = {m: np.load(data / f"{m}.npz") for m in B.ids}
tmp = Path(tempfile.mkdtemp()); res = dict(protocol_sha256=open(Path(__file__).parents[1] / "PROTOCOL_1c_A1.sha256").read().split()[0])
VARS = ("base_proj", "R3_proj", "A1_proj", "A1", "oracle6_proj")


def compare(P, valid, base_v, var_v):
    base = valid[base_v]; mask = base & valid[var_v]; s = B.score(P[var_v][0], mask)["metrics"]; b = B.score(P[base_v][0], mask)["metrics"]
    def diff(k, sc=1.0):
        dr = (s[k]["_draws"] - b[k]["_draws"]) * sc; lo, hi = np.nanpercentile(dr, [2.5, 97.5]); return dict(value=float((s[k]["value"] - b[k]["value"]) * sc), ci95=[float(lo), float(hi)])
    Gv = E0[:, None] + P[var_v][1]; Gb = E0[:, None] + P[base_v][1]
    e_v = (Gv[B.pi] - Gv[B.pj]) - B.dref; e_b = (Gb[B.pi] - Gb[B.pj]) - B.dref
    worse = mask & ((np.abs(e_v[:, it4]) - np.abs(e_b[:, it4])).max(1) > 10)
    fv = B.score(P[var_v][0], "model_valid")["metrics"]; fb = B.score(P[base_v][0], "model_valid")["metrics"]
    dsucc = (fv["crossing_success50"]["_draws"] - fb["crossing_success50"]["_draws"])
    r = dict(n_base=int(base.sum()), n_compared=int(mask.sum()), new_failures=int((base & ~valid[var_v]).sum()), newly_valid=int((B.ref_valid & ~base & valid[var_v]).sum()),
             var_metrics={k: s[k]["value"] for k in ("pair_dFth_rms_1000K", "pair_dS_rms_1000K", "near_Tc_rms_e", "wrong_ordering_mean", "crossing_hit50")},
             base_metrics={k: b[k]["value"] for k in ("pair_dFth_rms_1000K", "pair_dS_rms_1000K", "near_Tc_rms_e", "wrong_ordering_mean", "crossing_hit50")},
             d_pair_dFth=diff("pair_dFth_rms_1000K"), d_pair_dS=diff("pair_dS_rms_1000K"), d_near_Tc=diff("near_Tc_rms_e"), d_wrong_points=diff("wrong_ordering_mean", 100),
             n_worsened_gt10=int(worse.sum()), success_var=fv["crossing_success50"]["value"], success_base=fb["crossing_success50"]["value"],
             d_success=dict(value=float(fv["crossing_success50"]["value"] - fb["crossing_success50"]["value"]), ci95=[float(x) for x in np.nanpercentile(dsucc, [2.5, 97.5])]),
             ordering_success_var=fv["ordering_success_mean"]["value"], failure_var=fv["failure_model_invalid"]["value"])
    c1 = all(r[k]["ci95"][1] < 0 for k in ("d_pair_dFth", "d_pair_dS", "d_near_Tc")); c2 = r["success_var"] >= 0.5
    c3 = r["d_wrong_points"]["value"] <= 0 and r["n_worsened_gt10"] <= 0.02 * base.sum()
    r["criteria"] = dict(c1=c1, c2=c2, c3=c3, verdict="pass" if (c1 and c2 and c3) else ("improves, insufficient for the crossing target" if c1 else "fail"))
    return r


for M in ("F0", "MS"):
    P, valid, blk = {}, {}, {}
    for v in VARS:
        F = np.array([Z[m][f"F_{M}_{v}"] for m in B.ids]); fm = np.array([float(Z[m][f"fmin_{M}_{v}"]) for m in B.ids])
        fm = np.where(np.isfinite(fm), fm, -np.inf); F = np.where(np.isfinite(F), F, 0.0)
        p = tmp / f"{M}_{v}.npz"; np.savez(p, ids=np.array(B.ids), T=T, E0=E0, F=F, fmin=fm, info=json.dumps(dict(model=M, variant=v, setting="S1"))); P[v] = (p, F, fm)
        valid[v] = B.ref_valid & (fm[B.pi] >= GATE) & (fm[B.pj] >= GATE); blk[f"entries_invalid_{v}"] = int((fm < GATE).sum())
    blk["A1_inapplicable_entries"] = int(sum(not bool(Z[m][f"A1_applicable_{M}"]) for m in B.ids))
    npA = np.array([int(Z[m]["npar_A1"]) for m in B.ids]); npR = np.array([int(Z[m][f"npar_{M}_R3"]) for m in B.ids])
    blk["npar"] = dict(A1_median=float(np.median(npA)), A1_p90=float(np.percentile(npA, 90)), A1_max=int(npA.max()), R3_median=float(np.median(npR)), ratio_median=float(np.median(npA / npR)))
    r2 = np.array([float(Z[m][f"r2_{M}_A1"]) for m in B.ids]); blk["A1_r2_median"] = float(np.nanmedian(r2))
    blk["A1_proj_vs_base_proj"] = compare(P, valid, "base_proj", "A1_proj")
    blk["R3_proj_vs_base_proj"] = compare(P, valid, "base_proj", "R3_proj")
    blk["A1_proj_vs_R3_proj"] = compare(P, valid, "R3_proj", "A1_proj")
    blk["oracle6_proj_vs_base_proj"] = compare(P, valid, "base_proj", "oracle6_proj")
    # DEVIATIONS_1c D2: same-cutoff sensitivity, pairs whose two phases both have d_img >= 12.1 A (r_c = 6 A for A1)
    dimg = np.array([float(Z[m]["dimg"]) for m in B.ids]); d12 = (dimg[B.pi] >= 12.1) & (dimg[B.pj] >= 12.1)
    rc = np.array([float(Z[m]["rc"]) for m in B.ids]); blk["rc_distribution"] = dict(p10=float(np.percentile(rc, 10)), median=float(np.median(rc)), n_rc_lt6=int((rc < 6).sum()))
    sub = {v: (valid[v] & d12) for v in VARS}; blk["n_pairs_dimg_ge12"] = int((B.ref_valid & d12).sum())
    blk["d12_A1_proj_vs_R3_proj"] = compare(P, sub, "R3_proj", "A1_proj")
    blk["d12_A1_proj_vs_base_proj"] = compare(P, sub, "base_proj", "A1_proj")
    blk["A1_unproj_vs_base_proj"] = compare(P, valid, "base_proj", "A1")
    res[M] = blk
s_gain = {M: res[M]["A1_proj_vs_base_proj"]["success_var"] - res[M]["R3_proj_vs_base_proj"]["success_var"] for M in ("F0", "MS")}
passed = all(res[M]["A1_proj_vs_base_proj"]["criteria"]["verdict"] == "pass" for M in ("F0", "MS"))
res["decision"] = dict(A1_pass_both=passed, success_gain_over_R3=s_gain,
                       outcome="A1 crosses the representation limit" if passed else
                       ("tensor freedom contributes but is insufficient" if all(g >= 0.05 for g in s_gain.values()) else "tensor freedom is not the main bottleneck"))


def strip(o):
    if isinstance(o, dict): return {k: strip(v) for k, v in o.items()}
    if isinstance(o, (np.floating, np.integer, np.bool_)): return o.item()
    return o


out.write_text(json.dumps(strip(res), indent=1) + "\n")
f = lambda d: f"{d['value']:+.2f} [{d['ci95'][0]:+.2f},{d['ci95'][1]:+.2f}]"
for M in ("F0", "MS"):
    b = res[M]; print(f"== {M}: invalid entries {[(v, b[f'entries_invalid_{v}']) for v in VARS]}; A1 inapplicable {b['A1_inapplicable_entries']}; npar {b['npar']}; A1 r2 med {b['A1_r2_median']:.2f}")
    print(f"  rc {b['rc_distribution']}; pairs with both d_img>=12.1: {b['n_pairs_dimg_ge12']} (success columns below are always full-sample)")
    for k in ("R3_proj_vs_base_proj", "A1_proj_vs_base_proj", "A1_unproj_vs_base_proj", "oracle6_proj_vs_base_proj", "A1_proj_vs_R3_proj",
              "d12_A1_proj_vs_base_proj", "d12_A1_proj_vs_R3_proj"):
        r = b[k]; print(f"  {k:26s} cmp {r['n_compared']} newfail {r['new_failures']} newvalid {r['newly_valid']} | dFth {f(r['d_pair_dFth'])} dS {f(r['d_pair_dS'])} eTc {f(r['d_near_Tc'])} "
                        f"wrong {f(r['d_wrong_points'])} worse>10 {r['n_worsened_gt10']} | succ {r['success_base']:.2f}->{r['success_var']:.2f} ({f(r['d_success'])}) "
                        f"| var {json.dumps({kk: round(v, 3) for kk, v in r['var_metrics'].items()})} -> {r['criteria']['verdict']}")
print("DECISION", json.dumps(res["decision"]))
