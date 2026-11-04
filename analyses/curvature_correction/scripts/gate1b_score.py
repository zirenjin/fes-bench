"""PROTOCOL_1b.md scoring: projected R3 vs projected baseline (and the unprojected pair for continuity), with the gate-1
criteria unchanged, plus the step-1 self-check distributions and minimum-frequency locations.
Usage: gate1b_score.py <data1b_dir> <out.json>
"""
import json, sys, tempfile
from pathlib import Path
import numpy as np
V01 = Path("C:/users/ziren/documents/dpa-adapt/outputs/fes-bench-v0.1-20261102"); sys.path.insert(0, str(V01))
from score import Bench, GATE  # noqa: E402

data, out = Path(sys.argv[1]), Path(sys.argv[2]); B = Bench(); T = B.T; it4 = [B.it[t] for t in (300, 600, 1000, 1500)]
ref = np.load(V01 / "predictions/reference_PBE.npz"); E0 = ref["E0"]; Z = {m: np.load(data / f"{m}.npz") for m in B.ids}
assert max(np.abs(Z[m]["F_DFT"] - ref["F"][n]).max() for n, m in enumerate(B.ids)) < 1e-6
tmp = Path(tempfile.mkdtemp()); res = dict(protocol_sha256=open(Path(__file__).parents[1] / "PROTOCOL_1b.sha256").read().split()[0])
q = lambda a: dict(median=float(np.median(a)), p90=float(np.percentile(a, 90)), max=float(np.max(a)))
for M in ("F0", "MS"):
    blk = {}
    blk["self_check"] = {k: q([float(Z[m][f"{k}_{M}"]) if k == "assembly_reldiff" else float(Z[m][f"{k}_{M}_{t}"]) for m in B.ids])
                         for k, t in (("assembly_reldiff", None),)}
    for k in ("asr", "perm", "sg"):
        for t in ("base", "R3"): blk["self_check"][f"{k}_{t}"] = q([float(Z[m][f"{k}_{M}_{t}"]) for m in B.ids])
    for t in ("dK",): blk["self_check"]["asr_dK"] = q([float(Z[m][f"asr_{M}_dK"]) for m in B.ids]); blk["self_check"]["perm_dK"] = q([float(Z[m][f"perm_{M}_dK"]) for m in B.ids])
    P, valid = {}, {}
    for v in ("base", "R3", "base_proj", "R3_proj"):
        F = np.array([Z[m][f"F_{M}_{v}"] for m in B.ids]); fm = np.array([float(Z[m][f"fmin_{M}_{v}"]) for m in B.ids]); qm = np.array([float(Z[m][f"qmin_{M}_{v}"]) for m in B.ids])
        p = tmp / f"{M}_{v}.npz"; np.savez(p, ids=np.array(B.ids), T=T, E0=E0, F=F, fmin=fm, info=json.dumps(dict(model=M, variant=v, setting="S1"))); P[v] = (p, F, fm, qm)
        valid[v] = B.ref_valid & (fm[B.pi] >= GATE) & (fm[B.pj] >= GATE)
        inv = fm < GATE; blk[f"entries_invalid_{v}"] = int(inv.sum()); blk[f"invalid_near_gamma_{v}"] = int((inv & (qm < 0.05)).sum())
        blk[f"invalid_fmin_quantiles_{v}"] = [float(x) for x in np.percentile(fm[inv], [10, 50, 90])] if inv.any() else None
    for base_v, var_v in (("base", "R3"), ("base_proj", "R3_proj")):
        base = valid[base_v]; mask = base & valid[var_v]; s = B.score(P[var_v][0], mask)["metrics"]; b = B.score(P[base_v][0], mask)["metrics"]
        def diff(k, sc=1.0):
            dr = (s[k]["_draws"] - b[k]["_draws"]) * sc; lo, hi = np.nanpercentile(dr, [2.5, 97.5]); return dict(value=float((s[k]["value"] - b[k]["value"]) * sc), ci95=[float(lo), float(hi)])
        Gv = E0[:, None] + P[var_v][1]; Gb = E0[:, None] + P[base_v][1]
        e_v = (Gv[B.pi] - Gv[B.pj]) - B.dref; e_b = (Gb[B.pi] - Gb[B.pj]) - B.dref
        worse = mask & ((np.abs(e_v[:, it4]) - np.abs(e_b[:, it4])).max(1) > 10)
        fv = B.score(P[var_v][0], "model_valid")["metrics"]; fb = B.score(P[base_v][0], "model_valid")["metrics"]
        r = dict(n_base_pairs=int(base.sum()), n_compared=int(mask.sum()), new_failures=int((base & ~valid[var_v]).sum()), newly_valid=int((B.ref_valid & ~base & valid[var_v]).sum()),
                 base_metrics={k: b[k]["value"] for k in ("pair_dFth_rms_1000K", "pair_dS_rms_1000K", "near_Tc_rms_e", "wrong_ordering_mean", "crossing_hit50")},
                 var_metrics={k: s[k]["value"] for k in ("pair_dFth_rms_1000K", "pair_dS_rms_1000K", "near_Tc_rms_e", "wrong_ordering_mean", "crossing_hit50")},
                 d_pair_dFth=diff("pair_dFth_rms_1000K"), d_pair_dS=diff("pair_dS_rms_1000K"), d_near_Tc=diff("near_Tc_rms_e"), d_wrong_points=diff("wrong_ordering_mean", 100),
                 n_worsened_gt10=int(worse.sum()), crossing_success50_variant=fv["crossing_success50"]["value"], crossing_success50_ci=fv["crossing_success50"]["ci95"],
                 crossing_success50_base=fb["crossing_success50"]["value"], failure_variant=fv["failure_model_invalid"]["value"], failure_base=fb["failure_model_invalid"]["value"])
        c1 = all(r[k]["ci95"][1] < 0 for k in ("d_pair_dFth", "d_pair_dS", "d_near_Tc")); c2 = r["crossing_success50_variant"] >= 0.5
        c3 = r["d_wrong_points"]["value"] <= 0 and r["n_worsened_gt10"] <= 0.02 * base.sum()
        r["criteria"] = dict(c1=c1, c2=c2, c3=c3, verdict="pass" if (c1 and c2 and c3) else ("improves, insufficient for the crossing target" if c1 else "fail"))
        blk[f"{var_v}_vs_{base_v}"] = r
    res[M] = blk


def strip(o):
    if isinstance(o, dict): return {k: strip(v) for k, v in o.items()}
    if isinstance(o, (np.floating, np.integer, np.bool_)): return o.item()
    return o


out.write_text(json.dumps(strip(res), indent=1) + "\n")
for M in ("F0", "MS"):
    b = res[M]; print(f"== {M} self-check:", json.dumps({k: round(v["median"], 8) for k, v in b["self_check"].items()}))
    print("   invalid entries:", {v: (b[f"entries_invalid_{v}"], b[f"invalid_near_gamma_{v}"], b[f"invalid_fmin_quantiles_{v}"]) for v in ("base", "R3", "base_proj", "R3_proj")})
    for k in ("R3_vs_base", "R3_proj_vs_base_proj"):
        r = b[k]; f = lambda d: f"{d['value']:+.2f} [{d['ci95'][0]:+.2f},{d['ci95'][1]:+.2f}]"
        print(f"   {k}: base {r['n_base_pairs']} compared {r['n_compared']} newfail {r['new_failures']} newvalid {r['newly_valid']} | dFth {f(r['d_pair_dFth'])} dS {f(r['d_pair_dS'])} "
              f"eTc {f(r['d_near_Tc'])} wrong {f(r['d_wrong_points'])} worse>10 {r['n_worsened_gt10']} | success base {r['crossing_success50_base']:.2f} var {r['crossing_success50_variant']:.2f} "
              f"| var metrics {json.dumps({k2: round(v, 3) for k2, v in r['var_metrics'].items()})} -> {r['criteria']['verdict']}")
