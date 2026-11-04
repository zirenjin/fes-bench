"""PROTOCOL.md (gate 1) scoring with the FES-Bench v0.1 scorer, setting S1 (E0 from the reference).

For each model M: base pairs = reference-valid and valid under none_R0. Paired comparisons (variant - none_R0) use base
pairs on which the variant is also valid; variant-invalid base pairs are counted as new failures. Crossing success is
over all reference-valid crossing pairs (variant-valid AND a crossing within +-50 K).
Usage: gate1_score.py <gate1_data_dir> <out.json>
"""
import json, sys, tempfile
from pathlib import Path
import numpy as np
V01 = Path("C:/users/ziren/documents/dpa-adapt/outputs/fes-bench-v0.1-20261102"); sys.path.insert(0, str(V01))
from score import Bench, GATE  # noqa: E402

data, out = Path(sys.argv[1]), Path(sys.argv[2]); B = Bench(); T = B.T; it4 = [B.it[t] for t in (300, 600, 1000, 1500)]
ref = np.load(V01 / "predictions/reference_PBE.npz"); E0 = ref["E0"]
Z = {m: np.load(data / f"{m}.npz") for m in B.ids}
assert max(np.abs(Z[m]["F_DFT"] - ref["F"][n]).max() for n, m in enumerate(B.ids)) < 1e-6, "oracle path does not reproduce the reference"
VAR = ("none_raw", "none_R0", "R1", "R2", "R3"); tmp = Path(tempfile.mkdtemp()); res = dict(protocol_sha256=open(Path(__file__).parents[1] / "PROTOCOL.sha256").read().split()[0])


def save(M, v):
    F = np.array([Z[m][f"F_{M}_{v}"] for m in B.ids]); fmin = np.array([float(Z[m][f"fmin_{M}_{v}"]) for m in B.ids])
    fmin = np.where(np.isfinite(fmin), fmin, -np.inf); F = np.where(np.isfinite(F), F, 0.0)
    p = tmp / f"{M}_{v}.npz"; np.savez(p, ids=np.array(B.ids), T=T, E0=E0, F=F, fmin=fmin, info=json.dumps(dict(model=M, variant=v, setting="S1"))); return p, F, fmin


for M in ("F0", "MS"):
    P = {v: save(M, v) for v in VAR}
    valid = {v: B.ref_valid & (P[v][2][B.pi] >= GATE) & (P[v][2][B.pj] >= GATE) for v in VAR}
    base = valid["none_R0"]; blk = dict(n_reference_valid=int(B.ref_valid.sum()), n_base_pairs=int(base.sum()))
    G = {v: E0[:, None] + P[v][1] for v in VAR}
    for v in VAR:
        mask = base & valid[v]; s = B.score(P[v][0], mask); b = B.score(P["none_R0"][0], mask); m = s["metrics"]; mb = b["metrics"]
        def diff(k, sc=1.0):
            dr = (m[k]["_draws"] - mb[k]["_draws"]) * sc; lo, hi = np.nanpercentile(dr, [2.5, 97.5]); return dict(value=float((m[k]["value"] - mb[k]["value"]) * sc), ci95=[float(lo), float(hi)])
        e_v = (G[v][B.pi] - G[v][B.pj]) - B.dref; e_b = (G["none_R0"][B.pi] - G["none_R0"][B.pj]) - B.dref
        worse = mask & ((np.abs(e_v[:, it4]) - np.abs(e_b[:, it4])).max(1) > 10)
        full = B.score(P[v][0], "model_valid")["metrics"]   # success metrics are over all reference-valid pairs regardless of mask
        r = dict(n_compared_pairs=int(mask.sum()), new_failures_among_base=int((base & ~valid[v]).sum()), newly_valid=int((B.ref_valid & ~base & valid[v]).sum()),
                 pair_dFth_rms=m["pair_dFth_rms_1000K"]["value"], pair_dS_rms=m["pair_dS_rms_1000K"]["value"], near_Tc_rms=m["near_Tc_rms_e"]["value"],
                 wrong_ordering_mean=m["wrong_ordering_mean"]["value"], within_eps50=m["near_Tc_within_energy_budget50"]["value"], hit50=m["crossing_hit50"]["value"],
                 crossing_success50=full["crossing_success50"]["value"], crossing_success50_ci=full["crossing_success50"]["ci95"],
                 ordering_success=full["ordering_success_mean"]["value"], failure=full["failure_model_invalid"]["value"],
                 d_pair_dFth=diff("pair_dFth_rms_1000K"), d_pair_dS=diff("pair_dS_rms_1000K"), d_near_Tc=diff("near_Tc_rms_e"),
                 d_wrong_points=diff("wrong_ordering_mean", 100), n_worsened_gt10=int(worse.sum()), events=m["events"])
        if v in ("R1", "R2", "R3"):
            c1 = all(r[k]["ci95"][1] < 0 for k in ("d_pair_dFth", "d_pair_dS", "d_near_Tc")); c2 = r["crossing_success50"] >= 0.5
            c3 = r["d_wrong_points"]["value"] <= 0 and r["n_worsened_gt10"] <= 0.02 * base.sum()
            r["criteria"] = dict(c1_three_errors_improve=c1, c2_crossing_success_ge_0p5=c2, c3_no_new_harm=c3,
                                 verdict="pass" if (c1 and c2 and c3) else ("improves, insufficient for the transition target" if c1 else "fail"))
            r2 = [float(Z[mm][f"r2_{M}_{v}"]) for mm in B.ids]; r2f = [float(Z[mm][f"r2_fullfc_{M}"]) for mm in B.ids]; npar = [int(Z[mm][f"npar_{M}_{v}"]) for mm in B.ids]
            r["fit"] = dict(median_r2=float(np.nanmedian(r2)), median_r2_over_r2_fullfc=float(np.nanmedian(np.array(r2) / np.maximum(r2f, 1e-9))),
                            median_npar=float(np.median(npar)), max_npar=int(max(npar)), n_inapplicable=int(np.sum(~np.isfinite(r2))))
        blk[v] = r
    if True:   # descriptive R0 control: none_R0 vs none_raw on pairs valid for both
        mask = base & valid["none_raw"]; a = B.score(P["none_R0"][0], mask)["metrics"]; b = B.score(P["none_raw"][0], mask)["metrics"]
        blk["R0_control_none_R0_minus_raw"] = {k: a[k]["value"] - b[k]["value"] for k in ("pair_dFth_rms_1000K", "pair_dS_rms_1000K", "near_Tc_rms_e", "wrong_ordering_mean")}
    res[M] = blk


def strip(o):
    if isinstance(o, dict): return {k: strip(v) for k, v in o.items()}
    if isinstance(o, (np.floating, np.integer, np.bool_)): return o.item()
    return o


out.write_text(json.dumps(strip(res), indent=1) + "\n")
for M in ("F0", "MS"):
    print(f"== {M}: base pairs {res[M]['n_base_pairs']} of {res[M]['n_reference_valid']}; R0 control {json.dumps({k: round(v, 3) for k, v in res[M]['R0_control_none_R0_minus_raw'].items()})}")
    for v in VAR:
        r = res[M][v]; f = lambda d: f"{d['value']:+.2f} [{d['ci95'][0]:+.2f},{d['ci95'][1]:+.2f}]"
        print(f"  {v:8s} n={r['n_compared_pairs']} newfail={r['new_failures_among_base']} newvalid={r['newly_valid']} pairFth {r['pair_dFth_rms']:.2f} ({f(r['d_pair_dFth'])}) "
              f"dS {r['pair_dS_rms']:.1f} ({f(r['d_pair_dS'])}) eTc {r['near_Tc_rms']:.2f} ({f(r['d_near_Tc'])}) wrong {100 * r['wrong_ordering_mean']:.1f} ({f(r['d_wrong_points'])}) "
              f"eps50 {r['within_eps50']:.2f} hit50 {r['hit50']:.2f} succ50 {r['crossing_success50']:.2f} worse>10 {r['n_worsened_gt10']}"
              + (f" -> {r['criteria']['verdict']} | fit r2 {r['fit']['median_r2']:.2f}, r2/r2full {r['fit']['median_r2_over_r2_fullfc']:.2f}, npar med {r['fit']['median_npar']:.0f}" if "criteria" in r else ""))
