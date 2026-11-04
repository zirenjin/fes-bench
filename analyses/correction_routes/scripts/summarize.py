"""Tables for RESULT.md from results/ROUTES.json (repeat 1 primary; range over repeats 1-3 for stability).
Usage: summarize.py <ROUTES.json>
"""
import json, sys
import numpy as np

R = json.load(open(sys.argv[1])); runs = R["runs"]
print("gate", json.dumps(R["gate"]), "C runnable", R["C_runnable"], R["n"])
for S in ("S0", "S1"):
    print(f"\n## {S}")
    print("| model | route | M1 RMS e(Tc) | Δ vs none [95%] | within ε50 | wrong mean | Δ points [95%] | pair RMS 300/1000/1500 | δΔS 1000 | p95 e1000 | n>20 | worsened>10 | hits/missed/false | cov 1000 | verdict | M1Δ range reps | M2Δ range reps |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for M in ("F0", "MS", "D"):
        routes = ("none", "A", "B", "C") if M != "D" else ("D1", "D2", "D3")
        for r in routes:
            key = f"rep1/{M}/{S}/{r}"
            if key not in runs: continue
            blk = runs[key]
            ms = [("", blk)] if r == "none" else [(b, m) for b, m in blk.items() if b.startswith("vs_")]
            for b, m in ms:
                ev = m["M5_events"]; v = m.get("vs_none")
                def rng(f):
                    vals = [f(runs[f"rep{k}/{M}/{S}/{r}"][b]) for k in (1, 2, 3)]; return f"{min(vals):+.2f}..{max(vals):+.2f}"
                print(f"| {M}{(' ' + b) if b else ''} | {r} | {m['M1_rms_e_tc']:.2f} | "
                      + ("—" if v is None else f"{v['M1_rms_diff']['value']:+.2f} [{v['M1_rms_diff']['ci95'][0]:+.2f}, {v['M1_rms_diff']['ci95'][1]:+.2f}]")
                      + f" | {m['M1_within_eps50']:.2f} | {100 * m['M2_wrong_mean']:.1f} | "
                      + ("—" if v is None else f"{v['M2_wrong_mean_diff_points']['value']:+.1f} [{v['M2_wrong_mean_diff_points']['ci95'][0]:+.1f}, {v['M2_wrong_mean_diff_points']['ci95'][1]:+.1f}]")
                      + f" | {m['M3_pair_rms']['300']:.2f}/{m['M3_pair_rms']['1000']:.2f}/{m['M3_pair_rms']['1500']:.2f} | {m['M4_dS_rms']['1000']:.1f} | {m['M6_p95_e1000']:.1f} | {m['M6_n_max_gt20']} | "
                      + f"{m.get('M6_n_worsened_gt10', '—')} | {ev['hits']}/{ev['missed']}/{ev['false']} | {m['M7_conformal']['1000']['coverage']:.2f} | "
                      + f"{m['verdict']['label'] if 'verdict' in m else '—'} | "
                      + ("" if v is None else rng(lambda x: x['vs_none']['M1_rms_diff']['value'])) + " | "
                      + ("" if v is None else rng(lambda x: x['vs_none']['M2_wrong_mean_diff_points']['value'])) + " |")
print("\n## M8 fit vs transfer (repeat 1)")
for key, blk in runs.items():
    if key.startswith("rep1/") and "M8" in blk:
        x = blk["M8"]
        print(key, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in x.items() if not k.endswith("cfgs")},
              "cfg2", sorted(set(map(tuple, x["stage2_cfgs"]))), "cfg1", sorted(set(map(tuple, x.get("stage1_cfgs", [])))))
    if key.startswith("rep1/D/"):
        print(key, "params", blk.get("D_params"))
