"""Paired F0 vs MatterSim contrast (post hoc, descriptive; requested after RESULT.md, not in PROTOCOL.md).

Both models are scored on the same formula draws of one cluster bootstrap (2000 draws, formulas resampled with
replacement), so ratios and differences carry the pairing. Metric definitions follow analyze_covariance.py:
- single-phase / pair dF_th RMS: pooled over entries / pairs (as in Q1 rms_single_vs_pair_meV);
- cancellation c: formula-weighted mean of the per-formula 1 - <dD^2>/<d_i^2 + d_j^2>;
- wrong ordering (none): formula-weighted, non-degenerate pairs (|reference gap| >= 1 meV/atom).
Also reported: formula-weighted RMS (sqrt of the mean over formulas of the per-formula mean square; the weighting
PROTOCOL.md names for the primary metric) and the fraction of formulas whose pair mean square is larger for MatterSim.
Subsets: 'all' (the 521 primary entries) and 'common_stable' (entries whose dense-mesh minimum frequency at the
DFT geometry is >= GATE for DFT, F0 and MatterSim alike; GATE = -0.05 THz, the threshold already used in RESULT.md,
fixed before this script was run). Pairs need both entries to pass; formulas need at least one surviving pair.
Usage: paired_model_contrast.py <F0_cov_dir> <MS_cov_dir> <FMIN_ENTRIES.csv> <out.json> [--gate -0.05]
"""
import csv, json, sys
from collections import defaultdict
from pathlib import Path
import numpy as np

f0_dir, ms_dir, fmin_csv, out = map(Path, sys.argv[1:5])
GATE = float(sys.argv[sys.argv.index("--gate") + 1]) if "--gate" in sys.argv else -0.05
TS = [300, 600, 1000, 1500]; NB = 2000


def load(d):
    ent = {r["mp_id"]: {"formula": r["formula"], **{t: float(r[f"dFth_{t}K"]) for t in TS}} for r in csv.DictReader(open(d / "ENTRIES.csv"))}
    prs = {(r["i"], r["j"]): dict(formula=r["formula"], a=float(r["a"]), **{f"b{t}": float(r[f"b{t}"]) for t in TS},
                                  **{f"ref{t}": float(r[f"ref{t}"]) for t in TS}) for r in csv.DictReader(open(d / "PAIRS.csv"))}
    return ent, prs


E = {"F0": None, "MS": None}; P = {}
E["F0"], P["F0"] = load(f0_dir); E["MS"], P["MS"] = load(ms_dir)
assert set(E["F0"]) == set(E["MS"]) and set(P["F0"]) == set(P["MS"])
fmin = {r["mp_id"]: r for r in csv.DictReader(open(fmin_csv))}
assert all(abs(float(r["fmin_dft"]) - float(r["fmin_dft_ms"])) < 1e-6 for r in fmin.values())
passes = {m: min(float(r["fmin_dft"]), float(r["fmin_f0"]), float(r["fmin_ms"])) >= GATE for m, r in fmin.items()}
fails = {k: sum(float(r[c]) < GATE for r in fmin.values()) for k, c in (("dft", "fmin_dft"), ("f0", "fmin_f0"), ("ms", "fmin_ms"))}


def stats(pair_keys):
    """Per-formula sufficient statistics for both models."""
    byf = defaultdict(list)
    for k in pair_keys: byf[P["F0"][k]["formula"]].append(k)
    fs = sorted(byf); S = {}
    for mdl in ("F0", "MS"):
        e, p = E[mdl], P[mdl]; S[mdl] = {}
        for f in fs:
            ks = byf[f]; ents = sorted({m for k in ks for m in k}); d = {}
            for t in TS:
                di = np.array([e[i][t] for i, _ in ks]); dj = np.array([e[j][t] for _, j in ks]); b = np.array([p[k][f"b{t}"] for k in ks])
                den = np.mean(di ** 2 + dj ** 2)
                nd = [k for k in ks if abs(p[k][f"ref{t}"]) >= 1.0]
                wr = np.mean([np.sign(p[k][f"ref{t}"] + p[k]["a"] + p[k][f"b{t}"]) != np.sign(p[k][f"ref{t}"]) for k in nd]) if nd else np.nan
                d[t] = dict(ss_e=float(sum(e[m][t] ** 2 for m in ents)), n_e=len(ents), ss_p=float(np.sum(b ** 2)), n_p=len(ks),
                            c=float(1 - np.mean(b ** 2) / den) if den > 0 else np.nan, wrong=float(wr))
            S[mdl][f] = d
    return fs, S


def metrics(S, fs, t):
    """fs may contain repeats (bootstrap draw)."""
    ss_e = sum(S[f][t]["ss_e"] for f in fs); n_e = sum(S[f][t]["n_e"] for f in fs)
    ss_p = sum(S[f][t]["ss_p"] for f in fs); n_p = sum(S[f][t]["n_p"] for f in fs)
    return dict(single_rms=np.sqrt(ss_e / n_e), pair_rms=np.sqrt(ss_p / n_p), c=np.nanmean([S[f][t]["c"] for f in fs]),
                wrong=np.nanmean([S[f][t]["wrong"] for f in fs]),
                single_rms_fw=np.sqrt(np.mean([S[f][t]["ss_e"] / S[f][t]["n_e"] for f in fs])),
                pair_rms_fw=np.sqrt(np.mean([S[f][t]["ss_p"] / S[f][t]["n_p"] for f in fs])))


def contrast(pair_keys, seed):
    fs, S = stats(pair_keys); rng = np.random.default_rng(seed)
    draws = [list(rng.choice(fs, len(fs))) for _ in range(NB)]
    res = dict(n_formulas=len(fs), n_pairs=len(pair_keys), n_entries=len({m for k in pair_keys for m in k}))
    for t in TS:
        pt = {mdl: metrics(S[mdl], fs, t) for mdl in ("F0", "MS")}
        bd = [{mdl: metrics(S[mdl], dr, t) for mdl in ("F0", "MS")} for dr in draws]
        def summ(fn):
            v = np.array([fn(x) for x in bd]); lo, hi = np.nanpercentile(v, [2.5, 97.5])
            return dict(value=float(fn(pt)), ci95=[float(lo), float(hi)])
        r_single = lambda x: x["MS"]["single_rms"] / x["F0"]["single_rms"]
        r_pair = lambda x: x["MS"]["pair_rms"] / x["F0"]["pair_rms"]
        res[f"T{t}"] = dict(
            F0={k: float(v) for k, v in pt["F0"].items()}, MS={k: float(v) for k, v in pt["MS"].items()},
            single_rms_ratio_MS_over_F0=summ(r_single), pair_rms_ratio_MS_over_F0=summ(r_pair),
            c_diff_MS_minus_F0=summ(lambda x: x["MS"]["c"] - x["F0"]["c"]),
            wrong_diff_MS_minus_F0_points=summ(lambda x: 100 * (x["MS"]["wrong"] - x["F0"]["wrong"])),
            frac_draws_inversion=float(np.mean([r_single(x) < 1 and r_pair(x) > 1 for x in bd])),
            single_rms_fw_ratio_MS_over_F0=summ(lambda x: x["MS"]["single_rms_fw"] / x["F0"]["single_rms_fw"]),
            pair_rms_fw_ratio_MS_over_F0=summ(lambda x: x["MS"]["pair_rms_fw"] / x["F0"]["pair_rms_fw"]),
            frac_draws_inversion_fw=float(np.mean([x["MS"]["single_rms_fw"] < x["F0"]["single_rms_fw"] and x["MS"]["pair_rms_fw"] > x["F0"]["pair_rms_fw"] for x in bd])),
            frac_formulas_MS_pair_worse=float(np.mean([S["MS"][f][t]["ss_p"] > S["F0"][f][t]["ss_p"] for f in fs])),
            frac_formulas_MS_single_better=float(np.mean([S["MS"][f][t]["ss_e"] < S["F0"][f][t]["ss_e"] for f in fs])))
    return res


all_keys = sorted(P["F0"]); stable_keys = [k for k in all_keys if passes[k[0]] and passes[k[1]]]
all_formulas = {P["F0"][k]["formula"] for k in all_keys}; st_formulas = {P["F0"][k]["formula"] for k in stable_keys}
res = dict(status="post hoc, descriptive; not in PROTOCOL.md", gate_THz=GATE, n_bootstrap=NB,
           entries_failing_gate=fails, entries_failing_any=sum(not v for v in passes.values()),
           all=contrast(all_keys, 20261102), common_stable=contrast(stable_keys, 20261102))
res["common_stable_coverage"] = dict(entries=res["common_stable"]["n_entries"] / res["all"]["n_entries"],
                                     pairs=len(stable_keys) / len(all_keys), formulas=len(st_formulas) / len(all_formulas))
out.write_text(json.dumps(res, indent=1) + "\n")
for sub in ("all", "common_stable"):
    r = res[sub]; print(f"== {sub}: {r['n_formulas']} formulas, {r['n_pairs']} pairs, {r['n_entries']} entries")
    for t in TS:
        x = r[f"T{t}"]; f = lambda d: f"{d['value']:.3f} [{d['ci95'][0]:.3f}, {d['ci95'][1]:.3f}]"
        print(f"  {t} K  single F0/MS {x['F0']['single_rms']:.2f}/{x['MS']['single_rms']:.2f} ratio {f(x['single_rms_ratio_MS_over_F0'])} | "
              f"pair {x['F0']['pair_rms']:.2f}/{x['MS']['pair_rms']:.2f} ratio {f(x['pair_rms_ratio_MS_over_F0'])} | "
              f"c {x['F0']['c']:.3f}/{x['MS']['c']:.3f} diff {f(x['c_diff_MS_minus_F0'])} | "
              f"wrong {100*x['F0']['wrong']:.1f}/{100*x['MS']['wrong']:.1f} diff {f(x['wrong_diff_MS_minus_F0_points'])} | inv {x['frac_draws_inversion']:.3f}")
        print(f"         fw single {x['F0']['single_rms_fw']:.2f}/{x['MS']['single_rms_fw']:.2f} ratio {f(x['single_rms_fw_ratio_MS_over_F0'])} | "
              f"fw pair {x['F0']['pair_rms_fw']:.2f}/{x['MS']['pair_rms_fw']:.2f} ratio {f(x['pair_rms_fw_ratio_MS_over_F0'])} | inv_fw {x['frac_draws_inversion_fw']:.3f} | "
              f"formulas: MS pair worse {x['frac_formulas_MS_pair_worse']:.2f}, MS single better {x['frac_formulas_MS_single_better']:.2f}")
print("fails", fails, "any", res["entries_failing_any"], "coverage", res["common_stable_coverage"])
