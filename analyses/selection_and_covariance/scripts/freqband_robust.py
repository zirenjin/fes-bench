"""Exploratory robustness checks for freqband_analyze.py (post-hoc).

For each variant, pair-level quantities at 300 / 1000 / 1500 K, formula-weighted, 95% cluster bootstrap over formulas:
  share_b       cross-term-inclusive attribution sum(c_b x) / sum(x^2) (as before, now formula-weighted)
  removal_<2THz counterfactual: 1 - <(x - c_low)^2> / <x^2>, the MS reduction if the < 2 THz error alone were removed
  alone_<2THz   <c_low^2> / <x^2>, the MS of the < 2 THz error on its own
Variants:
  sorted        modes matched by sorted index per q-point, assigned to the band of nu_DFT (the original analysis)
  cut0.1/cut0.3 as sorted, but modes with nu_DFT below 0.1 / 0.3 THz dropped on both sides
  no_acoustic   as sorted, but the three lowest branches at every q-point dropped
  matching_free each side's modes binned by their own frequency (no mode matching): c_b = F_b(MLIP) - F_b(DFT)
Usage: freqband_robust.py <FINAL_GROUPS.csv> <freq_dir> <out.json>
"""
import csv, itertools, json, sys
from pathlib import Path
import numpy as np

groups_csv, fdir, out = map(Path, sys.argv[1:4])
H, KB = 4.135667696e-3, 8.617333262e-5
EDGES = [0, 1, 2, 4, 8, 16, np.inf]; NAMES = ["<1", "1-2", "2-4", "4-8", "8-16", ">16"]; TS = [300, 1000, 1500]; CUT = 1e-3
rng = np.random.default_rng(20261102)
groups = {r["formula"]: r["mp_ids"].split(";") for r in csv.DictReader(open(groups_csv))}


def fth(nu, T):
    ok = nu > CUT
    return np.where(ok, KB * T * np.log1p(-np.exp(-np.where(ok, H * nu / (KB * T), 1.0))), 0.0)


def contributions(z, variant, T):
    wd, wm = np.sort(z["freq_dft"], axis=1), np.sort(z["freq_mlip"], axis=1); w = (z["weights"] / z["weights"].sum())[:, None]
    scale = 1000 / int(z["n_prim"]); keep = np.ones_like(wd, bool)
    if variant.startswith("cut"):
        keep = wd >= float(variant[3:])
    elif variant == "no_acoustic":
        keep[:, :3] = False
    if variant == "matching_free":
        bm, bd = np.digitize(wm, EDGES) - 1, np.digitize(wd, EDGES) - 1
        fm, fd = fth(wm, T) * w * scale, fth(wd, T) * w * scale
        return np.array([fm[bm == b].sum() - fd[bd == b].sum() for b in range(len(NAMES))])
    band = np.digitize(wd, EDGES) - 1; dm = (fth(wm, T) - fth(wd, T)) * w * scale
    return np.array([dm[(band == b) & keep].sum() for b in range(len(NAMES))])


Z = {m: np.load(fdir / f"{m}.npz") for ms in groups.values() for m in ms if (fdir / f"{m}.npz").exists()}
pairs = [(f, i, j) for f, ms in groups.items() for i, j in itertools.combinations([m for m in ms if m in Z], 2)]
forms = sorted({f for f, _, _ in pairs}); pidx = {f: [k for k, p in enumerate(pairs) if p[0] == f] for f in forms}
draws = [rng.choice(forms, len(forms)) for _ in range(2000)]
res = dict(n_pairs=len(pairs), n_formulas=len(forms), bands_THz=NAMES)
for variant in ("sorted", "cut0.1", "cut0.3", "no_acoustic", "matching_free"):
    vb = {}
    for t in TS:
        C = {m: contributions(z, variant, t) for m, z in Z.items()}
        A = np.array([C[i] - C[j] for _, i, j in pairs]); x = A.sum(1); low = A[:, :2].sum(1)
        per = {f: dict(cx=(A[pidx[f]] * x[pidx[f], None]).mean(0), x2=(x[pidx[f]] ** 2).mean(),
                       rem=((x[pidx[f]] - low[pidx[f]]) ** 2).mean(), lo2=(low[pidx[f]] ** 2).mean()) for f in forms}

        def stats(fs):
            cx = np.mean([per[f]["cx"] for f in fs], 0); x2 = np.mean([per[f]["x2"] for f in fs])
            return np.r_[cx / x2, 1 - np.mean([per[f]["rem"] for f in fs]) / x2, np.mean([per[f]["lo2"] for f in fs]) / x2]
        point = stats(forms); bs = np.array([stats(s) for s in draws]); lo, hi = np.percentile(bs, [2.5, 97.5], axis=0)
        fmt = lambda k: dict(value=float(point[k]), ci95=[float(lo[k]), float(hi[k])])
        vb[f"T{t}"] = dict(share={n: fmt(b) for b, n in enumerate(NAMES)}, removal_lt2THz=fmt(len(NAMES)), alone_lt2THz=fmt(len(NAMES) + 1),
                           pair_rms_retained_meV=float(np.sqrt(np.mean(x ** 2))))
    res[variant] = vb
out.write_text(json.dumps(res, indent=1) + "\n")
for v in ("sorted", "cut0.1", "cut0.3", "no_acoustic", "matching_free"):
    b = res[v]["T1000"]
    print(v, "rms", round(b["pair_rms_retained_meV"], 2), "share<1,1-2", [round(b["share"][n]["value"], 2) for n in NAMES[:2]],
          "removal<2", round(b["removal_lt2THz"]["value"], 2), [round(c, 2) for c in b["removal_lt2THz"]["ci95"]],
          "alone<2", round(b["alone_lt2THz"]["value"], 2))
