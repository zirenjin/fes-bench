"""Exploratory (post-hoc): which frequency bands carry the F0 thermal free-energy error, single-phase and pair-level.

Modes are matched by sorted index at each q-point (phonopy order); each mode's thermal free-energy error
f_th(nu_MLIP) - f_th(nu_DFT), with f_th(nu, T) = kT ln(1 - exp(-h nu / kT)), is assigned to the band of nu_DFT.
Variance shares: share_b = sum(c_b * x) / sum(x^2) with x = sum_b c_b, so the shares add to 1.
Usage: freqband_analyze.py <FINAL_GROUPS.csv> <freq_dir> <fixed_dir> <out.json>
"""
import csv, itertools, json, sys
from pathlib import Path
import numpy as np

groups_csv, fdir, fixed, out = map(Path, sys.argv[1:5])
H, KB = 4.135667696e-3, 8.617333262e-5          # eV/THz, eV/K
EDGES = [0, 1, 2, 4, 8, 16, np.inf]; NAMES = ["<1", "1-2", "2-4", "4-8", "8-16", ">16"]; TS = [300, 1000, 1500]; CUT = 1e-3
groups = {r["formula"]: r["mp_ids"].split(";") for r in csv.DictReader(open(groups_csv))}


def fth(nu, T):
    x = np.where(nu > CUT, H * nu / (KB * T), np.inf)
    return np.where(nu > CUT, KB * T * np.log1p(-np.exp(-x)), 0.0)            # eV per mode


C, rel, closure, nmode = {}, {b: [] for b in NAMES}, [], {}
for ms in groups.values():
    for m in ms:
        p = fdir / f"{m}.npz"
        if not p.exists():
            continue
        z = np.load(p); wd, wm, w = np.sort(z["freq_dft"], axis=1), np.sort(z["freq_mlip"], axis=1), z["weights"] / z["weights"].sum()
        band = np.digitize(wd, EDGES) - 1; scale = 1000 / int(z["n_prim"])      # eV per primitive cell -> meV/atom
        C[m] = {}
        for t in TS:
            dm = (fth(wm, t) - fth(wd, t)) * w[:, None] * scale
            C[m][t] = np.array([dm[band == b].sum() for b in range(len(NAMES))])
        nmode[m] = np.array([(w[:, None] * (band == b)).sum() for b in range(len(NAMES))]) / wd.shape[1]
        ok = wd > 0.05
        for b, name in enumerate(NAMES):
            sel = ok & (band == b)
            if sel.any():
                rel[name].append(((wm[sel] - wd[sel]) / wd[sel], np.broadcast_to(w[:, None], wd.shape)[sel]))
        f = np.load(fixed / f"{m}.npz"); T = f["T"]; k = int(np.flatnonzero(T == 1000)[0])
        ref = (f["F_mlip"][k] - f["F_mlip"][0]) - (f["F_dft"][k] - f["F_dft"][0])
        closure.append(abs(C[m][1000].sum() - ref))

res = dict(n_entries=len(C), bands_THz=NAMES, closure_vs_fixed_run_1000K_meV=dict(median=float(np.median(closure)), p95=float(np.percentile(closure, 95))),
           mode_fraction_by_band=[float(v) for v in np.mean(list(nmode.values()), axis=0)])
res["relative_freq_error_by_band"] = {b: dict(weighted_mean=float(np.average(np.concatenate([r for r, _ in v]), weights=np.concatenate([w for _, w in v]))),
                                              weighted_mean_abs=float(np.average(np.abs(np.concatenate([r for r, _ in v])), weights=np.concatenate([w for _, w in v]))))
                                      for b, v in rel.items() if v}
pairs = [(i, j) for ms in groups.values() for i, j in itertools.combinations([m for m in ms if m in C], 2)]
for t in TS:
    S = np.array([C[m][t] for m in C]); Pm = np.array([C[i][t] - C[j][t] for i, j in pairs])
    share = lambda A: [float(v) for v in (A * A.sum(1, keepdims=True)).sum(0) / (A.sum(1) ** 2).sum()]
    res[f"T{t}"] = dict(single_rms=float(np.sqrt((S.sum(1) ** 2).mean())), pair_rms=float(np.sqrt((Pm.sum(1) ** 2).mean())),
                        single_share=share(S), pair_share=share(Pm),
                        band_cancellation_c=[float(1 - np.mean(Pm[:, b] ** 2) / np.mean([C[i][t][b] ** 2 + C[j][t][b] ** 2 for i, j in pairs]))
                                             if np.mean([C[i][t][b] ** 2 + C[j][t][b] ** 2 for i, j in pairs]) > 0 else None for b in range(len(NAMES))])
res["n_pairs"] = len(pairs)
out.write_text(json.dumps(res, indent=1) + "\n"); print(json.dumps(res, indent=1))
