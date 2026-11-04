import csv, itertools, json, sys
from pathlib import Path
import numpy as np
groups = {r["formula"]: r["mp_ids"].split(";") for r in csv.DictReader(open("metadata/selection/FINAL_GROUPS.csv"))}
TS = [300, 600, 1000, 1500]
for head in ["MP_traj_v024_alldata_mixu", "Omat24"]:
    E = {}
    for f, ms in groups.items():
        for m in ms:
            p = Path(f"mlip/{head}/{m}.npz")
            if not p.exists(): continue
            z = np.load(p); T = z["T"]; Fd, Fm = z["F_dft"], z["F_mlip"]
            E[m] = {t: (float(Fd[T == t][0] - Fd[0]), float(Fm[T == t][0] - Fm[0])) for t in TS}
    print(head)
    for t in TS:
        sd = np.array([E[m][t][0] for m in E]); sm = np.array([E[m][t][1] for m in E])
        D, M = [], []
        for f, ms in groups.items():
            ok = [m for m in ms if m in E]
            for i, j in itertools.combinations(ok, 2):
                D.append(E[i][t][0] - E[j][t][0]); M.append(E[i][t][1] - E[j][t][1])
        D, M = np.array(D), np.array(M)
        o = lambda y, x: float(np.sum(x * y) / np.sum(x * x))   # through-origin, orientation invariant
        r = float(np.sum(D * M) / np.sqrt(np.sum(D * D) * np.sum(M * M)))
        print(f" {t}K single-phase slope MLIP~DFT {o(sm, sd):.3f} | pair: MLIP~DFT {o(M, D):.3f}  DFT~MLIP {o(D, M):.3f}  r {r:.3f}  rms D {np.sqrt(np.mean(D**2)):.2f} M {np.sqrt(np.mean(M**2)):.2f} err {np.sqrt(np.mean((M-D)**2)):.2f}")
