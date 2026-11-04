"""Supplementary (not in PROTOCOL.md): all Alexandria arms restricted to the Tc-applicable pairs, same temperature exclusion
(|T - Tc| >= 100 K), so the Tc arm is compared on the same pairs and points.
Usage: supp_same_subset.py <FINAL_GROUPS.csv> <F0_mlip_dir> <out.json>
"""
import csv, itertools, json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).parent))
from anchor_common import tc_anchor  # noqa: E402

groups = {r["formula"]: r["mp_ids"].split(";") for r in csv.DictReader(open(sys.argv[1]))}; mdir = Path(sys.argv[2])
E = {m: np.load(mdir / f"{m}.npz") for ms in groups.values() for m in ms}
T = E[next(iter(E))]["T"]; TE = [600, 1000, 1500]; it = {t: int(np.flatnonzero(T == t)[0]) for t in TE + [300]}
res = {arm: {str(t): [0, 0] for t in TE} for arm in ("none", "static", "Tc", "T0ref")}
dS = []
for f, ms in groups.items():
    for i, j in itertools.combinations(ms, 2):
        a, b = E[i], E[j]; ed = lambda z: float(z["E_dft_eV_atom"]) * 1000; em = lambda z: float(z["E_mlip_eV_atom"]) * 1000
        dref = (ed(a) + a["F_dft"]) - (ed(b) + b["F_dft"]); d0 = (em(a) + a["F_mlip"]) - (em(b) + b["F_mlip"])
        tc, dtc = tc_anchor(T, dref, d0, 50.0)
        if tc is None or tc > 1500:
            continue
        arms = {"none": d0, "static": (ed(a) + a["F_mlip"]) - (ed(b) + b["F_mlip"]), "Tc": dtc, "T0ref": d0 - d0[it[300]] + dref[it[300]]}
        dS.append(float(np.mean(np.abs(np.gradient(d0 - dref, T)[T >= 50])) * 1000))
        for t in TE:
            k = it[t]
            if abs(dref[k]) < 1.0 or abs(t - tc) < 100:
                continue
            for arm, d in arms.items():
                res[arm][str(t)][0] += int(np.sign(d[k]) != np.sign(dref[k])); res[arm][str(t)][1] += 1
out = {arm: {t: dict(wrong=v[0] / v[1], n=v[1]) for t, v in d.items()} for arm, d in res.items()}
out["mean_abs_dS_on_subset_ueV"] = float(np.mean(dS)); out["n_pairs"] = len(dS)
Path(sys.argv[3]).write_text(json.dumps(out, indent=1) + "\n"); print(json.dumps(out, indent=1))
