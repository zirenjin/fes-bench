"""Descriptive addendum (post hoc, requested in review): descriptor baseline and DPA-3.1 on the same pairs, i.e. pairs
valid for both. Usage: python analyses/descriptor_baseline/common_subset.py   (from the repository root)
"""
import json, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).parent; DB = HERE.parents[1] / "database_pairs"; sys.path.insert(0, str(DB))
from score import Bench, GATE  # noqa: E402

B = Bench(); res = {}
KEYS = ("near_Tc_rms_e", "pair_dFth_rms_1000K", "pair_dS_rms_1000K", "wrong_ordering_mean", "crossing_hit50", "near_Tc_within_energy_budget50")
for S in ("S1", "S0"):
    fb, fd = HERE / f"predictions/bartel_{S}.npz", DB / f"predictions/F0_DPA-3.1-3M_{S}.npz"
    vb, vd = np.load(fb)["fmin"], np.load(fd)["fmin"]
    mask = B.ref_valid & (vb[B.pi] >= GATE) & (vb[B.pj] >= GATE) & (vd[B.pi] >= GATE) & (vd[B.pj] >= GATE)
    blk = dict(n_pairs=int(mask.sum()), n_crossing_pairs=int((mask & B.cross).sum()))
    for name, f in (("descriptor", fb), ("DPA-3.1", fd)):
        m = B.score(f, mask)["metrics"]; blk[name] = {k: dict(value=m[k]["value"], ci95=m[k]["ci95"]) for k in KEYS}
    res[S] = blk
(HERE / "results/BARTEL_COMMON_SUBSET.json").write_text(json.dumps(res, indent=1) + "\n")
for S, b in res.items():
    print(S, b["n_pairs"], b["n_crossing_pairs"], {k: (round(b["descriptor"][k]["value"], 3), round(b["DPA-3.1"][k]["value"], 3)) for k in KEYS})
