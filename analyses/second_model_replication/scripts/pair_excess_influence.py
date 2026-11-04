"""Which formulas carry MatterSim's pooled pair dF_th RMS excess over F0 (post hoc, descriptive influence diagnostic).

Two drop orders are reported:
- by the excess itself (MS - F0 sum of squares; selected on the contrast, so an influence check only);
- by the larger of the two models' per-formula sums of squares (symmetric; does not select on the sign of the contrast).
Usage: pair_excess_influence.py <F0_cov_dir> <MS_cov_dir> <FMIN_ENTRIES.csv> > out.txt
"""
import csv, sys
from collections import defaultdict
import numpy as np

f0_dir, ms_dir, fmin_csv = sys.argv[1:4]
load = lambda d: {(r["i"], r["j"]): r for r in csv.DictReader(open(d + "/PAIRS.csv"))}
F, M = load(f0_dir), load(ms_dir)
fmin = {r["mp_id"]: r for r in csv.DictReader(open(fmin_csv))}
gate_min = lambda m: min(float(fmin[m][c]) for c in ("fmin_dft", "fmin_f0", "fmin_ms"))
for t in (300, 600, 1000, 1500):
    b = f"b{t}"; per = defaultdict(lambda: [0.0, 0.0, 0, 0.0])
    for k in F:
        v = per[F[k]["formula"]]; v[0] += float(F[k][b]) ** 2; v[1] += float(M[k][b]) ** 2; v[2] += 1; v[3] = min(v[3], gate_min(k[0]), gate_min(k[1]))
    s0 = sum(v[0] for v in per.values()); s1 = sum(v[1] for v in per.values()); ex = s1 - s0
    print(f"T = {t} K: pooled pair sum of squares F0 {s0:.0f}, MS {s1:.0f}, excess {ex:.0f}")
    by_ex = sorted(per, key=lambda f: -(per[f][1] - per[f][0])); cum = 0.0
    for n, f in enumerate(by_ex[:8], 1):
        v = per[f]; cum += v[1] - v[0]
        print(f"  {n} {f:10s} pairs {v[2]:2d}  pair RMS F0 {np.sqrt(v[0] / v[2]):6.2f}  MS {np.sqrt(v[1] / v[2]):6.2f}  "
              f"cumulative share of excess {cum / ex:4.2f}  min fmin (DFT/F0/MS) {v[3]:.2f} THz")
    def ratio(keep):
        return np.sqrt(sum(per[f][1] for f in keep) / sum(per[f][0] for f in keep))
    print("  pair RMS ratio MS/F0 after dropping the top-k by excess: " + ", ".join(f"k={k}: {ratio(by_ex[k:]):.3f}" for k in (0, 1, 2, 3, 5, 10)))
    by_max = sorted(per, key=lambda f: -max(per[f][0], per[f][1]))
    print("  after dropping the k largest-error formulas (either model): " + ", ".join(f"k={k}: {ratio(by_max[k:]):.3f}" for k in (1, 3, 5, 10)))
    r = [np.sqrt(v[1] / v[0]) for v in per.values() if v[0] > 0]
    print(f"  per-formula pair RMS ratio MS/F0: median {np.median(r):.3f}; MS worse in {np.mean([x > 1 for x in r]):.2f} of {len(r)} formulas")
