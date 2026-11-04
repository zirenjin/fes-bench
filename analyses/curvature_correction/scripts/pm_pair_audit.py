"""Post-hoc audit of the RE1 random subsets: how many +-u partners are both selected.
Reproduces the subset rule of re1.py (seeded permutation of individual displacement configurations, first k kept) and
the partner rule (same atom, opposite displacement). Descriptive only; no fit is changed.
Usage: python pm_pair_audit.py <alexandria raw dir> <ENTRIES.csv> <out.json>
"""
import bz2, csv, hashlib, json, sys
from pathlib import Path
import numpy as np
import yaml

raw, entries, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
BUDGETS = (0.10, 0.25, 0.50, 1.00); REPEATS = (1, 2)
ids = [r["mp_id"] for r in csv.DictReader(open(entries))]
tot = {b: dict(conf=0, conf_with_partner_in_archive=0, conf_with_partner_selected=0) for b in BUDGETS}; n_pairs_archive = 0; n_disp = 0
for m_id in ids:
    d = yaml.load(bz2.open(raw / f"{m_id}.yaml.bz2"), Loader=yaml.CSafeLoader); disp = d["displacements"]; nd = len(disp); n_disp += nd
    partner = {i: j for i, xi in enumerate(disp) for j, xj in enumerate(disp)
               if i != j and xi["atom"] == xj["atom"] and np.allclose(np.array(xi["displacement"]), -np.array(xj["displacement"]))}
    n_pairs_archive += len(partner) // 2
    for r in REPEATS:
        perm = np.random.default_rng(int(hashlib.sha256(f"{m_id}|{r}".encode()).hexdigest()[:8], 16)).permutation(nd)
        for b in BUDGETS:
            k = nd if b == 1.0 else max(1, int(round(b * nd))); s = set(perm[:k].tolist())
            tot[b]["conf"] += len(s); tot[b]["conf_with_partner_in_archive"] += sum(1 for i in s if i in partner)
            tot[b]["conf_with_partner_selected"] += sum(1 for i in s if i in partner and partner[i] in s)
res = dict(n_entries=len(ids), n_displacements=n_disp, n_pm_pairs_in_archive=n_pairs_archive,
           fraction_of_archive_in_pm_pairs=2 * n_pairs_archive / n_disp, budgets={})
for b, t in tot.items():
    res["budgets"][str(b)] = dict(t, fraction_with_partner_selected=t["conf_with_partner_selected"] / t["conf"],
                                  fraction_of_partnered_with_partner_selected=t["conf_with_partner_selected"] / max(t["conf_with_partner_in_archive"], 1))
out.write_text(json.dumps(res, indent=1) + "\n"); print(json.dumps(res, indent=1))
