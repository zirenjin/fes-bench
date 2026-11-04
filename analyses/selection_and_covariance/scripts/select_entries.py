"""Apply the PROTOCOL.md selection rules using DFT-side information only (entry summaries from summarize_entries.py).

Rules: minimum DFT mesh frequency >= -0.05 THz; duplicates within a formula (same space group, |dE0| < 1 meV/atom,
|dV/V| < 0.5%) collapse to the lowest mp-id number; keep formulas with >= 2 distinct surviving entries.
Also requires the stored free energy to be reproduced from the stored force constants (<= 0.1 meV/atom), so that the
DFT reference used downstream is internally consistent; failures are listed, not repaired.
Usage: select_entries.py <summary.csv ...> --groups <POLYMORPH_GROUPS.csv> --out <dir>
"""
import argparse, csv, json
from collections import defaultdict
from pathlib import Path

ap = argparse.ArgumentParser(); ap.add_argument("summaries", nargs="+"); ap.add_argument("--groups"); ap.add_argument("--out")
ap.add_argument("--fmin", type=float, default=-0.05, help="stability threshold in THz (primary -0.05; D4 sensitivity -0.2)")
a = ap.parse_args(); out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
rows = {r["mp_id"]: r for f in a.summaries for r in csv.DictReader(open(f))}
group_of = {m: g["formula"] for g in csv.DictReader(open(a.groups)) for m in g["mp_ids"].split(";")}
num = lambda m: int(m.split("-")[1])
log = defaultdict(list); keep_by_formula = defaultdict(list)
for mid, formula in sorted(group_of.items(), key=lambda kv: num(kv[0])):
    r = rows.get(mid)
    if r is None:
        log["missing_summary"].append(mid); continue
    if float(r["Fvib_rebuild_max_abs_diff_meV_atom"]) > 0.1:
        log["reference_not_reproduced"].append(mid); continue
    if float(r["min_mesh_freq_THz"]) < a.fmin:
        log["dynamically_unstable"].append(mid); continue
    dup = next((k for k in keep_by_formula[formula]
                if rows[k]["spacegroup"] == r["spacegroup"]
                and abs(float(rows[k]["energy_eV_atom"]) - float(r["energy_eV_atom"])) < 1e-3
                and abs(float(rows[k]["volume_A3_atom"]) / float(r["volume_A3_atom"]) - 1) < 5e-3), None)
    if dup:
        log["duplicate_of"].append(f"{mid}->{dup}"); continue
    keep_by_formula[formula].append(mid)
final = {f: m for f, m in keep_by_formula.items() if len(m) >= 2}
dropped_formulas = sorted(f for f, m in keep_by_formula.items() if len(m) < 2)
with open(out / "FINAL_GROUPS.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["formula", "n_entries", "mp_ids"])
    for f, m in sorted(final.items(), key=lambda kv: (-len(kv[1]), kv[0])):
        w.writerow([f, len(m), ";".join(m)])
summary = dict(groups_in=len(set(group_of.values())), entries_in=len(group_of), final_groups=len(final),
               final_entries=sum(len(m) for m in final.values()), final_pairs=sum(len(m) * (len(m) - 1) // 2 for m in final.values()),
               groups_by_size={k: sum(1 for m in final.values() if len(m) == k) for k in sorted({len(m) for m in final.values()})},
               excluded={k: len(v) for k, v in log.items()}, formulas_dropped_below_two=len(dropped_formulas))
(out / "SELECTION_SUMMARY.json").write_text(json.dumps(dict(summary, details={k: v for k, v in log.items()}), indent=1) + "\n")
print(json.dumps(summary, indent=1))
