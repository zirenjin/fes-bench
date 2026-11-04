"""Step 1 of the thermodynamic error-covariance audit: find same-formula polymorph groups in the
Alexandria PBE phonon set (Loew et al. 2025) from metadata only. No phonon data is read here.

Pre-declared rules (written before any group was inspected):
- group entries by reduced formula (Materials Project OPTIMADE `chemical_formula_reduced`);
- exclude formulas containing O or F together with an element that MP-compatible settings treat with
  Hubbard U (Co, Cr, Fe, Mn, Mo, Ni, V, W): the U/no-U gauge is not recorded in the files and must
  not be mixed between polymorphs;
- keep groups with at least two entries; entries without a resolved formula are reported, not guessed;
- dynamical stability, duplicate structures and the static-energy gauge are checked later, per file.
Usage: group_polymorphs.py <ids_sizes.txt> <mp_formulas.json> <out_dir>
"""
import csv, json, re, sys
from collections import defaultdict
from pathlib import Path

ids_sizes, formulas, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
out.mkdir(parents=True, exist_ok=True)
STABLE_ONLY = True  # formulas were resolved only for entries flagged stable in pbe.csv
U_ELEMENTS = {"Co", "Cr", "Fe", "Mn", "Mo", "Ni", "V", "W"}
size = {l.split()[0]: int(l.split()[1]) for l in open(ids_sizes) if l.strip()}
form = json.load(open(formulas))
size = {k: v for k, v in size.items() if k in form or not STABLE_ONLY}


def elements(f):
    return set(re.findall(r"[A-Z][a-z]?", f))


groups, unresolved, excluded_u = defaultdict(list), [], set()
for mid in size:
    if mid not in form:
        unresolved.append(mid); continue
    f = form[mid] if isinstance(form[mid], str) else form[mid]["chemical_formula_reduced"]; el = elements(f)
    if el & U_ELEMENTS and el & {"O", "F"}:
        excluded_u.add(f); continue
    groups[f].append(mid)
poly = {f: sorted(m) for f, m in groups.items() if len(m) >= 2}
rows = [dict(formula=f, n_entries=len(m), mp_ids=";".join(m), n_elements=len(elements(f)), download_MB=round(sum(size[x] for x in m) / 1e6, 2))
        for f, m in sorted(poly.items(), key=lambda kv: (-len(kv[1]), kv[0]))]
with open(out / "POLYMORPH_GROUPS.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
summary = dict(entries=len(size), resolved=len(size) - len(unresolved), unresolved=len(unresolved), unresolved_ids=unresolved[:50],
               formulas=len(groups), excluded_U_formulas=len(excluded_u), polymorph_groups=len(poly),
               polymorph_entries=sum(len(m) for m in poly.values()),
               groups_by_size={k: sum(1 for m in poly.values() if len(m) == k) for k in sorted({len(m) for m in poly.values()})},
               download_MB_all_groups=round(sum(r["download_MB"] for r in rows), 1))
(out / "POLYMORPH_GROUPS_SUMMARY.json").write_text(json.dumps(summary, indent=1) + "\n")
print(json.dumps({k: v for k, v in summary.items() if k != "unresolved_ids"}, indent=1))
print("largest groups:", [(r["formula"], r["n_entries"]) for r in rows[:15]])
