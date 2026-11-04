"""Resolve reduced formulas for Alexandria PBE phonon entries flagged stable in pbe.csv, via MP OPTIMADE.
Polite: 20 ids per request, 1 s pause, incremental save and resume. Progress in <out>.progress.
Usage: resolve_formulas.py <pbe.csv> <out.json>"""
import csv, json, sys, time, urllib.parse, urllib.request
from pathlib import Path
src, out = Path(sys.argv[1]), Path(sys.argv[2])
ids = [r["mp_id"] for r in csv.DictReader(open(src)) if r["stable"] == "True"]
done = json.load(open(out)) if out.exists() else {}
todo = [i for i in ids if i not in done]
for k in range(0, len(todo), 20):
    chunk = todo[k:k + 20]
    q = urllib.parse.urlencode({"filter": " OR ".join(f'id="{i}"' for i in chunk), "response_fields": "chemical_formula_reduced", "page_limit": 50})
    for wait in (1, 10, 60):
        try:
            req = urllib.request.Request("https://optimade.materialsproject.org/v1/structures?" + q, headers={"User-Agent": "fes-bench-polymorph-audit/0.1 (academic research; metadata only)"})
            with urllib.request.urlopen(req, timeout=60) as r:
                for d in json.load(r)["data"]: done[d["id"]] = d["attributes"]["chemical_formula_reduced"]
            break
        except Exception as e:
            Path(str(out) + ".progress").write_text(f"error {e} at batch {k // 20}; waiting {wait}s\n"); time.sleep(wait)
    json.dump(done, open(out, "w"))
    Path(str(out) + ".progress").write_text(f"{len(done)} / {len(ids)} resolved, batch {k // 20 + 1} / {(len(todo) + 19) // 20}\n")
    time.sleep(1)
print("resolved", len(done), "of", len(ids), "; unresolved:", [i for i in ids if i not in done][:20])
