"""Download the Alexandria PBE phonon files (CC BY 4.0, Loew et al. 2025) for the polymorph groups.
Sequential, short pause, resume-safe; progress in <out_dir>/DOWNLOAD.progress.
Usage: download_groups.py <POLYMORPH_GROUPS.csv> <out_dir> [reverse]
With 'reverse', the list is processed from the end (used for a second, parallel downloader)."""
import csv, sys, time, urllib.request
from pathlib import Path
groups, out = Path(sys.argv[1]), Path(sys.argv[2]); out.mkdir(parents=True, exist_ok=True)
ids = [m for r in csv.DictReader(open(groups)) for m in r["mp_ids"].split(";")]
if len(sys.argv) > 3 and sys.argv[3] == "reverse":
    ids = ids[::-1]
UA = {"User-Agent": "fes-bench-polymorph-audit/0.1 (academic research)"}
for k, mid in enumerate(ids, 1):
    dst = out / f"{mid}.yaml.bz2"
    if not dst.exists():
        for wait in (2, 20, 120):
            try:
                req = urllib.request.Request(f"https://alexandria.icams.rub.de/data/phonon_benchmark/pbe/{mid}.yaml.bz2", headers=UA)
                tmp = dst.with_suffix(".part"); tmp.write_bytes(urllib.request.urlopen(req, timeout=120).read()); tmp.rename(dst); break
            except Exception as e:
                (out / "DOWNLOAD.progress").write_text(f"error {e} on {mid}; retry in {wait}s\n"); time.sleep(wait)
        time.sleep(0.2)
    if k % 20 == 0 or k == len(ids): (out / "DOWNLOAD.progress").write_text(f"{k} / {len(ids)}\n")
print("done", sum(1 for m in ids if (out / f"{m}.yaml.bz2").exists()), "of", len(ids))
