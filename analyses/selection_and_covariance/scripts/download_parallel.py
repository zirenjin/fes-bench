"""Parallel, resume-safe download of the Alexandria PBE phonon files for the polymorph groups (6 connections).
Usage: download_parallel.py <POLYMORPH_GROUPS.csv> <out_dir> [n_threads]"""
import csv, sys, time, urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
groups, out = Path(sys.argv[1]), Path(sys.argv[2]); out.mkdir(parents=True, exist_ok=True)
NT = int(sys.argv[3]) if len(sys.argv) > 3 else 6
ids = [m for r in csv.DictReader(open(groups)) for m in r["mp_ids"].split(";")]
UA = {"User-Agent": "fes-bench-polymorph-audit/0.1 (academic research)"}


def fetch(mid):
    dst = out / f"{mid}.yaml.bz2"
    if dst.exists():
        return True
    for wait in (2, 20, 120):
        try:
            req = urllib.request.Request(f"https://alexandria.icams.rub.de/data/phonon_benchmark/pbe/{mid}.yaml.bz2", headers=UA)
            tmp = out / f"{mid}.part"; tmp.write_bytes(urllib.request.urlopen(req, timeout=300).read()); tmp.rename(dst); return True
        except Exception:
            time.sleep(wait)
    return False


with ThreadPoolExecutor(NT) as ex:
    for k, ok in enumerate(ex.map(fetch, ids), 1):
        if k % 20 == 0 or k == len(ids):
            (out / "DOWNLOAD.progress").write_text(f"{k} / {len(ids)} processed; present {sum(1 for _ in out.glob('mp-*.yaml.bz2'))}\n")
print("present", sum(1 for m in ids if (out / f"{m}.yaml.bz2").exists()), "of", len(ids))
