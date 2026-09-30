"""Audit E_DPA/F_QH provenance and comparison-side consistency — plan Table 2b."""
from __future__ import annotations
import argparse,csv,json,hashlib,subprocess
from pathlib import Path
CANON_SHA='86dd3a804d78ca5d203ebf98747e8f16dff9713ba8950097ceb760b161e19907'
CANON_HEAD='Domains_Alloy'
def main():
 p=argparse.ArgumentParser(); p.add_argument('--repo-root',type=Path,default=Path('.')); a=p.parse_args(); root=a.repo_root.resolve(); rows=[]
 for path in sorted(list((root/'data/processed').rglob('*.json'))+list((root/'result/experiments').rglob('*.json'))):
  try:d=json.loads(path.read_text())
  except Exception:continue
  txt=json.dumps(d)
  if not any(k in txt for k in ('E_DPA','F_QH','old_E_DPA','new_E_DPA','formal_energy_eV_per_atom')): continue
  prov=d.get('provenance',{}) if isinstance(d,dict) else {}; head=prov.get('head',d.get('head','') if isinstance(d,dict) else ''); sha=prov.get('checkpoint_sha256',d.get('checkpoint_sha256','') if isinstance(d,dict) else '')
  if not head and isinstance(d,dict): head=d.get('representative',{}).get('head','') or d.get('formal_energy_provenance',{}).get('head','')
  if not sha and isinstance(d,dict): sha=d.get('formal_energy_provenance',{}).get('checkpoint_sha256','')
  rel=str(path.relative_to(root)); expected_relax = rel.startswith('data/processed/sio2/') and isinstance(d,dict) and 'formal_energy_provenance' in d
  status = 'expected_noncanonical_relaxation' if expected_relax else ('ok' if sha==CANON_SHA and head==CANON_HEAD else ('missing_provenance' if not sha or not head else 'mismatch'))
  rows.append({'artifact':rel,'checkpoint_sha256':sha or '','head':head or '','canonical_checkpoint':(sha==CANON_SHA or not sha),'canonical_head':(head==CANON_HEAD or not head),'status':status})
 out=root/'result/experiments/checkpoint_consistency'; out.mkdir(parents=True,exist_ok=True)
 with (out/'findings.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=['artifact','checkpoint_sha256','head','canonical_checkpoint','canonical_head','status']); w.writeheader(); w.writerows(rows)
 meta={'git_commit':subprocess.check_output(['git','-C',str(root),'rev-parse','HEAD'],text=True).strip(),'canonical_checkpoint_sha256':CANON_SHA,'canonical_head':CANON_HEAD,'rows':len(rows),'mismatches':sum(r['status']=='mismatch' for r in rows),'missing':sum(r['status']=='missing_provenance' for r in rows)}
 (out/'findings.meta.json').write_text(json.dumps(meta,indent=2)+'\n')
 print(json.dumps(meta)); return 0
if __name__=='__main__': raise SystemExit(main())
