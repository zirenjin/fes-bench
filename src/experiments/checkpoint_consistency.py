"""Audit E_DPA/F_QH provenance and checkpoint/head consistency — plan Table 2b."""
from __future__ import annotations
import argparse,csv,json,hashlib,subprocess
from pathlib import Path
CANON_SHA='86dd3a804d78ca5d203ebf98747e8f16dff9713ba8950097ceb760b161e19907'
CANON_HEAD='Domains_Alloy'
def _load_policy(root):
    return json.loads((root/'configs/models/head_policy.yaml').read_text())['systems']

def _phase_meta(root):
    for p in sorted((root/'data/processed').glob('*/ */meta.json')):
        yield p
    # pathlib's glob does not allow a space-free recursive expression above;
    # keep this explicit so the audit only considers phase metadata.
    for p in sorted((root/'data/processed').glob('*/*/meta.json')):
        yield p

def policy_warnings(root):
    """Return human-readable warnings for phase-level head mismatches."""
    policy = _load_policy(root)
    warnings=[]
    for p in _phase_meta(root):
        try: d=json.loads(p.read_text())
        except Exception: continue
        rep=d.get('representative') or {}
        if not rep: continue
        system=p.parts[-3]
        expected=(policy.get(system) or {}).get('energy_head')
        relax=rep.get('head')
        energy=(d.get('formal_energy_provenance') or {}).get('head') or relax
        # QH provenance is stored in the canonical summary, not phase meta.
        qh_head=''
        for q in sorted((root/'result/experiments').glob(f'quasi_harmonic*/raw_runs/{system}/qh_summary.json')):
            try:
                qd=json.loads(q.read_text())
            except Exception: continue
            if qd.get('system')==system:
                qh_head=qd.get('head',''); break
        if expected and relax != expected:
            warnings.append(f'{p.relative_to(root)} relaxation head={relax!r}, policy={expected!r}')
        if expected and energy != expected:
            warnings.append(f'{p.relative_to(root)} E_DPA head={energy!r}, policy={expected!r}')
        if qh_head and expected and qh_head != expected:
            warnings.append(f'{p.relative_to(root)} QH head={qh_head!r}, policy={expected!r}')
        if len({h for h in (relax,energy,qh_head) if h}) > 1:
            warnings.append(f'{p.relative_to(root)} relaxation/E_DPA/QH heads disagree: {relax!r}, {energy!r}, {qh_head!r}')
    # Comparison tables carrying ΔE/ΔF values must identify both sides with
    # the same immutable model provenance.
    for p in sorted((root/'result/experiments').rglob('*.csv')):
        try:
            import csv as _csv
            with p.open(newline='', encoding='utf-8') as fh:
                reader=_csv.DictReader(fh); fields=set(reader.fieldnames or [])
                if not {'old_checkpoint_sha256','new_checkpoint_sha256','old_head','new_head'} <= fields:
                    continue
                for i,row in enumerate(reader, start=2):
                    if row['old_checkpoint_sha256'] != row['new_checkpoint_sha256'] or row['old_head'] != row['new_head']:
                        warnings.append(f'{p.relative_to(root)}:{i} comparison sides use different checkpoint/head')
        except (OSError, UnicodeError):
            continue
    return warnings

def main():
 p=argparse.ArgumentParser(); p.add_argument('--repo-root',type=Path,default=Path('.')); a=p.parse_args(); root=a.repo_root.resolve(); rows=[]
 policy=_load_policy(root)
 for path in sorted(list((root/'data/processed').rglob('*.json'))+list((root/'result/experiments').rglob('*.json'))):
  try:d=json.loads(path.read_text())
  except Exception:continue
  txt=json.dumps(d)
  if not any(k in txt for k in ('E_DPA','F_QH','old_E_DPA','new_E_DPA','formal_energy_eV_per_atom')): continue
  prov=d.get('provenance',{}) if isinstance(d,dict) else {}; head=prov.get('head',d.get('head','') if isinstance(d,dict) else ''); sha=prov.get('checkpoint_sha256',d.get('checkpoint_sha256','') if isinstance(d,dict) else '')
  if not head and isinstance(d,dict): head=d.get('representative',{}).get('head','') or d.get('formal_energy_provenance',{}).get('head','')
  if not sha and isinstance(d,dict): sha=d.get('formal_energy_provenance',{}).get('checkpoint_sha256','')
  rel=str(path.relative_to(root)); system = rel.split('/')[2] if rel.startswith('data/processed/') and len(rel.split('/')) > 2 else None
  expected_head=(policy.get(system) or {}).get('energy_head', CANON_HEAD)
  status = 'ok' if sha==CANON_SHA and head==expected_head else ('missing_provenance' if not sha or not head else 'mismatch')
  rows.append({'artifact':rel,'checkpoint_sha256':sha or '','head':head or '','canonical_checkpoint':(sha==CANON_SHA or not sha),'canonical_head':(head==expected_head or not head),'status':status})
 out=root/'result/experiments/checkpoint_consistency'; out.mkdir(parents=True,exist_ok=True)
 with (out/'findings.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=['artifact','checkpoint_sha256','head','canonical_checkpoint','canonical_head','status']); w.writeheader(); w.writerows(rows)
 pw=policy_warnings(root)
 meta={'git_commit':subprocess.check_output(['git','-C',str(root),'rev-parse','HEAD'],text=True).strip(),'canonical_checkpoint_sha256':CANON_SHA,'canonical_heads':{k:v.get('energy_head') for k,v in policy.items()},'rows':len(rows),'mismatches':sum(r['status']=='mismatch' for r in rows),'missing':sum(r['status']=='missing_provenance' for r in rows),'policy_warnings':pw}
 (out/'findings.meta.json').write_text(json.dumps(meta,indent=2)+'\n')
 for w in pw: print('⚠ '+w)
 print(json.dumps(meta)); return 0
if __name__=='__main__': raise SystemExit(main())
