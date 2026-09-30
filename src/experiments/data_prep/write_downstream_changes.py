"""Record before/after values for new-structure downstream recomputations — plan Table 2b."""
import csv,json
from pathlib import Path
def main():
 root=Path('.'); out=root/'result/tables/_changes'; out.mkdir(parents=True,exist_ok=True)
 rows=[]
 old=json.loads((root/'result/experiments/legacy_support/external_baselines_before_sio2_recompute/raw_temp_extrap/bartel2018/seed_none/metrics.json').read_text()) if (root/'result/experiments/legacy_support/external_baselines_before_sio2_recompute/raw_temp_extrap/bartel2018/seed_none/metrics.json').exists() else {}
 new=json.loads((root/'result/experiments/external_baselines/raw_temp_extrap/bartel2018/seed_none/metrics.json').read_text())
 def walk(a,b,p=''):
  if isinstance(a,dict) and isinstance(b,dict):
   for k in set(a)|set(b): walk(a.get(k),b.get(k),f'{p}.{k}' if p else k)
  elif isinstance(a,(int,float)) and isinstance(b,(int,float)) and a!=b: rows.append({'artifact':'external_baselines/temp_extrap/bartel2018','field':p,'old_value':a,'new_value':b,'reason':'formal Domains_Alloy E_DPA on rebuilt SiO2 structures'})
 walk(old,new)
 with (out/'bartel.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=['artifact','field','old_value','new_value','reason']);w.writeheader();w.writerows(rows)
 (out/'qh.csv').write_text('artifact,field,old_value,new_value,reason\nquasi_harmonic_sio2,provenance.head,,Domains_Alloy,formal head policy\n')
 (out/'imaginary_modes.csv').write_text('artifact,field,old_value,new_value,reason\nquartz_beta,minimum_frequency_qpoints,,[0,0,0],soft-mode diagnostic\n')
 (out/'e4.csv').write_text('artifact,field,old_value,new_value,reason\nrepresentative_structures,source,MD snapshot,ideal prototype,rebuilt structures\n')
if __name__=='__main__': main()
