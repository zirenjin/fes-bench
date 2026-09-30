"""Attach formal E_DPA provenance to SiO2 phase metadata — plan Table 2b."""
import argparse,json
from pathlib import Path
def main():
 p=argparse.ArgumentParser(); p.add_argument('--repo-root',type=Path,default=Path('.')); p.add_argument('--energies',required=True); a=p.parse_args(); root=a.repo_root.resolve(); d=json.loads((root/a.energies).read_text()); by={r['phase']:r for r in d['rows']}
 for phase,row in by.items():
  path=root/'data/processed/sio2'/phase/'meta.json'; m=json.loads(path.read_text()); m['formal_energy_eV_per_atom']=row['E_DPA_eV_per_atom']; m['formal_energy_provenance']={'checkpoint':row['checkpoint'],'checkpoint_sha256':row['checkpoint_sha256'],'head':row['head'],'source':'result/experiments/formal_sio2_energies/energies.json','structure_sha256':row['structure_sha256']}; path.write_text(json.dumps(m,indent=2,ensure_ascii=False)+'\n')
if __name__=='__main__': main()
