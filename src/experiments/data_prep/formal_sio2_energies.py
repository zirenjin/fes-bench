"""Compute canonical SiO2 E_DPA inputs for downstream baselines — plan Table 2b."""
from __future__ import annotations
import argparse, hashlib, json, subprocess
from pathlib import Path
from ase.io import read

CHECKPOINT_SHA = "86dd3a804d78ca5d203ebf98747e8f16dff9713ba8950097ceb760b161e19907"

def sha256(path: Path) -> str:
    h=hashlib.sha256(); h.update(path.read_bytes()); return h.hexdigest()

def main() -> int:
    p=argparse.ArgumentParser(); p.add_argument('--repo-root',type=Path,default=Path('.')); p.add_argument('--checkpoint',required=True); p.add_argument('--head',default='Domains_Alloy'); p.add_argument('--output',default='result/experiments/formal_sio2_energies/energies.json'); a=p.parse_args()
    from deepmd.calculator import DP
    root=a.repo_root.resolve(); calc=DP(model=a.checkpoint,head=a.head); rows=[]
    for phase in ('quartz_beta','cristobalite_beta','tridymite_p63mmc'):
        path=root/'data/processed/sio2'/phase/'structure.extxyz'; atoms=read(path); atoms.calc=calc; total=float(atoms.get_potential_energy()); rows.append({'phase':phase,'structure':str(path.relative_to(root)),'structure_sha256':sha256(path),'n_atoms':len(atoms),'E_DPA_eV_per_atom':total/len(atoms),'checkpoint':'<external-checkpoint>/DPA-3.1-3M.pt','checkpoint_sha256':CHECKPOINT_SHA,'head':a.head})
    try: commit=subprocess.check_output(['git','-C',str(root),'rev-parse','HEAD'],text=True).strip()
    except Exception: commit='unavailable_external_workspace'
    out=root/a.output; out.parent.mkdir(parents=True,exist_ok=True); out.write_text(json.dumps({'provenance':{'git_commit':commit,'checkpoint_sha256':CHECKPOINT_SHA,'checkpoint':'<external-checkpoint>/DPA-3.1-3M.pt','head':a.head},'rows':rows},indent=2)+'\n'); print(out)
    return 0
if __name__=='__main__': raise SystemExit(main())
