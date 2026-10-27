#!/usr/bin/env python3
"""Audit the historical Hf sign convention against original phase-resolved files.

Figure contract: explain the large cancellation, without inferring a physical cause.
Quantitative 2x2 grid: phase identity -> signed pair -> residual -> decomposition.
All 1252 temperatures retained. One historical seed/checkpoint; no error bars.
"""
from pathlib import Path
import csv,json,hashlib,shutil,sys
import numpy as np
R=Path(__file__).resolve().parents[1];W=R.parent.parent
E=Path('[local archive path]')
I=R/'data/hf_sign_audit';I.mkdir(exist_ok=True)
orig={}
for ph in ['bcc','hcp']:
 for n in ['fqh.csv','reference_G.csv','structure.extxyz','meta.json','phonon_report.json']:
  orig[ph+'_'+n]=W/f'outputs/fes-bench-decisive-hf-20261019/inputs/data/processed/hf/{ph}/{n}'
 orig[ph+'_raw_fqh.csv']=E/f'result/experiments/legacy_support/phase2_domains_alloy/hf/{ph}_fqh.csv'
 orig[ph+'_raw_reference.tab']=E/f'data/raw/hf/table-thermodynamic-Hf-{ph}-PBE.tab'
for n in ['run.py','materialize.py']:orig['qh_'+n]=E/'src/lib/fes_bench/qh'/n
orig['import_thermo.py']=E/'src/lib/fes_bench/data/import_thermo.py'
orig['raw_qh_summary.json']=E/'result/experiments/legacy_support/phase2_domains_alloy/hf/qh_summary.json'
provenance={}
for n,p in orig.items():
 if p.is_file():shutil.copyfile(p,I/n)
 assert (I/n).is_file(),p
 provenance[n]={'original':str(p),'sha256':hashlib.sha256((I/n).read_bytes()).hexdigest()}
def read(p):
 with p.open(encoding='utf-8-sig') as f:return list(csv.DictReader(f))
def index(rs):return {float(r['T_K']):r for r in rs}
fq={ph:index(read(I/(ph+'_fqh.csv'))) for ph in ['bcc','hcp']}
g={ph:index(read(I/(ph+'_reference_G.csv'))) for ph in ['bcc','hcp']}
raw={ph:index(read(I/(ph+'_raw_fqh.csv'))) for ph in ['bcc','hcp']}
assert all(float(r['volume_scale'])==1.02 for ph in raw for r in raw[ph].values())
refraw={}
for ph in ['bcc','hcp']:
 with (I/(ph+'_raw_reference.tab')).open() as f:
  next(f);refraw[ph]={float(r['# T (K)']):float(r['Gibbs_energy'])*.001 for r in csv.DictReader(f,delimiter='\t')}
D=json.loads((R/'data/review_data.json').read_text());comp=index(D['hfComponent'])
assert set(fq['bcc'])==set(fq['hcp'])==set(g['bcc'])==set(g['hcp'])==set(comp)
checks={};points=[]
for t in sorted(comp):
 q={ph:float(fq[ph][t]['F_QH_eV_per_atom']) for ph in fq};v={ph:float(fq[ph][t]['F_vib_eV_per_atom']) for ph in fq};e={ph:float(raw[ph][t]['E_static_eV_per_atom']) for ph in fq};gr={ph:float(g[ph][t]['G_eV_per_atom']) for ph in fq}
 dq=q['bcc']-q['hcp'];dg=gr['bcc']-gr['hcp'];dv=v['bcc']-v['hcp'];de=e['bcc']-e['hcp'];c=comp[t]
 points.append(dict(T_K=t,domain=c['domain'],F_QH_bcc_eV_atom=q['bcc'],F_QH_hcp_eV_atom=q['hcp'],delta_G_ref_eV_atom=dg,delta_F_QH_eV_atom=dq,negative_delta_F_QH_eV_atom=-dq,R_plus_eV_atom=dg-dq,R_minus_eV_atom=dg+dq,delta_static_eV_atom=de,delta_F_vib_eV_atom=dv,learned_correction_eV_atom=float(c['delta_correction_bcc_minus_hcp_eV_atom'])))
 for name,val in [('pair_QH_vs_archive',dq-float(c['delta_FQH_bcc_minus_hcp_eV_atom'])),('pair_ref_vs_archive',dg-float(c['deltaG_ref_bcc_minus_hcp_eV_atom'])),('target_identity',dg-dq-float(c['delta_R_target_bcc_minus_hcp_eV_atom'])),('QH_static_plus_vib',dq-de-dv),('learned_reconstruction',dq+float(c['delta_correction_bcc_minus_hcp_eV_atom'])-float(c['deltaG_T2_eV_atom']))]:checks[name]=max(checks.get(name,0),abs(val))
 for ph in fq:
  for name,val in [('canonical_vs_raw_QH',q[ph]-float(raw[ph][t]['F_QH_eV_per_atom'])),('reference_vs_raw_table',gr[ph]-refraw[ph][t])]:checks[name]=max(checks.get(name,0),abs(val))
assert max(checks.values())<1e-12,checks
fields=list(points[0]);csvpath=R/'figures/source_data/hf_qh_sign_audit.csv'
with csvpath.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(points)
summary={'n_paired_temperatures':len(points),'checks_max_abs_eV_atom':checks,'at_1400_K':next(r for r in points if r['T_K']==1400),'static_gap_meV_atom':1000*points[0]['delta_static_eV_atom'],'MAE_gap_meV_atom':float(np.mean([abs(r['delta_F_QH_eV_atom']-r['delta_G_ref_eV_atom'])*1000 for r in points])),'counterfactual_reversed_MAE_meV_atom':float(np.mean([abs(-r['delta_F_QH_eV_atom']-r['delta_G_ref_eV_atom'])*1000 for r in points])),'input_provenance':provenance}
(I/'AUDIT.json').write_text(json.dumps(summary,indent=2)+'\n')
# Reuse the established Python style and measured alignment gate.
import plot_evidence as p
p.SELECT.clear();plt=p.plt
fig,axs=plt.subplots(2,2,figsize=(p.W,7.3));fig.subplots_adjust(left=.13,right=.97,bottom=.10,top=.85,wspace=.40,hspace=.85)
T=np.array([r['T_K'] for r in points])
def line(ax,k,name,color,ls='-',scale=1000):ax.plot(T,[r[k]*scale for r in points],label=name,color=color,ls=ls)
a,b,c,d=axs.flat
line(a,'F_QH_bcc_eV_atom','bcc',p.C['T2'],scale=1);line(a,'F_QH_hcp_eV_atom','hcp',p.C['T3'],scale=1);a.set_ylabel('F_QH (eV/atom)')
line(b,'delta_G_ref_eV_atom','Reference ΔG',p.C['ink']);line(b,'delta_F_QH_eV_atom','ΔF_QH',p.C['QH'],'--');line(b,'negative_delta_F_QH_eV_atom','−ΔF_QH (sign test)',p.C['residual'],':');b.set_ylabel('Pair energy (meV/atom)')
line(c,'R_plus_eV_atom','R+ = ΔG − ΔF_QH',p.C['T3']);line(c,'R_minus_eV_atom','R− = ΔG + ΔF_QH',p.C['residual'],'--');line(c,'learned_correction_eV_atom','Historical correction',p.C['T2']);c.set_ylabel('Pair correction (meV/atom)')
line(d,'delta_F_QH_eV_atom','Total ΔF_QH',p.C['ink']);line(d,'delta_static_eV_atom','ΔE_static',p.C['QH'],'--');line(d,'delta_F_vib_eV_atom','ΔF_vib',p.C['T3']);d.set_ylabel('Pair energy (meV/atom)')
for ax,letter,title in zip(axs.flat,'abcd',['Original phase curves','Independent bcc − hcp','Residual sign test','Static and vibrational parts']):
 p.panel(ax,letter,title);ax.set_xlabel('Temperature (K)');p.clean(ax,False);ax.legend(loc='lower left',bbox_to_anchor=(-.02,1.12),frameon=False,ncol=1,fontsize=7)
p.export('hf_qh_sign_audit','Hf: phase-order and residual audit',fig,points,'All 1252 paired temperatures. Historical seed 11 / step 1500. QH denotes the legacy diagnostic curve, not validated QHA. R− reverses the baseline only as a counterfactual. No reference alignment or parameters are fitted.','diagnostics')
# Add to the existing supplemental figure registry without changing other assets.
path=R/'figures/explainer_figures.json';m=json.loads(path.read_text());m['figures']=[x for x in m['figures'] if x['key']!='hf_qh_sign_audit']+p.FIGS
path.write_text(json.dumps(m,ensure_ascii=False,indent=2));(R/'figures/explainer_figures.js').write_text('window.FES_EXPLAINER_FIGURES='+json.dumps(m,ensure_ascii=False)+';\n')
r=summary['at_1400_K'];report=f'''# Hf historical QH sign and convention audit

## Result

No bcc/hcp subtraction reversal was found from the phase-resolved source exports through the plotted pairs. All {len(points)} temperatures were checked. Both reference and QH use bcc minus hcp. The target is ΔR = ΔG_ref − ΔF_QH. Residual cancellation is real in these data. It is not evidence that the QH curve is physically valid.

At 1400 K (meV/atom):

| Quantity | Value |
| --- | ---: |
| Reference ΔG | {r['delta_G_ref_eV_atom']*1000:.6f} |
| ΔF_QH | {r['delta_F_QH_eV_atom']*1000:.6f} |
| Target ΔR | {r['R_plus_eV_atom']*1000:.6f} |
| ΔE_static | {r['delta_static_eV_atom']*1000:.6f} |
| ΔF_vib | {r['delta_F_vib_eV_atom']*1000:.6f} |

The static bcc−hcp gap is about 152.03 meV/atom; the thermal part adds another positive gap. A sign flip is not a legal repair. Across this full support, baseline gap MAE is {summary['MAE_gap_meV_atom']:.3f} meV/atom; reversing it gives {summary['counterfactual_reversed_MAE_meV_atom']:.3f} meV/atom. These are descriptive sign checks, not held-out model scores.

## Four checks

1. Independently subtract the two original QH exports. They match canonical QH and the plotted pair to numerical precision. No phase-specific sign transform occurs in materialization.
2. Independently import the original phase-specific reference tables (meV/atom × 0.001). They match canonical reference values at the common temperatures. The import applies no phase-specific additive alignment. The reference sources and structure metadata identify bcc and hcp separately. Structure metadata reports Im-3m (one-atom primitive cell) and P6_3/mmc (two atoms), respectively. This audit does not rerun symmetry finding or the original phonon force calculations.
3. Verify ΔR = ΔG_ref − ΔF_QH and ΔG_pred = ΔF_QH + Δcorrection. Maximum absolute check errors are in `data/hf_sign_audit/AUDIT.json`; all are below 1e-12 eV/atom.
4. Inspect the baseline construction. The actual canonical arrays match the early 2×2×2, three-volume prototype under legacy support, not the later enlarged-supercell run. F_QH = E_static + F_vib, per atom. Static energies are −7.0414142013 (bcc) and −7.1934467107 (hcp) eV/atom, constant along this selected-volume curve. Both select the upper sampled volume scale 1.02 throughout the temperature support. This is a discrete-volume diagnostic, not a converged EOS minimization.

## Energy content and limits

The archived QH runner adds the potential static energy and Phonopy thermal free energy. The harmonic thermal convention includes zero-point energy; this audit checks the implementation path but does not recompute its phonon integral or isolate a numerical ZPE from the high-temperature table. There is no explicit PV term or pressure-dependent minimization in this runner; materialization inherits the reference P = 0.0001 GPa column. At the stored volumes PV is approximately 0.014 meV/atom, far below the observed hundreds of meV/atom gap. Electronic and magnetic thermal terms are not explicitly included in this baseline.

The reference and model have different absolute energy zeros. The historical reconstruction applies one train-only scalar offset, common to both phases; it cancels in ΔG. No additional phase-specific offset is introduced here. Common-offset invariance cannot establish that the model static phase gap is physically accurate, or independently certify the original reference publisher's cross-phase energy convention.

The bcc legacy phonon report has negative frequencies, including about −1.963 THz at the selected volume. It explicitly marks this baseline diagnostic_only=true and qh_reliable=false. Ordinary harmonic thermodynamics is not valid for an imaginary oscillator. The exact historical mode filtering/weighting needs its original runtime or a mode-resolved cache to be reconstructed; it cannot be inferred from this aggregate CSV alone. Therefore residuals here must not be called pure anharmonic, electronic or magnetic corrections.

## Interpretation

The evidence supports a large mismatch between a diagnostic numerical baseline and the reference, plus failure of this historical learned correction to reproduce the needed cancellation. It does not identify which physical contribution causes the mismatch. The baseline generator uses a DPA alloy head; the separate historical residual-model initializer had a Si/O mapping problem. These are distinct checkpoints and should not be conflated.

[Four-panel audit figure](figures/hf_qh_sign_audit.pdf) · [All paired values](figures/source_data/hf_qh_sign_audit.csv) · [Checks and input hashes](data/hf_sign_audit/AUDIT.json).
'''
(R/'HF_QH_SIGN_AUDIT.md').write_text(report)
print(json.dumps({k:v for k,v in summary.items() if k!='input_provenance'},indent=2))
