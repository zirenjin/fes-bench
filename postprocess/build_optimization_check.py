#!/usr/bin/env python3
"""Frozen-prediction fit audit. No model refit, target alignment or layer probe.

Figure contract: compare source-fit and held-out error distributions in one
protocol. Both panels retain all paired temperature observations and seeds.
Seven independent views; temperatures are repeated measurements, not materials.
"""
from pathlib import Path
from collections import defaultdict,Counter
import csv,json,hashlib,itertools
import numpy as np
import plot_evidence as p
R=p.ROOT;W=R.parent.parent;D=p.D
ROOT=W/'outputs';CROOT=ROOT/'fes-bench-cross-system-20261020/results'
SPECS=[('core_learned',CROOT/'learned_models/T2_T3_combined_point_predictions.csv','cross_system_60_20_20_batch8_3000'),('core_controls',CROOT/'remote_baselines/BASELINE_POINT_PREDICTIONS.csv','cross_system_60_20_20_batch8_3000'),('sio2_matched',CROOT/'T2_SiO2_QH_residual_20261020/MATCHED_COMPARISON_POINT_PREDICTIONS.csv','cross_system_60_20_20_batch8_3000'),('hf_multiseed',ROOT/'fes-bench-decisive-hf-20261019/MULTISEED_POINT_PREDICTIONS.csv','hf_60_20_20_batch2_3000'),('hf_controls',ROOT/'fes-bench-decisive-hf-20261019/POINT_PREDICTIONS.csv','hf_60_20_20_batch2_3000'),('hf_long',ROOT/'fes-bench-t2-hf-hfo2-init-contrast-20261020/results/HF_BATCH8_LONG/SELECTED_STEP_8000_ALL_REGIONS.csv','hf_60_20_20_batch8_10000_selected8000')]
def read(path):
 with path.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def frozen(path):
 sha=hashlib.sha256(path.read_bytes()).hexdigest();s=next(s for s in D['sources'] if s['sha256']==sha);return read(R/s['href']),s
phase=defaultdict(dict);direct={};provenance=[];skip=Counter()
for tag,path,protocol in SPECS:
 rows,src=frozen(path);provenance.append(dict(tag=tag,id=src['id'],sha256=src['sha256'],rows=len(rows)))
 for r in rows:
  if tag.startswith('core') and r['system'] not in ['ti','zr']:skip['earlier_sio2_superseded_by_matched_run']+=1;continue
  if tag=='hf_controls' and r['model'] not in ['thermo_form_3term','same_form_4term_analytic']:skip['hf_models_covered_by_multiseed_export']+=1;continue
  model=r.get('model') or r['method'];sys=r['system'];seed=r.get('seed','none') or 'none';domain=r['region'];ph=r['phase'];t=float(r.get('T_K') or r['temperature_K'])
  if tag=='hf_long':model='T2_Hf_checkpoint_8000'
  base=(protocol,sys,model,seed,domain)
  if ph in ['bcc_minus_hcp','bcc-hcp']:
   y=r.get('delta_G_reference_eV_per_atom') or r.get('delta_G_reference_eV_atom');hat=r.get('delta_G_prediction_eV_per_atom') or r.get('delta_G_pred_eV_atom')
   if y and hat:direct[base+('bcc−hcp',t)]=(float(y),float(hat),src['id'])
   continue
  y=r.get('G_ref_eV_atom') or r.get('G_reference_eV_per_atom')
  hat=r.get('G_pred_train_offset_eV_atom') or r.get('G_pred_after_train_offset_eV_atom') or r.get('G_prediction_eV_per_atom_train_calibrated') or r.get('G_prediction_eV_per_atom')
  if y is None or hat in [None,'']:skip['missing_phase_output']+=1;continue
  key=base+(ph,);assert t not in phase[key];phase[key][t]=(float(y),float(hat),src['id'])
phases=defaultdict(dict)
for k,v in phase.items():phases[k[:-1]][k[-1]]=v
out=[]
def add(case,protocol,system,model,seed,domain,pair,t,y,hat,sid):
 out.append(dict(case=case,protocol=protocol,system=system,model=model,seed=seed,domain=domain,pair=pair,T_K=t,reference_eV_atom=y,prediction_eV_atom=hat,error_meV_atom=(hat-y)*1000,absolute_error_meV_atom=abs(hat-y)*1000,source_id=sid))
for base,ps in phases.items():
 protocol,sys,model,seed,domain=base
 for a,b in itertools.combinations(sorted(ps),2):
  for t in sorted(set(ps[a])&set(ps[b])):
   ya,pa,s=ps[a][t];yb,pb,_=ps[b][t];add(protocol+'__'+sys,protocol,sys,model,seed,domain,a+'−'+b,t,ya-yb,pa-pb,s)
for k,(y,hat,s) in direct.items():
 protocol,sys,model,seed,domain,pair,t=k;add(protocol+'__'+sys,protocol,sys,model,seed,domain,pair,t,y,hat,s)
CP=CROOT/'decisive_transfer_C_lopo_loso_20261020';cr,cs=frozen(CP/'point_predictions.csv');provenance.append(dict(tag='C_transfer',id=cs['id'],sha256=cs['sha256'],rows=len(cr)))
for r in cr:
 if r['region']=='outer_test' and r['pair_role']!='target_related_pair':skip['C_outer_source_pair_diagnostic_not_target_evaluation']+=1;continue
 hat=float(r['delta_G_pred_eV_atom'])
 if r['method']=='source_only_family_mean' and r['region']=='outer_test':
  assert (r['phase_a'],r['phase_b'])==('bcc','hcp');hat=-hat
 add(r['condition'],r['condition'],r['system'],r['method'],'deterministic',r['region'],r['phase_a']+'−'+r['phase_b'],float(r['T_K']),float(r['delta_G_ref_eV_atom']),hat,cs['id'])
seen=set()
for r in out:
 k=tuple(r[x] for x in ['case','model','seed','domain','system','pair','T_K']);assert k not in seen,k;seen.add(k)
groups=defaultdict(list)
for r in out:groups[(r['case'],r['model'],r['seed'],r['domain'])].append(r)
summary=[]
for (case,model,seed,domain),rs in groups.items():
 e=np.array([r['error_meV_atom'] for r in rs]);ab=abs(e)
 summary.append(dict(case=case,model=model,seed=seed,domain=domain,n_pair_temperature_points=len(e),MAE_meV_atom=float(ab.mean()),median_abs_meV_atom=float(np.median(ab)),P90_abs_meV_atom=float(np.quantile(ab,.9)),max_abs_meV_atom=float(ab.max()),bias_meV_atom=float(e.mean()),source_ids=';'.join(sorted({r['source_id'] for r in rs}))))
def savecsv(name,rs):
 with (R/name).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rs[0]));w.writeheader();w.writerows(rs)
savecsv('data/OPTIMIZATION_FIT_SUMMARY.csv',summary)
# Exact source design ranks: no new coefficient fitting.
ranks=[]
selection=json.loads((CP/'model_selection.json').read_text())
for fn in sorted((CP/'results').glob('*readout.npz')):
 digest=hashlib.sha256(fn.read_bytes()).hexdigest();assert any(m.get('readout_sha256')==digest for m in selection['models']);z=np.load(fn);X=z['source_pair_features'];ranks.append(dict(readout=fn.name,readout_sha256=digest,source_pair_vectors=X.shape[0],feature_dimension=X.shape[1],source_pair_rank=int(np.linalg.matrix_rank(X))))
labels={'cross_system_60_20_20_batch8_3000__sio2':'SiO₂ · matched temperature tail','cross_system_60_20_20_batch8_3000__ti':'Ti · matched temperature tail','cross_system_60_20_20_batch8_3000__zr':'Zr · matched temperature tail','hf_60_20_20_batch2_3000__hf':'Hf · batch 2 · matched tail','hf_60_20_20_batch8_10000_selected8000__hf':'Hf · selected 8000 · tail','sio2_lopo_cristobalite_beta':'SiO₂ · cristobalite LOPO','metals_loso_ti':'Ti · family LOSO (Hf/Zr sources)'}
assert set(labels)=={r['case'] for r in out}
p.SELECT.clear();plt=p.plt;options=[]
for i,(case,title) in enumerate(labels.items(),1):
 rs=[r for r in out if r['case']==case];training='fit' if case in ['sio2_lopo_cristobalite_beta','metals_loso_ti'] else 'train';evaluation='outer_test' if training=='fit' else 'test'
 curves=sorted({(r['model'],r['seed']) for r in rs if r['domain']==training});key=f'optimization_{i:02d}'
 nrows=int(np.ceil(len(curves)/3));bottom=.21+.037*nrows
 fig,axes=plt.subplots(1,2,figsize=(p.W,4.55+.15*nrows),sharey=True);fig.subplots_adjust(left=.10,right=.97,bottom=bottom,top=.85,wspace=.26)
 for ax,domain,letter,name in zip(axes,[training,evaluation],'ab',['Source fit' if training=='fit' else 'Training fit','Held-out target' if evaluation=='outer_test' else 'Temperature test']):
  for model,seed in curves:
   rr=groups[(case,model,seed,domain)];assert rr,(case,model,seed,domain)
   x=np.sort([r['absolute_error_meV_atom'] for r in rr]);y=np.arange(1,len(x)+1)/len(x)
   name_model=p.label(model)+(' · '+seed if seed not in ['none','deterministic'] else '')
   ls={'11':'-','23':'--','37':':','none':'-','deterministic':'-'}.get(seed,'-')
   if '3term' in model:ls='--'
   ax.step(x,y,where='post',label=name_model,color=p.color(model),ls=ls,lw=1.2,alpha=.9)
  ax.set_xscale('symlog',linthresh=.0001);ax.set(xlabel='Absolute ΔG error (meV/atom)',ylim=(0,1.025),yticks=[0,.5,1]);ax.set_xlim(left=0);p.clean(ax,False);p.panel(ax,letter,name)
 xmax=max(r['absolute_error_meV_atom'] for r in rs if r['domain'] in [training,evaluation]);
 for ax in axes:ax.set_xlim(0,xmax*1.08)
 axes[0].set_ylabel('Empirical fraction of paired points')
 fig.legend(*axes[0].get_legend_handles_labels(),loc='lower center',bbox_to_anchor=(.52,.012),ncol=3,frameon=False,fontsize=6.8,columnspacing=1.0,handlelength=1.6)
 plotted=[r for r in rs if r['domain'] in [training,evaluation]]
 p.export(key,title,fig,plotted,'All paired points in the source-fit and declared test domains. Curves keep seeds separate. The horizontal axis is linear below 0.0001 meV/atom and logarithmic above it. The C folds score only target-related outer pairs. Source and target supports are stated in the table.','optimization')
 options.append(dict(key=key,case=case,label=title,training_domain=training,evaluation_domain=evaluation,models=len(curves)))
meta=dict(options=options,summary=summary,source_design=ranks,provenance=provenance,exclusions=dict(skip),observations=len(out),score_revision='v2_phase_order_and_ties_20261027',policy='Family outer predictions align to bcc minus hcp after the independent audit. No pass threshold or new fitting. Read source-fit error before interpreting transfer. Errors at correlated temperatures are descriptive distributions. QH/zero/family controls have source-domain errors but no optimized training loss.')
(R/'data/optimization_check.json').write_text(json.dumps(meta,indent=2));(R/'data/optimization_check.js').write_text('window.FES_OPTIMIZATION='+json.dumps(meta)+';\n')
xp=R/'figures/explainer_figures.json';X=json.loads(xp.read_text());X['figures']=[f for f in X['figures'] if not f['key'].startswith('optimization_')]+p.FIGS
xp.write_text(json.dumps(X,ensure_ascii=False,indent=2));(R/'figures/explainer_figures.js').write_text('window.FES_EXPLAINER_FIGURES='+json.dumps(X,ensure_ascii=False)+';\n')
print('Audited',len(out),'points;',len(summary),'model/seed/domain distributions;',len(options),'cases; ranks',ranks)
