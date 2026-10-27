#!/usr/bin/env python3
"""Independent pair/envelope audit. Correct derived scores; keep frozen inputs."""
from pathlib import Path
from collections import defaultdict
import csv,json,hashlib,itertools
import numpy as np
R=Path(__file__).resolve().parents[1];W=R.parent.parent;D=json.loads((R/'data/review_data.json').read_text());P=W/'outputs/fes-bench-cross-system-20261020/results/decisive_transfer_C_lopo_loso_20261020'
s=next(s for s in D['sources'] if s['name']=='point_predictions.csv' and 'decisive_transfer_C' in s['path']);path=R/s['href'];assert hashlib.sha256(path.read_bytes()).hexdigest()==s['sha256']
with path.open() as f:rows=list(csv.DictReader(f))
TIE=1e-12
# The saved family control is explicitly hcp minus bcc. Verify it independently.
z=np.load(P/'results/metals_loso_ti__source_only.npz');family=np.asarray(z['delta_coeff']);saved=np.asarray(z['source_deltas']);assert np.allclose(family,saved.mean(0),atol=1e-15)
def basis(t):
 x=np.asarray(t)/1000.;return np.column_stack([np.ones_like(x),x,x*np.log(x),x*x])
original_family_check=[];corrected_points=[];fix_count=0
for r in rows:
 q=dict(r);y=float(r['delta_G_ref_eV_atom']);hat=float(r['delta_G_pred_eV_atom'])
 if r['method']=='source_only_family_mean':
  assert (r['phase_a'],r['phase_b']) in [('hcp','bcc'),('bcc','hcp')]
  declared=float((basis([float(r['T_K'])])@family)[0]);expected=declared if (r['phase_a'],r['phase_b'])==('hcp','bcc') else -declared
  if r['region']=='outer_test':original_family_check.append(abs(hat-declared));hat=expected;fix_count+=1
  else:assert abs(hat-expected)<1e-12
 q.update(delta_G_pred_eV_atom=hat,delta_G_ref_eV_atom=y);corrected_points.append(q)
assert max(original_family_check)<1e-12

def roots(t,y):
 out=[];i=0
 while i<len(t):
  if abs(y[i])<=TIE:
   j=i
   while j+1<len(t) and abs(y[j+1])<=TIE:j+=1
   out.append(float((t[i]+t[j])/2));i=j+1
  elif i+1<len(t) and y[i]*y[i+1]<0:
   out.append(float(t[i]+(t[i+1]-t[i])*abs(y[i])/(abs(y[i])+abs(y[i+1]))));i+=1
  else:i+=1
 return out

def labels(g):
 order=np.argsort(g,axis=1);valid=(g[np.arange(len(g)),order[:,1]]-g[np.arange(len(g)),order[:,0]])>TIE;return order[:,0],valid

def transitions(t,g,names):
 lab,_=labels(g);out=[]
 for i in range(len(t)-1):
  if lab[i]!=lab[i+1]:
   a,b=lab[i],lab[i+1];y0=g[i,a]-g[i,b];y1=g[i+1,a]-g[i+1,b]
   if y0*y1<0:tc=t[i]+(t[i+1]-t[i])*abs(y0)/(abs(y0)+abs(y1))
   else:tc=(t[i]+t[i+1])/2
   out.append(float(tc))
 return out

def event_stats(ref,pred):
 # Ordered one-to-one matching, maximize count then minimize total distance.
 states={(0,0):(0,0.,[])}
 def put(k,v):
  old=states.get(k)
  if old is None or (-v[0],v[1])<(-old[0],old[1]):states[k]=v
 for i in range(len(ref)+1):
  for j in range(len(pred)+1):
   if (i,j) not in states:continue
   n,c,ms=states[i,j]
   if i<len(ref):put((i+1,j),(n,c,ms))
   if j<len(pred):put((i,j+1),(n,c,ms))
   if i<len(ref) and j<len(pred) and abs(ref[i]-pred[j])<=100:put((i+1,j+1),(n+1,c+abs(ref[i]-pred[j]),ms+[(i,j)]))
 ms=states[len(ref),len(pred)][2]
 return dict(n_reference_events=len(ref),n_predicted_events=len(pred),n_matched_events=len(ms),missed_events=len(ref)-len(ms),false_events=len(pred)-len(ms),hit_conditioned_Tc_MAE_K=float(np.mean([abs(ref[i]-pred[j]) for i,j in ms])) if ms else '',reference_roots_K=';'.join(map(str,ref)),predicted_roots_K=';'.join(map(str,pred)),matched_Tc_abs_errors_K=';'.join(str(abs(ref[i]-pred[j])) for i,j in ms))

pairgroups=defaultdict(list);envgroups=defaultdict(list)
for r in corrected_points:
 pairgroups[tuple(r[k] for k in ['condition','method','region','system','phase_a','phase_b'])].append(r)
 envgroups[tuple(r[k] for k in ['condition','method','region','system'])].append(r)
oldpair={tuple(r[k] for k in ['condition','method','region','system','phase_a','phase_b']):r for r in D['transferPair']}
pair_metrics=[]
for key,rr in pairgroups.items():
 rr=sorted(rr,key=lambda r:float(r['T_K']));t=np.array([float(r['T_K']) for r in rr]);y=np.array([r['delta_G_ref_eV_atom'] for r in rr]);hat=np.array([r['delta_G_pred_eV_atom'] for r in rr]);valid=abs(y)>TIE;abstain=abs(hat)<=TIE;model=key[1]
 old=oldpair[key];q=dict(old);q.update(n_identifiable_reference_points=int(valid.sum()),reference_tie_points=int((~valid).sum()),predicted_abstentions_on_identifiable_points=int((abstain&valid).sum()),sign_accuracy=float(np.mean((np.sign(y[valid])==np.sign(hat[valid]))&~abstain[valid])) if model!='zero_delta_G' else '',delta_G_MAE_eV_atom=float(abs(y-hat).mean()),delta_G_RMSE_eV_atom=float(np.sqrt(np.mean((y-hat)**2))),score_revision='v2_phase_order_and_ties_20261027',original_sign_accuracy=old['sign_accuracy'],original_delta_G_MAE_eV_atom=old['delta_G_MAE_eV_atom'])
 # Keep the declared exact-zero root convention for event continuity with the archive.
 q.update(event_stats(roots(t,y),[] if model=='zero_delta_G' else roots(t,hat)));pair_metrics.append(q)

envelope=[];participation=[];phase_points=[];consistency=[];model_shapes=[]
oldenv={tuple(r[k] for k in ['condition','method','region','system']):r for r in D['transferEnvelope']}
for key,rr in envgroups.items():
 condition,model,region,system=key
 if model=='zero_delta_G':continue
 names=sorted({r['phase_a'] for r in rr}|{r['phase_b'] for r in rr});anchor=names[0];maps={}
 for r in rr:maps.setdefault((r['phase_a'],r['phase_b']),{})[float(r['T_K'])]=(r['delta_G_ref_eV_atom'],r['delta_G_pred_eV_atom'])
 temp=sorted(set.intersection(*[set(m) for m in maps.values()]));assert all(set(temp)==set(m) for m in maps.values());t=np.array(temp)
 def gap(a,b,k):
  if a==b:return np.zeros(len(t))
  if (a,b) in maps:return np.array([maps[a,b][tt][k] for tt in t])
  return -np.array([maps[b,a][tt][k] for tt in t])
 gr=np.column_stack([gap(n,anchor,0) for n in names]);gp=np.column_stack([gap(n,anchor,1) for n in names]);cycle=0.
 for a,b in maps:
  cycle=max(cycle,float(np.max(abs(gp[:,names.index(a)]-gp[:,names.index(b)]-gap(a,b,1)))))
 assert cycle<1e-10,(key,cycle)
 lr,valid=labels(gr);lp,pred_valid=labels(gp);accuracy=float(np.mean((lr[valid]==lp[valid])&pred_valid[valid]));refroots=transitions(t,gr,names);predroots=transitions(t,gp,names);events=event_stats(refroots,predroots)
 old=oldenv[key];q=dict(old);q.update(stable_envelope_accuracy=accuracy,n_identifiable_reference_points=int(valid.sum()),reference_tie_points=int((~valid).sum()),predicted_abstentions_on_identifiable_points=int((~pred_valid&valid).sum()),reference_transition_count=events['n_reference_events'],predicted_transition_count=events['n_predicted_events'],matched_transition_count=events['n_matched_events'],missed_transition_count=events['missed_events'],false_transition_count=events['false_events'],hit_conditioned_transition_Tc_MAE_K=events['hit_conditioned_Tc_MAE_K'],reference_transitions_K=events['reference_roots_K'],predicted_transitions_K=events['predicted_roots_K'],score_revision='v2_phase_order_and_ties_20261027',original_stable_envelope_accuracy=old['stable_envelope_accuracy']);envelope.append(q)
 if len(names)==2:
  pair=next(p for p in pair_metrics if tuple(p[k] for k in ['condition','method','region','system'])==key);assert abs(accuracy-pair['sign_accuracy'])<1e-12
  consistency.append(dict(condition=condition,method=model,region=region,system=system,pair_sign_accuracy=pair['sign_accuracy'],stable_envelope_accuracy=accuracy,identifiable_points=int(valid.sum()),reference_ties=int((~valid).sum()),original_pair_sign_accuracy=pair['original_sign_accuracy'],original_envelope_accuracy=old['stable_envelope_accuracy']))
 if condition=='sio2_lopo_cristobalite_beta' and region=='outer_test':
  source=[names.index(n) for n in names if n!='cristobalite_beta'];source_labels=np.array(source)[np.argmin(gp[:,source],axis=1)];target=names.index('cristobalite_beta');tr=(lr==target)&valid;tp=(lp==target)&pred_valid
  participation.append(dict(method=model,n_states=len(t),reference_target_stable_states=int(tr.sum()),predicted_target_stable_states=int(tp.sum()),decisions_changed_by_target=int((source_labels!=lp).sum()),source_only_accuracy_against_full_reference=float(np.mean(source_labels[valid]==lr[valid])),full_envelope_accuracy=accuracy,target_stable_recall=float(np.mean(tp[tr])),target_phase_precision=float(np.mean(tr[tp])) if tp.any() else '',false_target_stable_states=int((tp&~tr).sum()),missed_target_stable_states=int((tr&~tp).sum())))
  for i,tt in enumerate(t):
   phase_points.append(dict(method=model,T_K=float(tt),reference_quartz_eV_atom=float(gr[i,names.index('quartz_beta')]-gr[i,names.index('quartz_beta')]),reference_tridymite_minus_quartz_eV_atom=float(gr[i,names.index('tridymite_p63mmc')]-gr[i,names.index('quartz_beta')]),reference_cristobalite_minus_quartz_eV_atom=float(gr[i,target]-gr[i,names.index('quartz_beta')]),predicted_tridymite_minus_quartz_eV_atom=float(gp[i,names.index('tridymite_p63mmc')]-gp[i,names.index('quartz_beta')]),predicted_cristobalite_minus_quartz_eV_atom=float(gp[i,target]-gp[i,names.index('quartz_beta')]),reference_stable_phase=names[lr[i]],predicted_source_only_stable_phase=names[source_labels[i]],predicted_full_stable_phase=names[lp[i]]))
  if model in ['C-simple','C-DPA']:
   gt=gp[:,names.index('tridymite_p63mmc')]-gp[:,names.index('quartz_beta')];gc=gp[:,target]-gp[:,names.index('quartz_beta')];factor=float(np.dot(gc,gt)/np.dot(gt,gt));res=float(np.max(abs(gc-factor*gt))*1000)
   w=np.load(P/'results'/f'{condition}__{model.replace("-","_")}_readout.npz')['W'];sv=np.linalg.svd(w,compute_uv=False);rr=[p for p in pair_metrics if p['condition']==condition and p['method']==model and p['region']=='outer_test'];pr=[float(v) for p in rr for v in p['predicted_roots_K'].split(';') if v]
   model_shapes.append(dict(method=model,cristobalite_minus_quartz_projection_on_source_gap=factor,max_projection_residual_meV_atom=res,readout_singular_values=sv.tolist(),secondary_to_primary_singular_ratio=float(sv[1]/sv[0]),predicted_pair_root_spread_K=max(pr)-min(pr),predicted_pair_roots_K=pr))

skill=[]
for condition in ['sio2_lopo_cristobalite_beta','metals_loso_ti']:
 target=[r for r in corrected_points if r['condition']==condition and r['region']=='outer_test' and r['pair_role']=='target_related_pair'];zero=[r for r in target if r['method']=='zero_delta_G'];zmae=float(np.mean([abs(r['delta_G_ref_eV_atom'])*1000 for r in zero]));zkeys={(r['phase_a'],r['phase_b'],r['T_K']) for r in zero}
 for model in sorted({r['method'] for r in target}):
  vals=[r for r in target if r['method']==model];assert {(r['phase_a'],r['phase_b'],r['T_K']) for r in vals}==zkeys;mae=float(np.mean([abs(r['delta_G_ref_eV_atom']-r['delta_G_pred_eV_atom'])*1000 for r in vals]));env=next((r for r in envelope if r['condition']==condition and r['region']=='outer_test' and r['method']==model),None)
  skill.append(dict(condition=condition,method=model,n_pair_temperature_points=len(vals),MAE_meV_atom=mae,zero_gap_MAE_meV_atom=zmae,relative_skill=1-mae/zmae if zmae>0 else '',stable_envelope_accuracy=env['stable_envelope_accuracy'] if env else '',decision_status='reference ties excluded' if env else 'NA: fully degenerate predictor'))

def writecsv(name,rr):
 with (R/name).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rr[0]));w.writeheader();w.writerows(rr)
writecsv('data/TRANSFER_PAIR_METRICS_V2.csv',pair_metrics);writecsv('data/TRANSFER_ENVELOPE_METRICS_V2.csv',envelope);writecsv('data/TI_PAIR_ENVELOPE_CONSISTENCY.csv',consistency);writecsv('data/SIO2_TARGET_PARTICIPATION.csv',participation);writecsv('data/TRANSFER_RELATIVE_SKILL.csv',skill);writecsv('figures/source_data/sio2_envelope_audit.csv',phase_points)
audit=dict(revision='v2_phase_order_and_ties_20261027',pair_metrics=pair_metrics,envelope_metrics=envelope,participation=participation,sio2_shapes=model_shapes,relative_skill=skill,two_phase_checks=consistency,family_outer_rows_corrected=fix_count,tie_tolerance_eV_atom=TIE,original_point_source_id=s['id'],original_point_source_sha256=s['sha256'],family_readout_sha256=hashlib.sha256((P/'results/metals_loso_ti__source_only.npz').read_bytes()).hexdigest(),policy='Family coefficients are hcp minus bcc; align every pair to its labels. Exclude reference ties from sign/envelope denominators; count predicted ties as abstentions. Fully degenerate zero-gap decisions are NA. Original frozen files remain unchanged.')
(R/'data/transfer_score_audit.json').write_text(json.dumps(audit,indent=2));(R/'data/transfer_score_audit.js').write_text('window.FES_TRANSFER_AUDIT='+json.dumps(audit)+';\n')
print(json.dumps(dict(participation=participation,sio2_shapes=model_shapes,skill=skill,family_outer_rows_corrected=fix_count,two_phase_checks=len(consistency)),indent=2))
