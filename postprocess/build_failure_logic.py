#!/usr/bin/env python3
"""Build a failure-first view without changing the frozen evidence snapshot.

OBS means a reported result, not a pass/fail threshold. NR means no matched result
in this view, not proof that a method has never run. NA is an input restriction.
RF is a missing compatible reference grid. NE is no reference event.
"""
from pathlib import Path
from collections import defaultdict
import json,csv,hashlib
R=Path(__file__).resolve().parents[1];D=json.loads((R/'data/review_data.json').read_text());S={s['id']:s for s in D['sources']}
# Current score layer; the original snapshot and source hashes remain unchanged.
transfer_audit=R/'data/transfer_score_audit.json'
if transfer_audit.exists():
 current=json.loads(transfer_audit.read_text());D['transferPair']=current['pair_metrics'];D['transferEnvelope']=current['envelope_metrics']
def source_rows(ident):
 with (R/S[ident]['href']).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
methods=[('thermo4','Direct\nThermo4'),('residual4','QH +\nThermo4'),('pid','Phase-ID'),('t3','T3'),('t2','T2'),('simple','C-simple'),('dpa','C-DPA'),('qh','QH-like'),('family','Family\nmean')]
map_method={'thermo_phase_id_4term':'thermo4','same_form_4term_analytic':'thermo4','qh_residual_thermo_phase_id_4term':'residual4','temperature_phase_id_mlp_pairloss':'pid','T3_direct_G':'t3','T2_qh_residual':'t2','QH_only':'qh','C-simple':'simple','C-DPA':'dpa','QH_like_numerical':'qh','source_only_family_mean':'family'}
cases=[]
for sys in ['ti','zr','sio2','hf']:
 cases.append(dict(id=sys+'_tail',system=sys,task='Temperature OOD',label=({'sio2':'SiO₂','hf':'Hf','ti':'Ti','zr':'Zr'}[sys])+' · temperature tail'+(' · b2' if sys=='hf' else ' · b8'),budget='Target training-temperature labels allowed',scope='Frozen 60/20/20 tail; test has no reference event.',reference='scorable',evidence=['E001','E004','E009'],protocol='hf_60_20_20_batch2_3000' if sys=='hf' else 'cross_system_60_20_20_batch8_3000'))
cases += [dict(id='sio2_lopo',system='sio2',task='Phase OOD',label='SiO₂ · phase OOD',budget='No target-phase free-energy labels',scope='One cristobalite outer fold; source-only selection.',reference='scorable',evidence=['E002','E008','E009'],condition='sio2_lopo_cristobalite_beta'),dict(id='ti_loso',system='ti',task='System OOD',label='Ti · system OOD',budget='No target-system free-energy labels',scope='One family outer fold; source-only selection.',reference='scorable',evidence=['E002','E008','E009'],condition='metals_loso_ti'),dict(id='c_pressure',system='C',task='Pressure transfer',label='C · pressure transfer',budget='Transfer split and label budget not frozen',scope='Paired QHA references exist. No matched pressure-transfer model score is recorded.',reference='scorable_not_evaluated',evidence=['E011']),dict(id='ca_new',system='CaSiO3',task='New paired grid',label='CaSiO₃ · new paired grid',budget='Compatible target grid not closed',scope='Single-volume candidates and MD do not close the new two-phase reference.',reference='blocked',evidence=['E017']),dict(id='mg_new',system='MgSiO3',task='New paired grid',label='MgSiO₃ · new paired grid',budget='Compatible target grid not closed',scope='The new Pv/PPv paired grid remains open.',reference='blocked',evidence=['E016']),dict(id='fe_new',system='Fe',task='New three-phase grid',label='Fe · new three-phase grid',budget='Compatible target grid not closed',scope='This row is the new full grid. It does not reject the twelve admitted source fcc/hcp pairs.',reference='blocked',evidence=['E018'])]
cells={(case['id'],mid):dict(case=case['id'],method=mid,state='RF' if case['reference']=='blocked' else 'NR',records=[],sources=[],note=case['scope']) for case in cases for mid,_ in methods}
def add(case,method,record,source):
 cell=cells[(case,method)];cell['state']='OBS';cell['records'].append(record);cell['sources']=sorted(set(cell['sources']+[source]));cell['note']='Observed score. The row retains its own task and information budget.'
for row in D['tail']:
 if row['region']!='test' or row['model'] not in map_method:continue
 sys=row['system'];case=next((c for c in cases if c.get('protocol')==row['protocol_id'] and c['system']==sys),None)
 if case:add(case['id'],map_method[row['model']],dict(seed=row['seed'],accuracy=float(row['stable_phase_accuracy']),mae=float(row['delta_G_MAE_meV_atom']),reference_events=int(row['reference_pair_roots'] or 0),missed_events=None,false_events=int(row['false_pair_roots']) if row['false_pair_roots'] else None), 'S00003')
# The later SiO2 residual run is independently scored and explicitly matched in S00001.
# It was absent from the earlier unified table; keep its direct source separate.
for ident,mid in [('S00611','t2'),('S00808','qh')]:
 for row in source_rows(ident):
  if row['region']=='test':add('sio2_tail',mid,dict(seed=row['seed'],accuracy=float(row['stable_envelope_accuracy']),mae=float(row['mean_pair_delta_G_MAE_eV_atom'])*1000,reference_events=0,missed_events=None,false_events=0),ident)
for row in source_rows('S00585'):
 if row['region']=='test' and row['method']=='qh_residual_thermo_phase_id_4term':add('sio2_tail','residual4',dict(seed='deterministic',accuracy=float(row['stable_envelope_accuracy']),mae=float(row['mean_pair_delta_G_MAE_eV_atom'])*1000,reference_events=0,missed_events=None,false_events=0),'S00585')
# This analytic Hf control is reported in the same matched Hf diagnostic.
add('hf_tail','residual4',dict(seed='deterministic',accuracy=1.0,mae=.234,reference_events=0,missed_events=None,false_events=0),'S00013')
for case in [c for c in cases if 'condition' in c]:
 for mid in ['thermo4','residual4','pid']:
  cells[(case['id'],mid)].update(state='NA',note='This named control needs labels or a phase ID for the held-out target. Those inputs are not allowed in this zero-target-label fold.')
 for row in D['transferEnvelope']:
  if row['condition']!=case['condition'] or row['region']!='outer_test' or row['method'] not in map_method:continue
  pairs=[p for p in D['transferPair'] if p['condition']==row['condition'] and p['region']=='outer_test' and p['method']==row['method'] and p['pair_role']=='target_related_pair']
  mae=sum(float(p['delta_G_MAE_eV_atom'])*1000 for p in pairs)/len(pairs) if pairs else None
  add(case['id'],map_method[row['method']],dict(seed='one frozen readout',accuracy=float(row['stable_envelope_accuracy']),mae=mae,reference_events=int(row['reference_transition_count']),missed_events=int(row['missed_transition_count']),false_events=int(row['false_transition_count'])),'S00006')
  cells[(case['id'],map_method[row['method']])]['sources']=sorted(set(cells[(case['id'],map_method[row['method']])]['sources']+['S00005']))
# Hf test-event counts are explicit in the matched diagnostic (S00013).
for (case,mid),cell in cells.items():
 if case=='hf_tail' and cell['state']=='OBS':
  cell['sources']=sorted(set(cell['sources']+['S00013']))
  for record in cell['records']:
   record['reference_events']=0
   record['false_events']=1 if mid=='t3' else 0
# Exact summaries; ranges retain seed extrema and are not uncertainty intervals.
for cell in cells.values():
 rs=cell['records']
 for metric in ['accuracy','mae','reference_events','missed_events','false_events']:
  vals=[r[metric] for r in rs if r.get(metric) is not None]
  cell[metric+'_min']=min(vals) if vals else None;cell[metric+'_max']=max(vals) if vals else None
 cell['n_records']=len(rs)
questions=[dict(id='where',title='Where do methods fail?',text='Compare tasks under a stated input budget. Record measured phase errors, event misses, applicability limits and missing evaluations.',result='I · Failure landscape',anchor='failure-landscape'),dict(id='what',title='What relative information is lost?',text='Inspect the unanchored phase gap, relative entropy slope and the margin between the lowest competing phases.',result='II · Thermodynamic information loss',anchor='information-loss'),dict(id='why',title='Why, and can physical priors help?',text='Use matched controls to separate missing inputs, fit quality, reference mismatch and baseline validity.',result='III · Physical baseline contribution',anchor='baseline-contribution')]
rules=[dict(observation='Single-phase shape is accurate, but ΔG is wrong.',candidate='Phase-dependent offsets or reference alignment.',test='Replace the static gap; evaluate unanchored ΔG.',state='Diagnostic hypothesis; not a cause assignment.',evidence=['E003','E010']),dict(observation='Train fit is good, but the temperature tail drifts.',candidate='Relative entropy or higher-order temperature response.',test='Compare slope and a matched temperature-function control.',state='Check train fit before calling this transfer failure.',evidence=['E001','E005','E007']),dict(observation='Mean ΔG MAE is low, but the stable phase is wrong.',candidate='A small margin between the lowest competing phases.',test='Inspect the responsible lower-envelope pair and error direction.',state='SiO₂ has a counterexample; its pair mechanism remains open.',evidence=['E001']),dict(observation='Performance drops when a phase is held out.',candidate='Structure-to-thermodynamics mapping for an unseen phase.',test='Compare C-simple, C-DPA and a matched few-shot extension.',state='Two current outer folds; no general cause established.',evidence=['E002']),dict(observation='QH residual does not improve a matched control.',candidate='The baseline may add no usable information for that budget.',test='Compare the same temperature prior and verify the QHA definition.',state='No consistent known-phase gain; causal transfer test open.',evidence=['E001','E012','E015']),dict(observation='Performance drops across systems.',candidate='Energy references, representation or thermal-response differences.',test='Use family folds, static-gap replacement and matched support.',state='Do not infer the cause from one Ti fold.',evidence=['E002','E014'])]
mechanisms=[dict(system='SiO₂',focus='Wrong stable phase despite lower mean pair error; unseen-phase prediction.',information='Critical low-phase gap and the relative location of an unseen phase.',state='Observed decision counterexample; one LOPO fold.',evidence=['E001','E002']),dict(system='Ti / Zr / Hf',focus='Learned fit gaps and family transfer.',information='Relative slope and cross-element phase gaps.',state='Tail diagnostics available; family transfer still limited.',evidence=['E001','E004','E005','E007']),dict(system='CaSiO₃',focus='High-pressure reference and soft-mode validity.',information='Pressure response and baseline applicability.',state='New compatible paired reference remains open.',evidence=['E017']),dict(system='MgSiO₃',focus='PPv soft mode and a compatible paired reference.',information='Anharmonic response and pressure-state reference validity.',state='New paired route remains open.',evidence=['E016']),dict(system='Fe',focus='Sparse comparable states and limited event coverage.',information='Electronic-temperature convention and relative stability.',state='Source pairs admitted; new full grid remains open.',evidence=['E018']),dict(system='C',focus='Pressure and bonding changes within a scoped reference.',information='Volume and pressure response under the QHA contract.',state='Paired QHA reference exists; model pressure-transfer test open.',evidence=['E011'])]
controls=[dict(id='M₀',inputs='Structure',model='Predict G directly.',purpose='Reference for the learned structural pipeline.'),dict(id='M₁',inputs='Structure + simple thermal prior',model='Use a matched-complexity temperature prior.',purpose='Control for ordinary smooth temperature information.'),dict(id='M₂',inputs='Structure + verified QHA',model='Use a physically valid baseline and a matched residual.',purpose='Test the added value of explicit vibrational information.')]
payload=dict(questions=questions,methods=[dict(id=i,label=l.replace('\n',' '),plot_label=l) for i,l in methods],cases=cases,cells=list(cells.values()),rules=rules,systems=mechanisms,controls=controls,states={'OBS':'Observed score','NR':'No matched reported evaluation','NA':'Inputs not allowed for this task','RF':'Compatible new reference grid missing','NE':'No reference event'},snapshot_sha256=hashlib.sha256((R/'data/review_data.json').read_bytes()).hexdigest(),notes=['A measured accuracy below 1 means at least one phase decision was wrong. No new MAE threshold is introduced.','Accuracy 1 in a tail does not establish boundary recovery.','Ranges show seed extrema, not confidence intervals.','The pressure row has QHA references but no frozen matched model test.','Reference-blocked rows refer only to the named new grids, not every source dataset for that material.'])
(R/'data/failure_logic.json').write_text(json.dumps(payload,ensure_ascii=False,indent=2));(R/'data/failure_logic.js').write_text('window.FES_FAILURE='+json.dumps(payload,ensure_ascii=False)+';\n')
with (R/'data/FAILURE_LANDSCAPE.csv').open('w',newline='') as f:
 rows=[{**{k:v for k,v in cell.items() if k!='records'},'sources':';'.join(cell['sources'])} for cell in payload['cells']];w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print('Failure landscape:',len(cases),'cases ×',len(methods),'methods;',sum(c['state']=='OBS' for c in payload['cells']),'observed cells')
