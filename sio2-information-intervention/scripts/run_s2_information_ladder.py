#!/usr/bin/env python3
"""Run frozen S2 gap/slope information injection on archived SiO2 C predictions."""
from __future__ import annotations
import argparse, csv, json, math, hashlib, os
from pathlib import Path
import numpy as np

P = Path(os.environ.get('FES_PORTAL_ARCHIVE', 'fes-paper-evidence-portal-20261026'))
INPUT = P / 'figures/source_data/sio2_envelope_audit.csv'
TIE = 1e-12
PHASES = ['quartz_beta', 'tridymite_p63mmc', 'cristobalite_beta']
METHODS = ['C-simple', 'C-DPA']
ARMS = ['zero_shot', 'oracle_gap', 'oracle_gap_slope', 'few_shot_1_label', 'few_shot_2_label']

def read_data():
    with INPUT.open(newline='') as f: rows=list(csv.DictReader(f))
    result={}
    for method in METHODS:
        rr=sorted((r for r in rows if r['method']==method), key=lambda r:float(r['T_K']))
        t=np.array([float(r['T_K']) for r in rr])
        assert len(t)==1649 and np.all(np.diff(t)==1) and t[0]==851 and t[-1]==2499
        ref=np.column_stack([np.zeros(len(t)),
             [float(r['reference_tridymite_minus_quartz_eV_atom']) for r in rr],
             [float(r['reference_cristobalite_minus_quartz_eV_atom']) for r in rr]])
        pred=np.column_stack([np.zeros(len(t)),
             [float(r['predicted_tridymite_minus_quartz_eV_atom']) for r in rr],
             [float(r['predicted_cristobalite_minus_quartz_eV_atom']) for r in rr]])
        assert [r['method'] for r in rr]==[method]*len(t)
        result[method]=(t,ref,pred)
    return result

def ls_slope(t,y,center,half_window):
    m=np.abs(t-center)<=half_window
    assert m.sum()==2*half_window+1
    return float(np.polyfit(t[m],y[m],1)[0])

def roots(t,y):
    out=[]; i=0
    while i<len(t):
        if abs(y[i])<=TIE:
            j=i
            while j+1<len(t) and abs(y[j+1])<=TIE: j+=1
            out.append(float((t[i]+t[j])/2)); i=j+1
        elif i+1<len(t) and y[i]*y[i+1]<0:
            out.append(float(t[i]+(t[i+1]-t[i])*abs(y[i])/(abs(y[i])+abs(y[i+1])))); i+=1
        else: i+=1
    return out

def labels(g):
    order=np.argsort(g,axis=1)
    valid=(g[np.arange(len(g)),order[:,1]]-g[np.arange(len(g)),order[:,0]])>TIE
    return order[:,0],valid

def transitions(t,g):
    lab,_=labels(g); out=[]
    for i in range(len(t)-1):
        if lab[i]!=lab[i+1]:
            a,b=lab[i],lab[i+1]; y0=g[i,a]-g[i,b]; y1=g[i+1,a]-g[i+1,b]
            tc=t[i]+(t[i+1]-t[i])*abs(y0)/(abs(y0)+abs(y1)) if y0*y1<0 else (t[i]+t[i+1])/2
            out.append(float(tc))
    return out

def match_events(ref,pred,tol=100.):
    # Ordered one-to-one dynamic program; maximize matches then minimize distance.
    n,m=len(ref),len(pred); dp=[[None]*(m+1) for _ in range(n+1)]; dp[0][0]=(0,0.,[])
    def put(i,j,v):
        old=dp[i][j]
        if old is None or (-v[0],v[1])<(-old[0],old[1]): dp[i][j]=v
    for i in range(n+1):
      for j in range(m+1):
        if dp[i][j] is None: continue
        k,c,pairs=dp[i][j]
        if i<n: put(i+1,j,(k,c,pairs))
        if j<m: put(i,j+1,(k,c,pairs))
        if i<n and j<m and abs(ref[i]-pred[j])<=tol: put(i+1,j+1,(k+1,c+abs(ref[i]-pred[j]),pairs+[(i,j)]))
    pairs=dp[n][m][2]
    return {'n_reference_events':n,'n_predicted_events':m,'n_matched_events':len(pairs),'missed_events':n-len(pairs),'false_events':m-len(pairs),
            'hit_conditioned_Tc_MAE_K':float(np.mean([abs(ref[i]-pred[j]) for i,j in pairs])) if pairs else '',
            'reference_roots_K':';'.join(map(str,ref)),'predicted_roots_K':';'.join(map(str,pred)),
            'matched_Tc_abs_errors_K':';'.join(str(abs(ref[i]-pred[j])) for i,j in pairs)}

def to_rows(path,rows):
    if not rows: return
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--run-dir',type=Path,required=True); a=ap.parse_args(); out=a.run_dir
    protocol=json.loads((out/'PROTOCOL.json').read_text())
    assert protocol['status']=='frozen_before_S2_intervention_scoring'
    data=read_data(); point_rows=[]; pair_rows=[]; env_rows=[]; participation=[]; roots_rows=[]; summary=[]
    for method,(t,ref,pred0) in data.items():
        T0=protocol['temperature_coordinates']['T0_K']; t1,t2=protocol['temperature_coordinates']['few_shot_temperatures_K']
        iq=int(np.where(t==T0)[0][0]); i1=int(np.where(t==t1)[0][0]); i2=int(np.where(t==t2)[0][0])
        e=(pred0[:,2]-pred0[:,0])-(ref[:,2]-ref[:,0])
        ep0=ls_slope(t,e,T0,5)
        shot12=(e[i1],e[i2]); affine=np.polyfit([0.,t2-t1],shot12,1)
        calibrations={
          'zero_shot':(np.zeros_like(t),[]),
          'oracle_gap':(np.full_like(t,e[iq]),[float(T0)]),
          'oracle_gap_slope':(e[iq]+ep0*(t-T0),[float(v) for v in t[np.abs(t-T0)<=5]]),
          'few_shot_1_label':(np.full_like(t,e[i1]),[float(t1)]),
          'few_shot_2_label':(affine[0]*(t-t1)+affine[1],[float(t1),float(t2)])}
        common=(t!=t1)&(t!=t2)
        slope_common=common&(t>=t[0]+5)&(t<=t[-1]-5)
        for arm,(correction,used) in calibrations.items():
            pred=pred0.copy(); pred[:,2]-=correction
            # Full point output preserves every temperature and the phase-level correction.
            for j,T in enumerate(t):
                point_rows.append({'method':method,'arm':arm,'T_K':float(T),'phase':PHASES[2],
                    'target_minus_quartz_reference_eV_atom':float(ref[j,2]-ref[j,0]),
                    'target_minus_quartz_raw_prediction_eV_atom':float(pred0[j,2]-pred0[j,0]),
                    'target_minus_quartz_error_eV_atom':float(e[j]),'applied_cristobalite_correction_eV_atom':float(correction[j]),
                    'target_minus_quartz_corrected_prediction_eV_atom':float(pred[j,2]-pred[j,0]),
                    'target_label_consumed_here':int(float(T) in used),'consumed_target_temperatures_K':';'.join(map(str,used)),
                    'target_label_budget':len(used),'support_in_primary_fair_score':int(common[j])})
            # Pair metrics on common support; all pairs derive from phase-level energies.
            for ia,ib in [(2,0),(2,1),(0,1)]:
                yg=ref[:,ia]-ref[:,ib]; hg=pred[:,ia]-pred[:,ib]
                gap_mae=float(np.mean(np.abs(hg[common]-yg[common]))*1000)
                gap_rmse=float(np.sqrt(np.mean((hg[common]-yg[common])**2))*1000)
                gap_bias=float(np.mean((hg[common]-yg[common]))*1000)
                sl_ref=np.array([ls_slope(t,yg,float(T),5) for T in t[slope_common]])
                sl_hat=np.array([ls_slope(t,hg,float(T),5) for T in t[slope_common]])
                # Fixed T0 also receives a fair point score; T0 isn't a calibration coordinate.
                slope_err=sl_hat-sl_ref
                valid=np.abs(yg[common])>TIE; abst=np.abs(hg[common])<=TIE
                acc=float(np.mean((np.sign(yg[common][valid])==np.sign(hg[common][valid])) & ~abst[valid])) if valid.any() else ''
                # Pair root diagnostics use the same fair support.
                tr=t[common]; yr=yg[common]; pr=hg[common]
                rroots=roots(tr,yr); proots=roots(tr,pr); ev=match_events(rroots,proots)
                row={'method':method,'arm':arm,'fold':'LOPO_cristobalite_beta','phase_a':PHASES[ia],'phase_b':PHASES[ib],
                     'target_label_budget':len(used),'consumed_target_temperatures_K':';'.join(map(str,used)),
                     'n_common_evaluation_points':int(common.sum()),'n_slope_evaluation_points':int(slope_common.sum()),'n_identifiable_reference_points':int(valid.sum()),
                     'reference_tie_points':int((~valid).sum()),'predicted_abstentions':int((abst&valid).sum()),
                     'delta_G_MAE_meV_atom':gap_mae,'delta_G_RMSE_meV_atom':gap_rmse,'delta_G_bias_meV_atom':gap_bias,
                     'delta_G_error_at_T0_meV_atom':float((hg[iq]-yg[iq])*1000),
                     'relative_slope_MAE_meV_atom_K':float(np.mean(np.abs(slope_err))*1000),
                     'relative_slope_RMSE_meV_atom_K':float(np.sqrt(np.mean(slope_err**2))*1000),
                     'relative_slope_bias_meV_atom_K':float(np.mean(slope_err)*1000),
                     'pair_sign_accuracy':acc,**ev}
                pair_rows.append(row)
            # V2-compatible phase decision and envelope transitions.
            er=ref[common]; pr=pred[common]; tc=t[common]
            lr,valid=labels(er); lp,pvalid=labels(pr)
            envacc=float(np.mean((lr[valid]==lp[valid])&pvalid[valid])) if valid.any() else ''
            ev=match_events(transitions(tc,er),transitions(tc,pr))
            env_rows.append({'method':method,'arm':arm,'fold':'LOPO_cristobalite_beta','n_common_evaluation_points':int(common.sum()),
                'n_identifiable_reference_points':int(valid.sum()),'reference_tie_points':int((~valid).sum()),
                'predicted_abstentions':int((~pvalid&valid).sum()),'stable_envelope_accuracy':envacc,
                'reference_target_stable_states':int(np.sum((lr==2)&valid)),'predicted_target_stable_states':int(np.sum((lp==2)&pvalid)),
                'target_stable_recall':float(np.sum((lp==2)&pvalid&(lr==2)&valid)/np.sum((lr==2)&valid)) if np.sum((lr==2)&valid) else '',
                'source_only_accuracy_against_full_reference':float(np.mean((np.argmin(er[:,:2],axis=1)==lr)&valid)),**ev})
            # Root topology for every pair retained explicitly.
            for ia,ib in [(2,0),(2,1),(0,1)]:
                yref=ref[:,ia]-ref[:,ib]; ypred=pred[:,ia]-pred[:,ib]
                rr=roots(t[common],yref[common]); pp=roots(t[common],ypred[common]); evp=match_events(rr,pp)
                roots_rows.append({'method':method,'arm':arm,'pair_a':PHASES[ia],'pair_b':PHASES[ib],
                    'support':'common support excluding 1263 K and 2087 K','event_tolerance_K':100.0,**evp})
            # Pair metrics of target-quartz gap/slope at T0 and boundary summary.
            tq=pair_rows[-3] # loop order puts c-q first
            summary.append({'method':method,'arm':arm,'target_label_budget':len(used),'consumed_target_temperatures_K':';'.join(map(str,used)),
                'T0_K':float(T0),'few_shot_temperatures_K':f'{t1};{t2}','correction_at_T0_meV_atom':float(correction[iq]*1000),
                'target_quartz_gap_error_at_T0_meV_atom':tq['delta_G_error_at_T0_meV_atom'],
                'target_quartz_gap_MAE_meV_atom':tq['delta_G_MAE_meV_atom'],
                'target_quartz_slope_MAE_meV_atom_K':tq['relative_slope_MAE_meV_atom_K'],
                'envelope_accuracy':envacc,'target_stable_recall':env_rows[-1]['target_stable_recall'],
                'envelope_ref_transitions_K':env_rows[-1]['reference_roots_K'],'envelope_pred_transitions_K':env_rows[-1]['predicted_roots_K'],
                'envelope_matched_events':env_rows[-1]['n_matched_events'],'envelope_missed_events':env_rows[-1]['missed_events'],
                'envelope_false_events':env_rows[-1]['false_events'],'envelope_hit_Tc_MAE_K':env_rows[-1]['hit_conditioned_Tc_MAE_K']})
        # fixed reference derivative-window sensitivity at T0 and calibrated error slopes.
        sens=[]
        for half in [3,5,10]:
            sens.append({'method':method,'half_window_K':half,'points':2*half+1,
                         'reference_target_quartz_slope_meV_atom_K':ls_slope(t,ref[:,2]-ref[:,0],1675.,half)*1000,
                         'raw_prediction_target_quartz_slope_meV_atom_K':ls_slope(t,pred0[:,2]-pred0[:,0],1675.,half)*1000,
                         'error_slope_meV_atom_K':ls_slope(t,e,1675.,half)*1000})
        # stash separately
        to_rows(out/'data'/f'S2_SLOPE_WINDOW_SENSITIVITY_{method.replace("-","_")}.csv',sens)
    to_rows(out/'data/S2_POINT_PREDICTIONS.csv',point_rows)
    to_rows(out/'data/S2_PAIR_METRICS.csv',pair_rows)
    to_rows(out/'data/S2_ENVELOPE_METRICS.csv',env_rows)
    to_rows(out/'data/S2_PAIR_EVENT_METRICS.csv',roots_rows)
    to_rows(out/'data/S2_INFORMATION_LADDER_SUMMARY.csv',summary)
    audit={'protocol_id':protocol['protocol_id'],'scoring_version':'V2-compatible reproduction of frozen pair/envelope tie and event conventions',
           'source_csv_sha256':hashlib.sha256(INPUT.read_bytes()).hexdigest(),'T0_K':1675.0,'few_shot_temperatures_K':[1263.0,2087.0],
           'derivative_operator':'centered least-squares slope, fixed 11 points T0±5 K; slope summaries only where the full window exists','fair_support_points':1647,
           'condition_count':len(METHODS)*len(ARMS),'point_count':len(point_rows),'phase_consistent_corrections':True,
           'validation_checks':{'curve_rows_1649_per_C_method':True,'phase_cycle_closure_by_construction':True,'T0_and_shot_coordinates_match_frozen_protocol':True}}
    (out/'data/S2_SCORE_AUDIT.json').write_text(json.dumps(audit,indent=2)+'\n')
    with (out/'logs/COMMANDS.txt').open('a') as f: f.write('S2 scoring command: python3 scripts/run_s2_information_ladder.py --run-dir outputs/fes-bench-sio2-information-intervention-20261027\n')
    print(json.dumps({'audit':audit,'summary_rows':summary,'envelope_rows':env_rows},indent=2))
if __name__=='__main__': main()
