#!/usr/bin/env python3
"""Reproducible, source-linked scientific figures for the frozen evidence snapshot.

Python/Matplotlib only. No synthetic observations or inferred confidence intervals.
All selectors partition source observations explicitly; raw source tables are retained.
"""
from pathlib import Path
from collections import Counter, defaultdict
import csv, json, re, sys, hashlib
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, MaxNLocator
from audit_panel_alignment import require_matplotlib_panel_alignment

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'figures'; QA=OUT/'qa'; DATA=OUT/'source_data'
for p in [OUT,QA,DATA]:p.mkdir(parents=True,exist_ok=True)
D=json.loads((ROOT/'data/review_data.json').read_text())
# Current score layer; the original snapshot and source hashes remain unchanged.
transfer_audit=ROOT/'data/transfer_score_audit.json'
if transfer_audit.exists():
 current=json.loads(transfer_audit.read_text());D['transferPair']=current['pair_metrics'];D['transferEnvelope']=current['envelope_metrics']
W=183/25.4
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['DejaVu Sans'],
 'font.size':8,'axes.labelsize':8,'axes.titlesize':9,'xtick.labelsize':8,'ytick.labelsize':8,
 'legend.fontsize':8,'axes.linewidth':.65,'axes.spines.top':False,'axes.spines.right':False,
 'text.color':'#26323D','axes.labelcolor':'#26323D','xtick.color':'#52606C','ytick.color':'#52606C',
 'svg.fonttype':'none','pdf.fonttype':42,'savefig.facecolor':'white','figure.facecolor':'white',
 'lines.linewidth':1.35,'lines.markersize':4,'axes.unicode_minus':True})
C={'analytic':'#637B7B','residual':'#A68B65','T2':'#BC7064','T3':'#547D9A','C-DPA':'#8877A5',
 'C-simple':'#3F8076','QH':'#A7AFB6','control':'#8A939B','ink':'#26323D','pale':'#EDF1F3','open':'#D9DFE3'}
LABELS={'thermo_phase_id_4term':'Thermo4','same_form_4term_analytic':'Thermo4 (matched)',
 'thermo_phase_id_3term':'Thermo3','thermo_form_3term':'Thermo3',
 'qh_residual_thermo_phase_id_4term':'QH + Thermo4','qh_residual_thermo_phase_id_3term':'QH + Thermo3',
 'temperature_phase_id_mlp_pairloss':'Phase-ID MLP','T3_direct_G':'T3: direct G',
 'T2_qh_residual':'T2: QH + residual','QH_only':'QH only',
 'T2_Hf_checkpoint_8000':'T2: selected 8000','T2_Hf_longrun_checkpoint_8000':'T2: selected 8000',
 'same_form_4term_T2_residual':'Residual Thermo4','thermo_form_3term_T2_residual':'Residual Thermo3',
 'C-DPA':'C-DPA','C-simple':'C-simple','QH_like_numerical':'QH-like control',
 'source_only_family_mean':'Source-family mean','zero_delta_G':'Zero gap'}
ORDER=['thermo_phase_id_4term','same_form_4term_analytic','thermo_phase_id_3term','thermo_form_3term',
 'qh_residual_thermo_phase_id_4term','qh_residual_thermo_phase_id_3term',
 'temperature_phase_id_mlp_pairloss','T3_direct_G','T2_qh_residual','T2_Hf_checkpoint_8000',
 'T2_Hf_longrun_checkpoint_8000','QH_only']
PHASE={'quartz_beta':'Quartz','cristobalite_beta':'Cristobalite','tridymite_p63mmc':'Tridymite',
 'bcc':'bcc','hcp':'hcp'}
SYS={'ti':'Ti','zr':'Zr','hf':'Hf','sio2':'SiO₂','SiO2':'SiO₂','CaSiO3':'CaSiO₃','MgSiO3':'MgSiO₃'}
FIGS=[]; VIEWS={}; EXCLUSIONS=[]; SELECT=set(sys.argv[1:])
def f(v):
 try:return float(v)
 except (ValueError,TypeError):return np.nan

def color(model):
 if model in ['C-DPA','C-simple']:return C[model]
 if model.startswith('T2'):return C['T2']
 if model.startswith('T3'):return C['T3']
 if model.startswith('qh_residual') or model.endswith('_T2_residual'):return C['residual']
 if model.startswith('thermo') or model=='same_form_4term_analytic':return C['analytic']
 if 'QH' in model:return C['QH']
 return C['control']
def label(m):return LABELS.get(m,m)
def slug(s):return re.sub(r'[^a-zA-Z0-9]+','_',str(s)).strip('_').lower()
def source_rows(name):
 s=next(x for x in D['sources'] if x['name']==name)
 with (ROOT/s['href']).open(encoding='utf-8-sig',newline='') as h:return list(csv.DictReader(h))
def csv_out(name,rows):
 path=DATA/(name+'.csv')
 if rows:
  fields=list(dict.fromkeys(k for r in rows for k in r))
  with path.open('w',newline='') as h:
   w=csv.DictWriter(h,fieldnames=fields);w.writeheader();w.writerows(rows)
 else:path.write_text('')
 return path.relative_to(ROOT).as_posix()
def panel(ax,letter,title):
 ax.set_title(title,loc='left',pad=13,fontweight='medium')
 ax.annotate(letter,xy=(0,1),xycoords='axes fraction',xytext=(-22,13),textcoords='offset points',fontsize=9,fontweight='bold',annotation_clip=False)
def clean(ax,grid=True):
 if grid:ax.grid(axis='x',color='#E7EBEE',lw=.55);ax.set_axisbelow(True)
 ax.tick_params(length=3,width=.6)
def logx(ax,values):
 assert all(v>0 for v in values if np.isfinite(v)), 'Log axis requires strictly positive observed values'
 ax.set_xscale('log');ax.xaxis.set_major_formatter(FuncFormatter(lambda v,p:f'{v:g}'))
 ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter());clean(ax)
def export(key,title,fig,rows,caption,section,view=None,exclude_axes=()):
 if SELECT and key not in SELECT:
  plt.close(fig);return key
 fig.canvas.draw()
 require_matplotlib_panel_alignment(fig,json_out=QA/(key+'.alignment.json'),
  tolerance_pt=1.5,gutter_tolerance_pt=1.5,strict=True,exclude_axes=exclude_axes)
 # Explicit final dimensions; no tight export that changes plot alignment.
 for ext in ['svg','pdf','png']:fig.savefig(OUT/(key+'.'+ext),dpi=300)
 item={'key':key,'title':title,'caption':caption,'section':section,'width_mm':183,
  'height_mm':round(fig.get_figheight()*25.4,2),'svg':f'figures/{key}.svg','pdf':f'figures/{key}.pdf',
  'png':f'figures/{key}.png','csv':csv_out(key,rows),'rows':len(rows)}
 FIGS.append(item)
 if view:VIEWS['|'.join(view)]=key
 plt.close(fig)
 print(key,flush=True)
 return key

def run():
 # Explicit claim-to-experiment relationships, not a strength score.
 fig,ax=plt.subplots(figsize=(W,3.8));fig.subplots_adjust(left=.29,right=.97,bottom=.20,top=.87)
 nr=len(D['claims']);nc=len(D['cards']);links=[]
 claimnames=['Transfer question','Temperature-tail controls','Two transfer folds','Metrics / events',
             'QH failure attribution','Scoped carbon QHA','Eight-system scope']
 for i,c in enumerate(D['claims']):
  for eid in c['cards']:
   j=int(eid[1:])-1;ax.scatter(j,i,s=43,color=C['C-simple'],edgecolors='white',linewidths=.5,zorder=3)
   links.append({'claim':c['title'],'evidence_state':c['state'],'experiment_id':eid})
 ax.set(xticks=range(nc),xticklabels=[x['id'] for x in D['cards']],yticks=range(nr),yticklabels=claimnames[:nr],xlim=(-.65,nc-.35),ylim=(nr-.5,-.5))
 ax.tick_params(axis='x',rotation=90,length=0,pad=8);ax.tick_params(axis='y',length=0,pad=8)
 ax.set_xticks(np.arange(-.5,nc,1),minor=True);ax.set_yticks(np.arange(-.5,nr,1),minor=True)
 ax.grid(which='minor',color='#E8ECEF',lw=.5);ax.tick_params(which='minor',length=0)
 for sp in ax.spines.values():sp.set_visible(False)
 ax.set_title('Every headline claim points to named experiments',loc='left',pad=16)
 export('evidence_map','主张与实验证据的对应关系',fig,links,'圆点表示当前登记中明确关联的实验；不代表证据强度评分。E001–E019 可在实验档案直接查看。','claims')

 # Reference layers kept separate. Zero means no row of this type in the frozen manifest.
 vals=[];rows=[]
 for s in D['systems']:
  pair=s['layers'].get('derived_source_deltaG_candidate',0)+s['layers'].get('source_deltaG_label',0)+s['layers'].get('derived_source_eos_DeltaG_candidate',0)
  newpair=sum(r['system']==s['system'] and r['admission_state'].startswith('ADMITTED') for r in D['carbon'])
  vals.append([s['sourceG'],pair,s['layers'].get('qha_comparison_G_label',0),s['acceptedG'],newpair,s['openG']])
  rows.append(dict(system=s['system'],source_G=s['sourceG'],source_pair_candidate_or_label=pair,QHA_comparison_G=s['layers'].get('qha_comparison_G_label',0),new_admitted_G=s['acceptedG'],new_admitted_pair=newpair,open_target_G=s['openG']))
 fig,ax=plt.subplots(figsize=(W,4.3));fig.subplots_adjust(left=.18,right=.98,bottom=.19,top=.85)
 a=np.array(vals,float);ax.imshow(np.log1p(a),cmap=LinearSegmentedColormap.from_list('inventory',['#F4F6F7','#B1C6CD']),aspect='auto',vmin=0,vmax=np.log1p(a.max()))
 for i in range(8):
  for j in range(6):ax.text(j,i,f'{int(a[i,j]):,}',ha='center',va='center',fontsize=8,color=C['ink'])
 ax.set(yticks=range(8),yticklabels=[SYS.get(s['system'],s['system']) for s in D['systems']],xticks=range(6),xticklabels=['Source G','Source pair\n(candidate/label)','QHA control G','New G\n(admitted)','New ΔG\n(admitted)','Target G\n(open)'])
 ax.tick_params(length=0,pad=10);ax.set_title('Reference inventory: counts retain their own layer and unit',loc='left',pad=18)
 for sp in ax.spines.values():sp.set_visible(False)
 export('reference_inventory','八体系参考层与标签缺口矩阵',fig,rows,'每格为对应类别的行／状态数，底色按 log(1+计数) 辅助阅读。源 pair 混有候选与已准入标签；有数据不等于通过合同。零表示该清单类别无记录，不表示体系无自由能。','coverage')
 fig,ax=plt.subplots(figsize=(W,3.1));fig.subplots_adjust(left=.18,right=.94,bottom=.22,top=.84)
 required=[s for s in D['systems'] if s['requiredG']>0]
 y=np.arange(len(required));accepted=[s['acceptedG'] for s in required];opened=[s['openG'] for s in required]
 ax.barh(y,accepted,color=C['C-simple'],height=.48,label='Admitted target G')
 ax.barh(y,opened,left=accepted,color=C['open'],height=.48,label='Open target G')
 ax.set(yticks=y,yticklabels=[SYS.get(s['system'],s['system']) for s in required],xlabel='Number of required phase-G states',xlim=(0,36));ax.invert_yaxis();clean(ax,False)
 for i,s in enumerate(required):ax.text(s['requiredG']+1,i,f"{s['acceptedG']}/{s['requiredG']}",va='center',fontsize=8)
 ax.legend(loc='lower left',bbox_to_anchor=(0,1.03),ncol=2,frameon=False)
 export('coverage_targets','新增相自由能目标的已接入／开放状态',fig,[dict(system=s['system'],admitted=s['acceptedG'],open=s['openG'],required=s['requiredG']) for s in required],'仅展示主清单明确列出新增 G 需求的四体系，标签状态数不是总工作量百分比。SiO₂/Ti/Zr/Hf 在该类别中无新增需求行。','coverage')

 # All source rows are partitioned by exactly the existing selector dimensions.
 groups=defaultdict(list)
 for r in D['tail']:groups[(r['protocol_id'],r['system'],r['region'])].append(r)
 for (protocol,system,region),rs in groups.items():
  models=[m for m in ORDER if any(r['model']==m for r in rs)]
  fig,axes=plt.subplots(1,2,figsize=(W,max(3.3,.34*len(models)+1.45)),sharey=True)
  fig.subplots_adjust(left=.29,right=.96,bottom=.19,top=.80,wspace=.36)
  for i,m in enumerate(models):
   mr=[r for r in rs if r['model']==m]
   for j,r in enumerate(mr):
    yy=i+(j-(len(mr)-1)/2)*.17;mark={'11':'o','23':'D','37':'^'}.get(r['seed'],'o')
    axes[0].scatter(f(r['delta_G_MAE_meV_atom']),yy,s=23,marker=mark,color=color(m),edgecolors='white',lw=.45,zorder=3)
    axes[1].scatter(f(r['stable_phase_accuracy']),yy,s=23,marker=mark,color=color(m),edgecolors='white',lw=.45,zorder=3)
  axes[0].set(yticks=range(len(models)),yticklabels=[label(m) for m in models],ylim=(len(models)-.5,-.5),xlabel='ΔG MAE (meV/atom)')
  logx(axes[0],[f(r['delta_G_MAE_meV_atom']) for r in rs])
  axes[1].set(xlim=(-.06,1.06),xticks=[0,.5,1],xlabel='Stable-phase accuracy');clean(axes[1]);axes[1].spines['left'].set_visible(False);axes[1].tick_params(axis='y',length=0)
  panel(axes[0],'a','Relative free energy');panel(axes[1],'b','Phase decision')
  fig.suptitle(f"{SYS.get(system,system)} · {region} · "+('batch 8 / 3000 steps' if protocol.startswith('cross') else 'Hf protocol'),fontsize=9,y=.98)
  cap=f'{system.upper()}／{region}，协议 {protocol}。圆点为 seed 11 或确定性控制，菱形为 seed 23，三角为 seed 37；展示所有归档结果，不将 seed 散布当置信区间。测试无参考根时 Tc 不适用。'
  export('tail_'+slug(protocol+'_'+system+'_'+region),'同协议模型比较 · '+system+' / '+region,fig,rs,cap,'models',['tail',protocol,system,region])
 # All models in the matched three-system test protocol. Missing methods remain blank.
 rs=[r for r in D['tail'] if r['protocol_id'].startswith('cross_system') and r['region']=='test']
 models=[m for m in ORDER if any(r['model']==m for r in rs)]
 fig,axes=plt.subplots(1,3,figsize=(W,4.2),sharex=True,sharey=True)
 fig.subplots_adjust(left=.25,right=.98,bottom=.18,top=.84,wspace=.18)
 low=min(f(r['delta_G_MAE_meV_atom']) for r in rs);high=max(f(r['delta_G_MAE_meV_atom']) for r in rs)
 for k,system in enumerate(['ti','zr','sio2']):
  ax=axes[k]
  for i,m in enumerate(models):
   mr=[r for r in rs if r['system']==system and r['model']==m]
   for j,r in enumerate(mr):ax.scatter(f(r['delta_G_MAE_meV_atom']),i+(j-(len(mr)-1)/2)*.17,s=22,color=color(m),marker={'11':'o','23':'D','37':'^'}.get(r['seed'],'o'),edgecolors='white',lw=.4,zorder=3)
  logx(ax,[f(r['delta_G_MAE_meV_atom']) for r in rs]);ax.set_xlim(low*.55,high*1.6);ax.set_ylim(len(models)-.5,-.5);ax.set_yticks(range(len(models)),[label(m) for m in models]);ax.set_xlabel('ΔG MAE (meV/atom)');panel(ax,chr(97+k),SYS[system])
  if k:ax.tick_params(axis='y',length=0)
 export('tail_overview','匹配测试尾部：三体系完整模型对照',fig,rs,'同一 batch=8、3000 步、60/20/20 尾部协议；每个点保留其 seed，空白表示该协议无归档结果。各体系同横轴；只比较相同监督信息。','models')

 # Transfer: matched pairs (no source diagnostic in headline comparison).
 target=[r for r in D['transferPair'] if r['region']=='outer_test' and r['pair_role']=='target_related_pair' and r['method'] in ['C-DPA','C-simple']]
 pairs=list(dict.fromkeys((r['condition'],r['phase_a'],r['phase_b']) for r in target))
 fig,axes=plt.subplots(1,2,figsize=(W,3.7));fig.subplots_adjust(left=.28,right=.96,bottom=.22,top=.76,wspace=.38)
 names=[]
 for i,p in enumerate(pairs):
  rr=[r for r in target if (r['condition'],r['phase_a'],r['phase_b'])==p]
  names.append(('SiO₂' if p[0].startswith('sio2') else 'Ti')+' · '+PHASE.get(p[1],p[1])+' − '+PHASE.get(p[2],p[2]))
  xx=[f(r['delta_G_MAE_eV_atom'])*1000 for r in rr];axes[0].plot(xx,[i]*len(xx),color='#CDD5DA',lw=1.2,zorder=1)
  for r,v in zip(rr,xx):axes[0].scatter(v,i,color=color(r['method']),marker='s' if r['method']=='C-DPA' else 'o',s=32,zorder=3,edgecolors='white',lw=.5)
 axes[0].set(yticks=range(len(pairs)),yticklabels=names,ylim=(len(pairs)-.45,-.55),xlabel='Target-pair ΔG MAE (meV/atom)');logx(axes[0],[f(r['delta_G_MAE_eV_atom'])*1000 for r in target])
 env=[r for r in D['transferEnvelope'] if r['region']=='outer_test' and r['method'] in ['C-DPA','C-simple']]
 folds=list(dict.fromkeys(r['condition'] for r in env))
 for i,fold in enumerate(folds):
  rr=[r for r in env if r['condition']==fold];xx=[f(r['stable_envelope_accuracy']) for r in rr];axes[1].plot(xx,[i]*len(xx),color='#CDD5DA',lw=1.2)
  for r,v in zip(rr,xx):axes[1].scatter(v,i,color=color(r['method']),marker='s' if r['method']=='C-DPA' else 'o',s=32,edgecolors='white',lw=.5,zorder=3)
 axes[1].set(yticks=range(len(folds)),yticklabels=['SiO₂ LOPO','Ti LOSO'],ylim=(len(folds)-.5,-.5),xlim=(0,1.07),xticks=[0,.5,1],xlabel='Envelope accuracy');clean(axes[1]);panel(axes[0],'a','Matched target pairs');panel(axes[1],'b','Stable envelope')
 fig.legend(handles=[Line2D([],[],color=C[m],marker='s' if m=='C-DPA' else 'o',linestyle='',label=m) for m in ['C-DPA','C-simple']],loc='upper center',bbox_to_anchor=(.6,.99),ncol=2,frameon=False)
 export('transfer_overview','两个冻结外折：结构表示的直接比较',fig,target,'只取 outer_test 中 target_related_pair；SiO₂ 两相对和 Ti 一相对分别展示，不把源相诊断当迁移成绩。右图为同外折稳定下包络准确率；两者是不同推断层级。','models')
 tg=defaultdict(list)
 for r in D['transferPair']:tg[(r['condition'],r['region'])].append(r)
 for (fold,region),rs in tg.items():
  fig,axes=plt.subplots(1,2,figsize=(W,max(3.5,len(rs)*.25+1.8)),sharey=True);fig.subplots_adjust(left=.38,right=.96,bottom=.17,top=.84,wspace=.35)
  names=[]
  for i,r in enumerate(rs):
   names.append(label(r['method'])+' · '+PHASE.get(r['phase_a'],r['phase_a'])+'−'+PHASE.get(r['phase_b'],r['phase_b']))
   axes[0].scatter(f(r['delta_G_MAE_eV_atom'])*1000,i,s=25,color=color(r['method']),edgecolors='white',lw=.45)
   axes[1].scatter(f(r['sign_accuracy']),i,s=25,color=color(r['method']),edgecolors='white',lw=.45)
  axes[0].set(yticks=range(len(rs)),yticklabels=names,ylim=(len(rs)-.5,-.5),xlabel='Pair ΔG MAE (meV/atom)');logx(axes[0],[f(r['delta_G_MAE_eV_atom'])*1000 for r in rs]);axes[1].set(xlim=(-.06,1.06),xticks=[0,.5,1],xlabel='Pair sign accuracy');clean(axes[1]);axes[1].tick_params(axis='y',length=0);panel(axes[0],'a','Relative free energy');panel(axes[1],'b','Pairwise phase sign')
  export('transfer_'+slug(fold+'_'+region),'迁移明细 · '+fold+' / '+region,fig,rs,'该外折／数据域全部方法与相对均展示。源相对诊断与目标相对的角色在下方原始表中标明；事件结果单独列出，未命中时 Tc 误差保持不适用。','models',['transfer',fold,region])

 def bound(r):
  m=re.search(r'(?:total|sum|envelope|bound)[^\d]*([\d.]+)\s*meV',r.get('convergence_error',''),re.I)
  return float(m[1]) if m else np.nan
 carbon=[dict(r,delta_G_meV_atom=f(r['label_value'])*1000,numerical_bound_meV_atom=bound(r)) for r in D['carbon']]
 pressures=sorted(set(f(r['P_GPa']) for r in carbon));temps=sorted(set(f(r['T_K']) for r in carbon))
 matrix=np.array([[next(r['delta_G_meV_atom'] for r in carbon if f(r['P_GPa'])==p and f(r['T_K'])==t) for t in temps] for p in pressures])
 fig,ax=plt.subplots(figsize=(W,3.8));fig.subplots_adjust(left=.14,right=.83,bottom=.20,top=.84)
 cmap=LinearSegmentedColormap.from_list('phase_gap',['#B66D5E','#F8F7F4','#547C97'])
 im=ax.imshow(matrix,cmap=cmap,norm=TwoSlopeNorm(vmin=-150,vcenter=0,vmax=150),aspect='auto')
 for i,p in enumerate(pressures):
  for j,t in enumerate(temps):
   r=next(r for r in carbon if f(r['P_GPa'])==p and f(r['T_K'])==t);v=r['delta_G_meV_atom'];b=r['numerical_bound_meV_atom']
   ax.text(j,i,f'{v:+.1f}',ha='center',va='center',color='white' if abs(v)>90 else C['ink'],fontsize=9)
   if abs(v)<=b:ax.add_patch(mpl.patches.Rectangle((j-.46,i-.43),.92,.86,fill=False,edgecolor=C['ink'],lw=1.3,linestyle='--'))
 ax.set(xticks=range(4),xticklabels=[int(t) for t in temps],yticks=range(4),yticklabels=[int(p) for p in pressures],xlabel='Temperature (K)',ylabel='Pressure (GPa)');ax.tick_params(length=0,pad=8)
 for sp in ax.spines.values():sp.set_visible(False)
 cax=fig.add_axes([.87,.22,.022,.57]);cb=fig.colorbar(im,cax=cax);cb.set_label('ΔG (meV/atom)');cb.set_ticks([-150,0,150])
 export('carbon_pt_grid','碳双相 QHA：温压网格与符号可辨识性',fig,carbon,'ΔG = G石墨 − G金刚石；负值偏向石墨。虚线框表示 |ΔG| 不超过归档数值界，符号未分辨；该界不包含全部系统误差，图不能视为完整相图。所有16个状态保留。','physics',exclude_axes=[cax])
 for p in pressures:
  rs=sorted([r for r in carbon if f(r['P_GPa'])==p],key=lambda r:f(r['T_K']))
  for show in [True,False]:
   fig,ax=plt.subplots(figsize=(W,3.5));fig.subplots_adjust(left=.14,right=.96,bottom=.21,top=.84)
   xx=np.array([f(r['T_K']) for r in rs]);yy=np.array([r['delta_G_meV_atom'] for r in rs]);bb=np.array([r['numerical_bound_meV_atom'] for r in rs]);co=C['C-simple'] if p<=5 else C['T3']
   ax.axhline(0,color='#BFC7CC',lw=.8,ls=(0,(3,3)))
   if show:ax.errorbar(xx,yy,yerr=bb,fmt='o-',color=co,lw=1.5,ms=4.5,capsize=3,elinewidth=1)
   else:ax.plot(xx,yy,'o-',color=co,lw=1.5,ms=4.5)
   ax.set(xlabel='Temperature (K)',ylabel='G(graphite) − G(diamond) (meV/atom)',xticks=temps);ax.set_title(f'Carbon · {p:g} GPa · scoped quantum-harmonic QHA',loc='left',pad=14);clean(ax,False)
   export('carbon_'+slug(p)+'_'+str(int(show)),'碳配对自由能 · '+f'{p:g} GPa',fig,rs,'四个温度点全部保留。点间线仅连接离散状态；误差条为归档数值敏感度，不是置信区间，也不覆盖泛函／超胞／非谐误差。','physics',['carbon',str(int(p)),str(int(show))])
 # Ti joint physical outputs, all temperatures and pressures.
 for pressure in ['all']+list(dict.fromkeys(r['P_GPa'] for r in D['tiG'])):
  rs=D['tiG'] if pressure=='all' else [r for r in D['tiG'] if r['P_GPa']==pressure]
  fig,axes=plt.subplots(1,2,figsize=(W,3.5));fig.subplots_adjust(left=.12,right=.97,bottom=.22,top=.80,wspace=.38)
  ps=list(dict.fromkeys(r['P_GPa'] for r in rs));cols=plt.colormaps['Blues'](np.linspace(.4,.9,len(ps)))
  for p,co in zip(ps,cols):
   rr=sorted([r for r in rs if r['P_GPa']==p],key=lambda r:f(r['T_K']))
   xx=[f(r['T_K']) for r in rr]
   axes[0].plot(xx,[f(r['G_eV_per_atom']) for r in rr],color=co,label=f'{float(p):g} GPa')
   axes[1].plot(xx,[f(r['Vmin_A3_per_atom']) for r in rr],color=co)
  axes[0].set(xlabel='Temperature (K)',ylabel='G (eV/atom)');axes[1].set(xlabel='Temperature (K)',ylabel='Equilibrium volume (Å³/atom)')
  for ax in axes:clean(ax,False);ax.xaxis.set_major_locator(MaxNLocator(4))
  panel(axes[0],'a','Single-phase free energy');panel(axes[1],'b','Volume response')
  fig.legend(*axes[0].get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.55,.99),ncol=5,frameon=False)
  export('ti_'+slug(pressure),'Ti hcp 单相 QHA：自由能与体积响应',fig,rs,'全部所选压力的31个温度点；固定 c/a、Te=0、15³ q 网格。155条 G 是单相对照，不支持 hcp–bcc ΔG 或相变温度。','physics',['ti',pressure])
 # Hf historical audit retains the complete 1252-temperature curve.
 for mode in ['prediction','components']:
  rs=sorted(D['hfComponent'],key=lambda r:f(r['T_K']))
  defs=[('deltaG_ref_bcc_minus_hcp_eV_atom','Reference',C['ink']),('deltaG_T2_eV_atom','T2 historical',C['T2']),('deltaG_thermo4_residual_eV_atom','QH + Thermo4',C['C-simple'])] if mode=='prediction' else [('deltaG_ref_bcc_minus_hcp_eV_atom','Reference ΔG',C['ink']),('delta_FQH_bcc_minus_hcp_eV_atom','ΔF_QH',C['QH']),('delta_R_target_bcc_minus_hcp_eV_atom','Target ΔR',C['T3']),('delta_correction_bcc_minus_hcp_eV_atom','Learned correction',C['T2'])]
  fig,ax=plt.subplots(figsize=(W,4.0));fig.subplots_adjust(left=.15,right=.96,bottom=.19,top=.77)
  domains=list(dict.fromkeys(r['domain'] for r in rs))
  for i,domain in enumerate(domains):
   xx=[f(r['T_K']) for r in rs if r['domain']==domain]
   ax.axvspan(min(xx),max(xx),color=['#F4F6F8','#EDF2EF','#F8F1EB'][i%3],zorder=0)
   ax.text((min(xx)+max(xx))/2,1.025,domain,transform=ax.get_xaxis_transform(),ha='center',fontsize=8,color='#64717D')
  for key,name,co in defs:ax.plot([f(r['T_K']) for r in rs],[f(r[key])*1000 for r in rs],label=name,color=co,lw=1.5 if co!=C['QH'] else 1.1,ls='--' if co==C['QH'] else '-')
  ax.set(xlabel='Temperature (K)',ylabel='Relative energy (meV/atom)');ax.legend(loc='lower center',bbox_to_anchor=(.5,1.15),ncol=2,frameon=False);clean(ax,False)
  export('hf_'+mode,'Hf 历史检查点：'+('重建预测' if mode=='prediction' else '逐分量诊断'),fig,rs,'Hf、seed 11、1500步历史检查点；所有1252温度行展示。背景为原始train/validation/test。此检查点初始化设置曾有错配，不能代表后续模型最优能力。','physics',['hf',mode])
 # Derivative sensitivity: source means and seed ranges, all groups/windows retained.
 shape=source_rows('RELATIVE_SHAPE_WINDOW_SENSITIVITY_20261020.csv');sg=defaultdict(list)
 for r in shape:
  r=dict(r,phase_pair_raw=r['phase_pair'])
  if r['phase_pair']=='bcc_minus_hcp':r['phase_pair']='bcc−hcp'
  sg[(r['protocol_id'],r['system'],r['region'],r['phase_pair'])].append(r)
 # Preserve existing asset identities; 19–21 were alias-only views, now merged into 16–18.
 old_options=json.loads((OUT/'figure_manifest.json').read_text())['shape_options']
 stable_keys={(x['protocol'],x['system'],x['region'],x['pair']):x['key'] for x in old_options if x['pair']!='bcc_minus_hcp'}
 shape_options=[]
 for idx,(g,rs) in enumerate(sg.items()):
  key=stable_keys.get(g,'shape_'+str(idx+1).zfill(2));assert len({(r['model'],r['smoother_window_points']) for r in rs})==len(rs);shape_options.append({'key':key,'protocol':g[0],'system':g[1],'region':g[2],'pair':g[3]})
  fig,axes=plt.subplots(1,2,figsize=(W,3.9));fig.subplots_adjust(left=.15,right=.97,bottom=.31,top=.84,wspace=.4)
  for m in dict.fromkeys(r['model'] for r in rs):
   rr=sorted([r for r in rs if r['model']==m],key=lambda r:f(r['smoother_window_points']))
   for ax,base in zip(axes,['relative_entropy_slope_MAE_meV_atom_K','delta_G_curvature_MAE_microeV_atom_K2']):
    xx=np.array([f(r['smoother_window_points']) for r in rr]);yy=np.array([f(r[base+'_mean']) for r in rr]);lo=np.array([f(r[base+'_min']) for r in rr]);hi=np.array([f(r[base+'_max']) for r in rr])
    ax.plot(xx,yy,color=color(m),marker='s' if '3term' in m else 'o',linestyle='--' if '3term' in m else '-',ms=3.2,label=label(m));ax.fill_between(xx,lo,hi,color=color(m),alpha=.1,lw=0)
  axes[0].set(ylabel='Slope MAE (meV/atom/K)');axes[1].set(ylabel='Curvature MAE (μeV/atom/K²)')
  for ax in axes:
   ax.set(xlabel='Smoother window (temperature points)',xticks=sorted(set(f(r['smoother_window_points']) for r in rs)));ax.set_yscale('log');ax.yaxis.set_major_formatter(FuncFormatter(lambda v,p:f'{v:g}'));ax.yaxis.set_minor_formatter(mpl.ticker.NullFormatter());clean(ax,False)
  panel(axes[0],'a','Relative slope');panel(axes[1],'b','Curvature sensitivity')
  fig.legend(*axes[0].get_legend_handles_labels(),loc='lower center',bbox_to_anchor=(.52,.005),ncol=3,frameon=False,columnspacing=1.5,handlelength=1.8)
  export(key,'相对斜率与曲率的平滑窗口敏感性',fig,rs,'每条线为归档模型的均值，淡带为归档 seed 最小–最大范围；不是预测区间。保留11/21/41窗口，导数单位分开标注。不同窗口不能用来挑选更好排名或物理因果归因。','diagnostics',['shape',key])
 # Job categories are an exact partition of ledger ROWS, not unique jobs or results.
 def statecat(s):
  if 'FAIL' in s or s in ['TIMEOUT','OUT_OF_MEMORY','CANCELLED','STOPPED']:return 'Failed / interrupted'
  if s.startswith('COMPLETED') or s.startswith('COMPUTED_'):return 'Completed / computed'
  if 'RUNNING' in s or 'PENDING' in s or 'PREPARED' in s or 'READY' in s:return 'Active / queued / staged'
  return 'Other / not admitted'
 cats=['Completed / computed','Failed / interrupted','Active / queued / staged','Other / not admitted']
 systems=list(dict.fromkeys(r['system'] for r in D['jobs']));counts=Counter((r['system'],statecat(r['state'])) for r in D['jobs'])
 fig,ax=plt.subplots(figsize=(W,max(3.5,len(systems)*.3+1.7)));fig.subplots_adjust(left=.23,right=.97,bottom=.18,top=.80)
 left=np.zeros(len(systems));cols=[C['T3'],C['T2'],C['residual'],C['open']]
 for cat,co in zip(cats,cols):
  values=np.array([counts[(s,cat)] for s in systems]);ax.barh(range(len(systems)),values,left=left,color=co,height=.55,label=cat);left+=values
 ax.set(yticks=range(len(systems)),yticklabels=systems,xlabel='Ledger rows (not unique jobs or accepted labels)');ax.invert_yaxis();ax.xaxis.set_major_locator(MaxNLocator(integer=True));clean(ax)
 ax.legend(loc='lower left',bbox_to_anchor=(-.06,1.05),ncol=2,frameon=False)
 export('compute_outcomes','计算账本：已记录的任务状态',fig,[dict(r,display_state_category=statecat(r['state'])) for r in D['jobs']],'全部253行账本按原state确定性分类。重复续跑／审计记录均保留，因此不是唯一作业数；COMPLETED不表示自由能已准入，状态也不代表实时队列。','jobs')
 # Exact pending slots, no percentage of manuscript readiness inferred.
 fig,ax=plt.subplots(figsize=(W,4.3));fig.subplots_adjust(left=.30,right=.95,bottom=.17,top=.87)
 ns=[];counts=[];pr=[]
 for i,r in enumerate(D['manuscript']):
  name=r['file'].replace('.tex','');ns.append(name);counts.append(len(r['pending']));pr.append({'section_file':r['file'],'pending_occurrences':len(r['pending']),'pending_ids':';'.join(r['pending'])})
 ax.barh(range(len(ns)),counts,color='#9CACB9',height=.55)
 for i,v in enumerate(counts):ax.text(v+.35,i,str(v),va='center',fontsize=8)
 ax.set(yticks=range(len(ns)),yticklabels=ns,xlabel='Explicit pending / result-slot occurrences',xlim=(0,max(counts)*1.2+1));ax.invert_yaxis();ax.xaxis.set_major_locator(MaxNLocator(integer=True));clean(ax,False)
 export('manuscript_gaps','论文章节：显式待填结果栏位',fig,pr,'逐文件统计 pending/resultslot 标记出现次数；不是独立未完成实验数，也不是论文完成度。所有11个章节文件保留，包括零标记章节。','paper')
 if SELECT:
  old=json.loads((OUT/'figure_manifest.json').read_text());updated={x['key']:x for x in FIGS};retired={'shape_19','shape_20','shape_21'};FIGS[:]=[updated.get(x['key'],x) for x in old['figures'] if x['key'] not in retired];VIEWS.update({k:v for k,v in old['views'].items() if v not in retired})
 result={'backend':'Python / Matplotlib','snapshot_sha256':hashlib.sha256((ROOT/'data/review_data.json').read_bytes()).hexdigest(),
  'palette':C,'figures':FIGS,'views':VIEWS,'shape_options':shape_options,'exclusions':EXCLUSIONS,
  'data_policy':'All existing observations preserved within explicit selector partitions; no synthetic data; no inferred confidence interval.'}
 (OUT/'figure_manifest.json').write_text(json.dumps(result,ensure_ascii=False,indent=2))
 (OUT/'figure_manifest.js').write_text('window.FES_FIGURES='+json.dumps(result,ensure_ascii=False)+';\n')
 print(f'Exported {len(FIGS)} scientific figures',flush=True)
if __name__=='__main__':run()
