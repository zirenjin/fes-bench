#!/usr/bin/env python3
"""Failure-first scientific landscape; build anew from explicit case contracts.

Each cell retains reported observations or an explicit NR/NA/RF/NE state.
No failure threshold, mean across systems, or missing-data imputation is used.
Colour uses the worst reported seed; text retains the seed range.
"""
from pathlib import Path
import json,csv
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,Normalize,LogNorm,to_rgba
from matplotlib.patches import Rectangle,Patch
from audit_panel_alignment import require_matplotlib_panel_alignment
R=Path(__file__).resolve().parents[1];OUT=R/'figures';Q=OUT/'qa';F=json.loads((R/'data/failure_logic.json').read_text())
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['DejaVu Sans'],'font.size':8,'axes.labelsize':8,'xtick.labelsize':8,'ytick.labelsize':8,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','pdf.fonttype':42})
cols=F['methods'];cases=F['cases'];cells={(c['case'],c['method']):c for c in F['cells']}
state_colors={'NR':'#F0F3F5','NA':'#F7F1E7','RF':'#DDE4ED','NE':'#E8F2F3'}
def num(v):return f'{v:.3g}'
def metric_cell(c,metric):
 if c['state']!='OBS':return c['state'],None,c['state']
 if metric=='events':
  if c['reference_events_max']==0:return 'NE',None,'NE'
  a=c['missed_events_min'];b=c['missed_events_max'];ref=c['reference_events_max']
  if a is None:return 'NR',None,'NR'
  return (f'{a}/{ref}' if a==b else f'{a}–{b}/{ref}'),b/ref,'OBS'
 lo=c[metric+'_min'];hi=c[metric+'_max']
 if lo is None:return 'NR',None,'NR'
 if metric=='accuracy':lo*=100;hi*=100
 text=num(lo) if lo==hi else num(lo)+('\n– ' if metric=='mae' else '–')+num(hi)
 return text,(lo/100 if metric=='accuracy' else hi),'OBS'
figs=[]
for metric in ['accuracy','mae','events']:
 fig,ax=plt.subplots(figsize=(183/25.4,5.1));fig.subplots_adjust(left=.29,right=.92,bottom=.18,top=.83)
 if metric=='accuracy':cmap=LinearSegmentedColormap.from_list('decision',['#D9A298','#FAF8F4','#81B4A5']);norm=Normalize(0,1);title='Stable-envelope accuracy (%)';ticks=[0,.5,1];ticklabels=['0','50','100']
 elif metric=='mae':cmap=LinearSegmentedColormap.from_list('error',['#ECF2F4','#7894A9','#C57C68']);vals=[c['mae_max'] for c in F['cells'] if c['mae_max'] is not None];assert min(vals)>0;norm=LogNorm(min(vals),max(vals));title='Mean pair ΔG MAE (meV/atom)';ticks=[.04,1,10,100];ticklabels=['0.04','1','10','100']
 else:cmap=LinearSegmentedColormap.from_list('miss',['#81B4A5','#D9A298']);norm=Normalize(0,1);title='Missed / reference envelope events';ticks=[0,1];ticklabels=['0','1']
 for i,case in enumerate(cases):
  for j,method in enumerate(cols):
   c=cells[(case['id'],method['id'])];text,value,state=metric_cell(c,metric);co=cmap(norm(value)) if value is not None else state_colors[state]
   rect=Rectangle((j-.5,i-.5),1,1,facecolor=co,edgecolor='white',linewidth=.8);ax.add_patch(rect)
   ax.text(j,i,text,ha='center',va='center',fontsize=7,color='#26323D')
 ax.set(xlim=(-.5,len(cols)-.5),ylim=(len(cases)-.5,-.5),yticks=range(len(cases)),yticklabels=[c['label'] for c in cases],xticks=range(len(cols)),xticklabels=[m['plot_label'] for m in cols])
 ax.tick_params(length=0,pad=8);ax.xaxis.tick_top();ax.tick_params(axis='x',labelsize=7.5)
 for spine in ax.spines.values():spine.set_visible(False)
 # Explicit divider between measured tail/outer-fold cases and open tasks.
 ax.axhline(5.5,color='#8FA4B1',lw=1)
 cax=fig.add_axes([.93,.31,.014,.39]);cb=fig.colorbar(mpl.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax);cb.set_ticks(ticks,labels=ticklabels);cb.ax.tick_params(labelsize=7);cb.outline.set_visible(False)
 fig.text(.29,.94,title,fontsize=10,fontweight='medium',ha='left')
 fig.legend(handles=[Patch(facecolor=state_colors[s],label=s+' · '+({'NR':'No matched result','NA':'Inputs not allowed','RF':'New reference missing','NE':'No reference event'}[s])) for s in ['NR','NA','RF','NE']],loc='lower left',bbox_to_anchor=(.07,.015),ncol=2,frameon=False,fontsize=7,columnspacing=2)
 key='failure_landscape_'+metric;fig.canvas.draw();require_matplotlib_panel_alignment(fig,json_out=Q/(key+'.alignment.json'),strict=True,tolerance_pt=1.5,exclude_axes=[cax])
 for ext in ['svg','pdf','png']:fig.savefig(OUT/(key+'.'+ext),dpi=300)
 csvpath=OUT/'source_data'/(key+'.csv');rows=[]
 for c in F['cells']:
  text,value,state=metric_cell(c,metric);rows.append({**{k:v for k,v in c.items() if k!='records'},'sources':';'.join(c['sources']),'display_metric':metric,'display_state':state,'display_value':text,'colour_value':value})
 with csvpath.open('w',newline='') as h:w=csv.DictWriter(h,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 caption='Each row keeps its own task, label budget and reference layer. Text shows seed extrema, not an uncertainty interval. Colour uses the worst reported seed. NR is no matched reported evaluation, NA is an input restriction, RF is a missing compatible new grid, and NE is no reference event. Reference limits are not model prediction errors. QH columns describe the declared numerical baselines, not certified QHA.'
 figs.append(dict(key=key,title='Failure landscape · '+title,caption=caption,section='failure-landscape',width_mm=183,height_mm=129.54,svg=f'figures/{key}.svg',pdf=f'figures/{key}.pdf',png=f'figures/{key}.png',csv=csvpath.relative_to(R).as_posix(),rows=len(rows)))
 plt.close(fig);print(key,flush=True)
meta={'figures':figs,'source':'data/failure_logic.json','schema':'Explicit observed and unscored states; no implicit model-failure labels.'}
(OUT/'failure_landscape_figures.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2));(OUT/'failure_landscape_figures.js').write_text('window.FES_FAILURE_FIGURES='+json.dumps(meta,ensure_ascii=False)+';\n')
