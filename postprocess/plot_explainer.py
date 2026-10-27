#!/usr/bin/env python3
"""Two source-backed figures for the English first-principles explainer.

Build anew; inherit only the existing Matplotlib typography and method palette.
Control values keep report precision. Hf uses a separate matched protocol.
The SiO2 example retains all four selected model/seed rows from the unified table.
"""
from pathlib import Path
import json
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from audit_panel_alignment import require_matplotlib_panel_alignment
import plot_evidence as p
R=Path(__file__).resolve().parents[1];E=json.loads((R/'data/explainer_content.json').read_text())
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['DejaVu Sans'],
 'font.size':8,'svg.fonttype':'none','pdf.fonttype':42})
def export(key,title,fig,rows,caption,section):
 fig.canvas.draw()
 require_matplotlib_panel_alignment(fig,json_out=p.QA/(key+'.alignment.json'),
  tolerance_pt=1.5,gutter_tolerance_pt=1.5,strict=True)
 for ext in ['svg','pdf','png']:
  fig.savefig(p.OUT/(key+'.'+ext),dpi=300)
 item={'key':key,'title':title,'caption':caption,'section':section,'width_mm':183,
  'height_mm':round(fig.get_figheight()*25.4,2),'svg':f'figures/{key}.svg',
  'pdf':f'figures/{key}.pdf','png':f'figures/{key}.png',
  'csv':p.csv_out(key,rows),'rows':len(rows)}
 p.FIGS.append(item);plt.close(fig)
# Pair each direct fit with its own residual control. Do not pool systems.
rows=E['controls'];fig,ax=plt.subplots(figsize=(p.W,3.55));fig.subplots_adjust(left=.20,right=.95,bottom=.21,top=.77)
for i,r in enumerate(rows):
 ax.plot([r['direct'],r['residual']],[i,i],color='#CBD6DC',lw=1.5,zorder=1)
 ax.scatter(r['direct'],i,s=32,color=p.C['analytic'],marker='o',edgecolors='white',lw=.5,zorder=3)
 ax.scatter(r['residual'],i,s=32,color=p.C['residual'],marker='s',edgecolors='white',lw=.5,zorder=3)
ax.set(yticks=range(len(rows)),yticklabels=[r['system']+(' *' if r['system']=='Hf' else '') for r in rows],ylim=(3.5,-.5),xlabel='Test ΔG MAE (meV/atom)');p.logx(ax,[v for r in rows for v in [r['direct'],r['residual']]])
ax.legend(handles=[p.Line2D([],[],marker='o',color=p.C['analytic'],linestyle='',label='Direct Thermo4'),p.Line2D([],[],marker='s',color=p.C['residual'],linestyle='',label='QH + residual Thermo4')],loc='lower center',bbox_to_anchor=(.5,1.08),ncol=2,frameon=False)
export('explainer_fit_controls','The same temperature function, with or without QH',fig,rows,'Each pair uses its own target labels and frozen temperature split. Hf (*) has a separate protocol. Values retain report precision. These deterministic fits do not define a cross-system uncertainty interval.','fitting')
# Four actual observations, no synthetic counterexample.
rs=[r for r in p.D['tail'] if r['system']=='sio2' and r['region']=='test' and r['model'] in ['T3_direct_G','temperature_phase_id_mlp_pairloss']]
rs=sorted(rs,key=lambda r:(int(r['seed']),r['model']!='temperature_phase_id_mlp_pairloss'))
assert len(rs)==4
fig,axes=plt.subplots(1,2,figsize=(p.W,3.6),sharey=True);fig.subplots_adjust(left=.24,right=.97,bottom=.21,top=.82,wspace=.32)
for i,r in enumerate(rs):
 for ax,key in zip(axes,['delta_G_MAE_meV_atom','stable_phase_accuracy']):ax.scatter(float(r[key]),i,s=33,color=p.color(r['model']),marker='D' if r['seed']=='23' else 'o',edgecolors='white',lw=.5,zorder=3)
axes[0].set(yticks=range(4),yticklabels=[('Phase-ID' if r['model'].startswith('temperature') else 'T3')+' · seed '+r['seed'] for r in rs],ylim=(3.5,-.5),xlim=(0,1.7),xlabel='Mean ΔG MAE (meV/atom)');axes[1].set(xlim=(-.06,1.06),xticks=[0,.5,1],xlabel='Stable-envelope accuracy')
for ax in axes:p.clean(ax)
p.panel(axes[0],'a','Average pair error');p.panel(axes[1],'b','Stable-phase decision')
export('explainer_decision_example','Lower mean error can still select the wrong phase',fig,rs,'Actual SiO₂ test-tail results, seeds 11 and 23. Each seed has a lower mean pair MAE for Phase-ID, but T3 selects the correct stable phase. The 330 temperatures are correlated evaluation points. This tail has no reference transition.','decisions')
existing=R/'figures/explainer_figures.json'
if existing.exists():p.FIGS.extend(x for x in json.loads(existing.read_text())['figures'] if x['key']=='hf_qh_sign_audit')
result={'figures':p.FIGS,'source':'Frozen evidence snapshot and directly linked diagnostic reports'}
(R/'figures/explainer_figures.json').write_text(json.dumps(result,ensure_ascii=False,indent=2));(R/'figures/explainer_figures.js').write_text('window.FES_EXPLAINER_FIGURES='+json.dumps(result,ensure_ascii=False)+';\n')
