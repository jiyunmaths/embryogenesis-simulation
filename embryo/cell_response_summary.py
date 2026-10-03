"""Assess frozen pulse responses without fitting cell-type labels."""
import argparse
import json
from pathlib import Path
import numpy as np
from .resolution import write_json


def predictive_comparison(records):
    """Exploratory leave-developmental-seed-out prediction; no trial-level p values."""
    y=np.log([r['response'] for r in records]);seeds=np.array([r['seed'] for r in records])
    context=np.array([r['context'] for r in records]);chemical=np.array([r['chemical'] for r in records])
    result={}
    for name,x in [('context',context),('context_and_chemistry',np.column_stack((context,chemical)))]:
        pred=np.zeros(len(y));folds=[]
        for seed in sorted(set(seeds)):
            train=seeds!=seed;test=~train
            mean=x[train].mean(axis=0);scale=x[train].std(axis=0);scale[scale<1e-12]=1.
            xt=(x[train]-mean)/scale;xe=(x[test]-mean)/scale
            ym=y[train].mean();coef=np.linalg.solve(xt.T@xt+np.eye(x.shape[1]),xt.T@(y[train]-ym))
            pred[test]=ym+xe@coef
            denom=float(np.sum((y[test]-y[test].mean())**2))
            folds.append(dict(held_out_seed=int(seed),log_rmse=float(np.sqrt(np.mean((y[test]-pred[test])**2))),
                              r2=float(1-np.sum((y[test]-pred[test])**2)/denom) if denom>0 else None))
        result[name]=dict(pooled_out_of_seed_r2=float(1-np.sum((y-pred)**2)/np.sum((y-y.mean())**2)),folds=folds)
    return result


def summarize(root):
    root=Path(root);all_rows=[];per_graph=[];predictors=[];pair_ratios=[]
    for seed in (7,8,9):
        folder=root/f'seed-{seed}'
        if json.loads((folder/'status.json').read_text())['state']!='completed':raise ValueError('Wait for all seeds')
        data=json.loads((folder/'results.json').read_text())
        for key,r in data.items():
            if not r['numerical_pass']:raise ValueError('Numerical gate failed; do not pool unvalidated results')
            for row in r['trials']:all_rows.append(dict(seed=seed,graph_family=key,**row))
            per_graph.append(dict(seed=seed,graph_family=key,counts={name:sum(t['classification']==name for t in r['trials']) for name in sorted({t['classification'] for t in r['trials']})},
                                  max_reference_error=max(t['reference_max_log_error'] for t in r['trials']),max_radau_error=r['independent_max_log_error']))
        for arm in ('switch_on','keep_off'):
            p=data[arm+'_pattern'];u=data[arm+'_uniform']
            for i,cell_id in enumerate(p['context']['ids']):
                def response(r):
                    rows=[t for t in r['trials'] if t['cell_index']==i and t['factor'] in (.9,1.1)]
                    assert len(rows)==2
                    return float(np.mean([t['target_activator_log_auc_per_log_pulse'] for t in rows]))
                pr,ur=response(p),response(u);ctx=p['context']
                predictors.append(dict(seed=seed,arm=arm,cell_id=cell_id,response=pr,
                    context=[np.log(ctx['volume'][i]),ctx['exposure'][i],np.log(ctx['weighted_degree'][i]),ctx['radius'][i]],
                    chemical=np.log(np.array(p['initial'])[:,i]).tolist()))
                pair_ratios.append(dict(seed=seed,arm=arm,cell_id=cell_id,pattern_response=pr,uniform_response=ur,ratio=pr/ur))
    models=predictive_comparison(predictors)
    result=dict(trials=len(all_rows),developmental_histories=3,graphs=6,equilibrium_backgrounds=12,
        counts={name:sum(t['classification']==name for t in all_rows) for name in sorted({t['classification'] for t in all_rows})},
        target_recovery_censored=sum(t['target_recovery_time'] is None for t in all_rows),
        maximum_reference_log_error=max(r['reference_max_log_error'] for r in all_rows),
        maximum_radau_log_error=max(r['max_radau_error'] for r in per_graph),
        per_background=per_graph,predictive_comparison=models,matched_background_response=pair_ratios,
        interpretation='Responses are chemistry and context dependent in a shared reaction system. No identity labels, functional phenotype, autonomy, inheritance, or live mechanics are established. Regression is exploratory; paired uniform/pattern comparison changes whole-network chemistry.')
    write_json(root/'assessment.json',result)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(14,4.4))
    colors={7:'#2166ac',8:'#2b8c65',9:'#b05b31'}
    for seed in (7,8,9):
        rr=[r for r in predictors if r['seed']==seed]
        axes[0].scatter([r['chemical'][0] for r in rr],[r['response'] for r in rr],color=colors[seed],label=f'Seed {seed}',alpha=.7,s=24)
        pairs=[r for r in pair_ratios if r['seed']==seed]
        axes[1].scatter([r['uniform_response'] for r in pairs],[r['pattern_response'] for r in pairs],color=colors[seed],alpha=.7,s=24)
    axes[0].set(xlabel='Baseline ln(activator), patterned state',ylabel='Normalized activator response integral',title='A  Continuous response differences');axes[0].legend(frameon=False)
    lo=min(min(r['pattern_response'],r['uniform_response']) for r in pair_ratios);hi=max(max(r['pattern_response'],r['uniform_response']) for r in pair_ratios)
    axes[1].plot([lo,hi],[lo,hi],':',color='gray');axes[1].set(xlabel='Response on uniform background',ylabel='Response on patterned background',title='B  Same cell position and geometry')
    for family,color in [('pattern','#2166ac'),('uniform','#b05b31')]:
        xx=[];yy=[]
        for factor in (.5,.75,.9,1.1,1.25,1.5):
            rr=[r for r in all_rows if r['graph_family'].endswith('_'+family) and r['factor']==factor]
            recovery=[r['target_recovery_time'] for r in rr if r['target_recovery_time'] is not None]
            xx.append((factor-1)*100);yy.append(np.median(recovery) if recovery else np.nan)
        axes[2].plot(xx,yy,'o-',label=family,color=color)
    axes[2].set(xlabel='Activator pulse (%)',ylabel='Median target recovery time',title='C  Recovery after intervention');axes[2].legend(frameon=False)
    for ax in axes:ax.spines[['top','right']].set_visible(False)
    fig.suptitle('Cell response assay: frozen geometry, evolving chemistry; no prescribed identities',fontsize=13)
    fig.text(.5,.01,'A–B: mean of −10% and +10% normalized responses; 96 cell positions across 3 histories. C: recovered trials only; censoring reported in assessment.',ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.045,1,.95));fig.savefig(root/'response.png',dpi=200);plt.close(fig)
    print(json.dumps({k:v for k,v in result.items() if k not in ('matched_background_response','per_background')},indent=2))
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=Path('outputs/cell-response'))
    summarize(p.parse_args().output)
