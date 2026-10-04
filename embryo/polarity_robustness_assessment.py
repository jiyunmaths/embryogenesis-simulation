"""Verify and summarize the fixed-ratio polarity study without altering its gates.

Recompute recorded spectra, contrast and refinement metrics. Preserve the
declared quantitative failure separately from concordant qualitative outcomes.
"""
import argparse
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from .feedback_long import digest
from .neighbor_context import read, validate_graph
from .parameter_robustness import frozen_spectrum, verify
from .polarity_robustness import checked_history, first_crossing, refinement_comparison
from .resolution import write_json


def contact_metrics(row):
    delta,masses=validate_graph(row['delta'],row['volumes'])
    conductance=masses[:,None]*delta
    indices=np.triu_indices(len(masses),1);g=conductance[indices]
    centers=np.array(row['centers']);distance=np.linalg.norm(centers[indices[0]]-centers[indices[1]],axis=1)
    positive=g>0;weights=g[positive]
    return dict(conductance_sum=float(weights.sum()),edges=int(positive.sum()),
        positive_edge_weight_cv=float(weights.std()/weights.mean()),
        area_proxy_sum=float(np.sum(g*distance)),mean_retained_centroid_distance=float(distance[positive].mean()))


def exceeded_intervals(times,error,limit):
    mask=np.asarray(error)>limit
    starts=np.flatnonzero(mask & ~np.r_[False,mask[:-1]])
    ends=np.flatnonzero(mask & ~np.r_[mask[1:],False])
    return [dict(first=float(times[a]),last=float(times[b]),observations=int(b-a+1)) for a,b in zip(starts,ends)]


def derive(root):
    root=Path(root).resolve();p=read(root/'protocol.json');verify(p);ph=digest(root/'protocol.json')
    status=read(root/'status.json')
    if status['state'] not in ('completed','completed_with_unresolved_checks') or status['completed']!=len(p['jobs']):
        raise ValueError('Requires a finished study, including explicitly unresolved gates')
    evidence={str(root/'protocol.json'):ph,str(root/'status.json'):digest(root/'status.json')}
    def include(path,expected=None):
        path=Path(path);actual=digest(path)
        if expected is not None and actual!=expected:raise ValueError('Changed evidence: '+str(path))
        evidence[str(path.resolve())]=actual
    refs=read(root/'frozen/results.json');include(root/'frozen/results.json')
    if refs['protocol_sha256']!=ph:raise ValueError('Changed frozen protocol')
    for path,sha in refs['paths_sha256'].items():include(path,sha)
    rows=[];histories={};prefixes=[];endpoints=[]
    spectrum_error=contrast_error=0.;spectrum_observations=0
    for job in p['jobs']:
        folder=root/job['key'];result=read(folder/'result.json');h=read(folder/'history.json')
        if result['protocol_sha256']!=ph or result['job']!=job or not result['quality_pass']:
            raise ValueError('Invalid moving result')
        include(folder/'result.json');include(folder/'history.json',result['history_sha256'])
        include(folder/'latest_state.npz',result['checkpoint_sha256']);checked_history(h,p,job,p['duration'])
        cp=read(folder/'prefix/comparison.json');include(folder/'prefix/comparison.json')
        if cp['protocol_sha256']!=ph or cp['job']!=job or not cp['passed']:raise ValueError('Invalid context gate')
        for path,sha in cp['evidence_sha256'].items():include(path,sha)
        prefixes.append(cp)
        ep=read(folder/'endpoint/protocol.json');verify(ep);include(folder/'endpoint/protocol.json')
        e=read(folder/'endpoint/assay/result.json');include(folder/'endpoint/assay/result.json')
        if (e['protocol_sha256']!=digest(folder/'endpoint/protocol.json') or not e['numerical_pass'] or
                not e['all_trials_settled'] or ep['moving_protocol_sha256']!=ph or
                ep['moving_history_sha256']!=result['history_sha256']):raise ValueError('Invalid endpoint assay')
        include(folder/'endpoint/assay/paths.npz',e['paths_sha256']);endpoints.append(e)
        times=np.array([r['elapsed'] for r in h]);spread=np.array([r['log_activator_sd'] for r in h])
        growth=np.array([r['uniform_growth_max'] for r in h]);axes=np.array([r['axis_ratio'] for r in h])
        for row in h:
            spec=frozen_spectrum(row['delta'],row['volumes'],2.,.02,27.5)
            error=max(abs(spec['uniform_jacobian_max_real']-row['uniform_growth_max']),
                      float(abs(np.array(spec['lambdas'])-row['laplacian_lambdas']).max()))
            if error>1e-12 or spec['unstable_modes']!=row['unstable_modes']:raise ValueError('Changed spectrum diagnostic')
            spectrum_error=max(spectrum_error,error);spectrum_observations+=1
            spread_error=abs(float(np.std(np.log(np.array(row['chemistry'])[0])))-row['log_activator_sd'])
            if spread_error>1e-12:raise ValueError('Changed contrast diagnostic')
            contrast_error=max(contrast_error,spread_error)
        late=spread[times>=p['duration']-p['late_window']-1e-9];persistent=bool(late.min()>p['criteria']['late_log_sd_min'])
        if persistent!=result['persistent_contrast'] or float(late.min())!=result['late_min_log_activator_sd']:
            raise ValueError('Changed persistence classification')
        first,last=contact_metrics(h[0]),contact_metrics(h[-1])
        crossing=first_crossing(times,growth)
        nonpositive=np.flatnonzero(growth<=0)
        index=int(nonpositive[0]) if len(nonpositive) else None
        at_crossing=None if index is None else dict(elapsed=float(times[index]),**contact_metrics(h[index]))
        row=dict(job=job,persistent_contrast=persistent,initial_log_sd=float(spread[0]),final_log_sd=float(spread[-1]),
            late_min_log_sd=float(late.min()),peak_log_sd=float(spread.max()),
            onset_elapsed=first_crossing(times,spread,p['criteria']['late_log_sd_min'],'up'),
            uniform_growth_zero_crossing_elapsed=crossing,initial_uniform_growth=float(growth[0]),
            final_uniform_growth=float(growth[-1]),initial_lambda_max=h[0]['laplacian_lambdas'][-1],
            final_lambda_max=h[-1]['laplacian_lambdas'][-1],later_positive_growth=bool(np.any(growth[times>=150]>0)),
            initial_axis_ratio=float(axes[0]),final_axis_ratio=float(axes[-1]),max_relative_axis_change=float(abs(axes/axes[0]-1).max()),
            contact_initial=first,contact_final=last,contact_at_first_nonpositive_growth=at_crossing,
            sampled_positive_growth_intervals=exceeded_intervals(times,growth,0.),
            conductance_fraction_change=last['conductance_sum']/first['conductance_sum']-1,
            area_proxy_fraction_change=last['area_proxy_sum']/first['area_proxy_sum']-1,
            mean_polarity_initial=float(np.linalg.norm(h[0]['polarity'],axis=1).mean()),
            mean_polarity_final=float(np.linalg.norm(h[-1]['polarity'],axis=1).mean()),
            audit=result['audit'],endpoint_phase=e['phase'],endpoint_local_bistability=e['local_bistability_supported'],
            endpoint_uniform_stable=e['uniform_stable'])
        rows.append(row);histories[job['key']]=h
        print(f"verified {job['key']}: persistence={persistent} SD={spread[-1]:.6g}",flush=True)
    refinements={};diagnostics=[]
    for horizon,name in ((p['pilot_duration'],'pilot-refinement.json'),(p['duration'],'long-refinement.json')):
        report=read(root/name);include(root/name)
        if report['protocol_sha256']!=ph:raise ValueError('Changed refinement protocol')
        for path,sha in report['source_histories_sha256'].items():include(path,sha)
        recomputed=[]
        for chi in p['contrasts']:
            pair=[next(j for j in p['jobs'] if j['seed']==p['refinement_seed'] and j['family']=='uniform' and j['point']['chi']==chi and j['level']==level) for level in ('coarse','fine')]
            part=[[r for r in histories[j['key']] if r['elapsed']<=horizon+1e-9] for j in pair]
            recomputed.append(dict(chi=chi,**refinement_comparison(*part,p,horizon)))
            if horizon==p['duration']:
                a,b=[np.array([r['chemistry'] for r in trace]) for trace in part]
                error=abs(np.log(a/b));q=np.unravel_index(error.argmax(),error.shape);time,species,cell=q
                t=np.array([r['elapsed'] for r in part[0]]);maximum=error.max(axis=(1,2))
                diagnostics.append(dict(chi=chi,maximum_chemical_log_error=float(error[q]),elapsed_at_maximum=float(t[time]),
                    species=('activator','inhibitor')[species],cell_id=part[0][time]['ids'][cell],
                    coarse_value=float(a[q]),fine_value=float(b[q]),concentration_ratio_discrepancy=float(np.expm1(error[q])),
                    final_chemical_log_max=float(maximum[-1]),late_chemical_log_max=float(maximum[t>=p['duration']-p['late_window']].max()),
                    observations_exceeding_chemical_tolerance=int(np.sum(maximum>p['refinement_criteria']['chemical_log_max'])),
                    exceeded_intervals=exceeded_intervals(t,maximum,p['refinement_criteria']['chemical_log_max']),
                    times=t.tolist(),max_chemical_log_error_series=maximum.tolist()))
        if recomputed!=report['comparisons'] or report['passed']!=all(r['passed'] for r in recomputed):
            raise ValueError('Changed numerical refinement decision')
        refinements[name]=report
    coarse=[r for r in rows if r['job']['level']=='coarse'];groups=[]
    for chi in p['contrasts']:
        group=[r for r in coarse if r['job']['point']['chi']==chi]
        uniform=[r for r in group if r['job']['family']=='uniform'];pattern=[r for r in group if r['job']['family']=='pattern']
        groups.append(dict(chi=chi,initiation_histories=[r['job']['seed'] for r in uniform if r['persistent_contrast']],
            maintenance_histories=[r['job']['seed'] for r in pattern if r['persistent_contrast']],
            developed_final_sd_range=[min(r['final_log_sd'] for r in pattern),max(r['final_log_sd'] for r in pattern)],
            developed_endpoint_bistability_histories=[r['job']['seed'] for r in pattern if r['endpoint_local_bistability']]))
    trials=[t for e in endpoints for t in e['trials']]
    quality=dict(moving_quality_pass=True,native_gpu_context_pass_count=len(prefixes),endpoint_numerical_pass_count=len(endpoints),
        max_volume_error=max(r['audit']['max_volume_error'] for r in rows),min_radius=min(r['audit']['min_radius'] for r in rows),
        max_clipping=max(r['audit']['max_clipping'] for r in rows),boundary_max=max(r['audit']['boundary_max'] for r in rows),
        dilution_error_max=max(r['audit']['dilution_error_max'] for r in rows),
        context_error_max={k:max(v['errors'][k] for v in prefixes) for k in prefixes[0]['errors']},
        context_phi_abs_max=max(v['phi_abs_max'] for v in prefixes),endpoint_solver_log_max=max(t['solver_log_error'] for t in trials),
        endpoint_rhs_max=max(t['rhs_max'] for t in trials),endpoint_trial_count=len(trials),
        spectrum_observations_recomputed=spectrum_observations,spectrum_recompute_max=spectrum_error,contrast_recompute_max=contrast_error)
    result=dict(assessed_utc=datetime.now(timezone.utc).isoformat(),protocol_sha256=ph,execution_state=status['state'],
        independent_histories=3,completed_moving_jobs=len(rows),groups=groups,rows=rows,frozen_references=refs['results'],
        pilot_refinement_pass=refinements['pilot-refinement.json']['passed'],long_refinement_pass=refinements['long-refinement.json']['passed'],
        refinement_diagnostics=diagnostics,refinement_reports=refinements,quality=quality,
        interpretation='Fixed-ratio directional-tension ablation permits initiation only in history 9, concordant at both tested timesteps. Nonzero contrasts suppress sampled initiation across three histories. Developed patterns persist at all contrasts. A transient raw chemical error fails the declared full-horizon gate for history-9 chi=0; retain unresolved status without changing tolerance.',
        scope='Mature t=150-to-390 continuations in three shared histories. No cell-identity, new-zygote, autonomy, inheritance, or continuum-convergence claim. Frozen spectra are conditional diagnostics; uniform-derived endpoint assays do not exclude every patterned attractor.',
        source_sha256={**p['source_sha256'],str(Path(__file__).resolve()):digest(__file__)},evidence_sha256=evidence)
    write_json(root/'assessment.json',result)
    return result,histories


def plot(root,result,histories,destination):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    root=Path(root);destination=Path(destination);destination.parent.mkdir(parents=True,exist_ok=True)
    colors={0.:'#2878b5',.35:'#e18d25',.7:'#b84349'}
    fig,axes=plt.subplots(2,3,figsize=(14,8.2),constrained_layout=True)
    for column,seed in enumerate((7,8,9)):
        ax=axes[0,column]
        for chi,color in colors.items():
            row=next(r for r in result['rows'] if r['job']['seed']==seed and r['job']['level']=='coarse' and r['job']['family']=='uniform' and r['job']['point']['chi']==chi)
            h=histories[row['job']['key']];ax.plot([r['elapsed'] for r in h],[r['log_activator_sd'] for r in h],color=color,label=rf'Moving $\chi={chi:g}$')
        with np.load(root/f'frozen/seed-{seed}_uniform.npz') as z:
            ax.plot(z['times'],np.std(np.log(z['trajectory'][:,0]),axis=1),'k--',lw=1.5,label='Frozen reference')
        ax.axhline(.1,color='.4',ls=':',label='Contrast threshold');ax.set_yscale('log');ax.set_ylim(2e-9,3.)
        ax.set_title(f'History {seed}: near-uniform start');ax.set_xlabel('Elapsed model time');ax.set_ylabel('Across-cell SD of log activator');ax.grid(alpha=.15)
        if column==0:ax.legend(fontsize=8,loc='lower left')
    ax=axes[1,0]
    for seed,marker in ((7,'o'),(8,'s'),(9,'^')):
        rows=sorted([r for r in result['rows'] if r['job']['seed']==seed and r['job']['level']=='coarse' and r['job']['family']=='uniform'],key=lambda r:r['job']['point']['chi'])
        ax.plot([r['job']['point']['chi'] for r in rows],[r['uniform_growth_zero_crossing_elapsed'] for r in rows],marker=marker,label=f'History {seed}')
    onset=next(r['onset_elapsed'] for r in result['rows'] if r['job']['seed']==9 and r['job']['level']=='coarse' and r['job']['family']=='uniform' and r['job']['point']['chi']==0.)
    ax.scatter([0.],[onset],marker='*',s=160,color='black',label='History 9 contrast onset',zorder=5)
    ax.set_xticks([0,.35,.7]);ax.set_xlabel('Directional polarity-tension contrast');ax.set_ylabel('Elapsed time');ax.set_title('First frozen growth zero crossing');ax.legend(fontsize=8);ax.grid(alpha=.15)
    ax=axes[1,1]
    for seed,marker in ((7,'o'),(8,'s'),(9,'^')):
        rows=sorted([r for r in result['rows'] if r['job']['seed']==seed and r['job']['level']=='coarse' and r['job']['family']=='pattern'],key=lambda r:r['job']['point']['chi'])
        ax.plot([r['job']['point']['chi'] for r in rows],[r['late_min_log_sd'] for r in rows],marker=marker,label=f'History {seed}')
    ax.axhline(.1,color='.4',ls=':');ax.set_ylim(0,1.85);ax.set_xticks([0,.35,.7]);ax.set_xlabel('Directional polarity-tension contrast');ax.set_ylabel('Minimum SD over final 24 units');ax.set_title('All developed starts retain contrast');ax.legend(fontsize=8);ax.grid(alpha=.15)
    ax=axes[1,2]
    for d in result['refinement_diagnostics']:
        ax.plot(d['times'],np.maximum(d['max_chemical_log_error_series'],1e-12),color=colors[d['chi']],label=rf'$\chi={d["chi"]:g}$')
    ax.axhline(.01,color='black',ls='--',label='Declared tolerance');ax.set_yscale('log');ax.set_ylim(1e-10,.08)
    ax.set_xlabel('Elapsed model time');ax.set_ylabel('Maximum raw chemical log error');ax.set_title('History 9: coarse/fine agreement');ax.legend(fontsize=8);ax.grid(alpha=.15)
    fig.suptitle('Initiation depends on polarity and history; maintenance persists\nThree histories, paired starts; one transient quantitative refinement failure retained',fontsize=14)
    fig.savefig(destination,dpi=180);plt.close(fig)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path('outputs/polarity-robustness'))
    parser.add_argument('--figure',type=Path,default=Path('docs/images/polarity-robustness-completed.png'))
    args=parser.parse_args();result,histories=derive(args.output);plot(args.output,result,histories,args.figure)
    print(dict(completed=result['completed_moving_jobs'],histories=result['independent_histories'],
               pilot_pass=result['pilot_refinement_pass'],long_refinement_pass=result['long_refinement_pass'],quality=result['quality']))
