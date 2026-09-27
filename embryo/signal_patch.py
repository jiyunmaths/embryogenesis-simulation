"""Forced chemical-patch response of existing mechanics, followed by release.

An external regulator clamp bypasses graph reaction/transport during forcing.
It does not demonstrate spontaneous symmetry breaking. No shape is prescribed.
"""
import argparse
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from .model import Simulation
from .shape import observe, summarize, axial_angle, plot_final_shapes

CASES=('patch_default','uniform_default','patch_strong','uniform_strong','patch_no_feedback')


def prescribed_patch(centers, weights, axis=(0.,0.,1.), amplitude=.8, width=.35):
    """Smooth material patch with target-volume-weighted mean activator one.

    Values are assigned once using initial cell centers, then follow cell IDs.
    The imposed laboratory direction is recorded, never called emergent.
    """
    centers=np.asarray(centers,dtype=float);weights=np.asarray(weights,dtype=float)
    axis=np.asarray(axis,dtype=float)
    if (centers.ndim!=2 or centers.shape[1]!=3 or weights.shape!=(len(centers),)
            or not np.isfinite(centers).all() or not np.isfinite(weights).all() or np.any(weights<=0)
            or axis.shape!=(3,) or not np.isfinite(axis).all() or np.linalg.norm(axis)==0
            or not np.isfinite(amplitude) or not 0<amplitude<1 or not np.isfinite(width) or width<=0):
        raise ValueError('invalid patch geometry, weights, axis, amplitude, or width')
    axis=axis/np.linalg.norm(axis)
    center=weights@centers/weights.sum()
    score=np.tanh(((centers-center)@axis)/width)
    score-=weights@score/weights.sum()
    spread=np.max(np.abs(score))
    if spread<1e-12:
        raise ValueError('cell centers do not span the imposed patch direction')
    score/=spread
    return 1+amplitude*score, np.ones(len(centers))


def patch_observe(sim, axis, start, release):
    row=observe(sim)
    covariance=np.asarray(row['covariance'])
    parallel=float(axis@covariance@axis)
    perpendicular=float((np.trace(covariance)-parallel)/2)
    row.update({'elapsed_from_intervention':sim.time-start,
                'phase':'clamp' if sim.time<=release+1e-10 else 'released',
                'angle_to_imposed_axis_degrees':axial_angle(row['principal_axis'],axis),
                'imposed_axis_ratio':float(np.sqrt(parallel/max(perpendicular,1e-30))),
                'mean_polarity_norm':float(np.linalg.norm(sim.polarity,axis=1).mean()),
                'mean_nominal_directional_tension_amplitude':float(sim.config.polarity_tension*np.linalg.norm(sim.polarity,axis=1).mean()),
                'mean_directional_tension_amplitude':float(sim.config.polarity_tension*np.linalg.norm(sim.polarity,axis=1).mean()) if sim.config.feedback and sim.config.polarity_enabled else 0.})
    return row


def run(output, checkpoint, case='patch_default', clamp_duration=30., release_duration=30.,
        sample_interval=.6, amplitude=.8, width=.35, axis=(0.,0.,1.), strong_tension=.75):
    output=Path(output);checkpoint=Path(checkpoint)
    if output.exists():raise FileExistsError('choose a fresh output directory')
    if case not in CASES:raise ValueError('unknown patch case')
    if any(not np.isfinite(x) or x<=0 for x in (clamp_duration,release_duration,sample_interval)):
        raise ValueError('positive finite clamp, release, and sample durations required')
    if sample_interval>min(clamp_duration,release_duration)/2:
        raise ValueError('record at least three samples per phase')
    sim=Simulation.restore(checkpoint)
    if sim.divisions or len(sim.phi)!=sim.config.max_cells:
        raise ValueError('patch experiment requires completed cleavage at the cell cap')
    if (not sim.config.signaling or not sim.config.differentiation or not sim.config.polarity_enabled
            or sim.config.fate_noise or sim.config.exposure_bias or sim.config.neighbor_inhibition):
        raise ValueError('patch protocol requires signaling, fate and polarity with no legacy bias or fate noise')
    original_config=asdict(sim.config)
    a,h=prescribed_patch(sim.centers(),sim.target,axis,amplitude,width)
    axis=np.array(axis,dtype=float,copy=True);axis/=np.linalg.norm(axis)
    if case.startswith('uniform'):a=np.ones_like(a)
    config={'feedback':case!='patch_no_feedback'}
    if case.endswith('strong'):config['fate_tension']=strong_tension
    sim.config=replace(sim.config,**config);sim.config.validate()
    dt=sim.config.dt
    counts=[round(x/dt) for x in (clamp_duration,release_duration,sample_interval)]
    if any(n<1 or not np.isclose(n*dt,x) for n,x in zip(counts,(clamp_duration,release_duration,sample_interval))):
        raise ValueError('durations and sample interval must be multiples of checkpoint dt')
    forcing,relaxing,every=counts;steps=forcing+relaxing
    if forcing%every or relaxing%every:raise ValueError('sample interval must divide both phase durations')
    start=sim.time;release=start+clamp_duration
    sim.config=replace(sim.config,steps=sim.step_number+steps,save_every=every)
    unforced_initial=observe(sim)
    original_signals={'activator':sim.activator.tolist(),'inhibitor':sim.inhibitor.tolist()}
    sim.activator=a.copy();sim.inhibitor=h.copy()
    report={'schema':1,'mode':case,'parameters':asdict(sim.config),
            'provenance':{'origin':'mature_checkpoint','path':str(checkpoint),
                          'sha256':hashlib.sha256(checkpoint.read_bytes()).hexdigest(),
                          'time':start,'original_config':original_config},
            'scope':'Forced-response diagnostic, not spontaneous symmetry breaking. No prescribed geometry or new force law.',
            'protocol':{'clamp_duration':clamp_duration,'release_duration':release_duration,
                        'release_time':release,'sample_interval':sample_interval,'axis':axis.tolist(),
                        'amplitude':amplitude,'width':width,'strong_tension':strong_tension,
                        'assignment':'fixed material cell IDs from initial centers; no divisions',
                        'regulators':'external activator/inhibitor clamp, then original GM reaction and graph transport',
                        'uniform_control':'activator=1, inhibitor=1; same target-weighted means',
                        'graph_spectra':'autonomous GM reference only; GM evolution is bypassed during clamp',
                        'initial_cell_ids':sim.ids.tolist(),'activator_by_id':a.tolist(),'inhibitor_by_id':h.tolist()},
            'unforced_initial':unforced_initial,'original_signals':original_signals}
    output.mkdir(parents=True)
    (output/'config.json').write_text(json.dumps(asdict(sim.config),indent=2)+'\n')
    (output/'protocol.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    history=[];frames=[];started=time.monotonic();ids=sim.ids.copy()
    for step in range(steps+1):
        if step%every==0 or step==steps:
            row=patch_observe(sim,axis,start,release);history.append(row)
            frames.append({'metrics':row,'cells':sim.surfaces(max_points=450),'graph':sim.graph_snapshot()})
            if step%(10*every)==0 or step==steps:
                (output/'history.json').write_text(json.dumps(history,indent=2,allow_nan=False)+'\n')
                print(f'{case}: t={sim.time:.2f}, {row["phase"]}, R={row["axis_ratio"]:.4f}, imposed-axis ratio={row["imposed_axis_ratio"]:.4f}, elapsed={time.monotonic()-started:.1f}s',flush=True)
        if step<steps:
            sim.step(prescribed_signals=(a,h) if step<forcing else None)
            if not np.array_equal(ids,sim.ids):raise RuntimeError('cell IDs changed during material patch test')
    last=max((r['division'] for r in sim.lineage if r['division'] is not None),default=None)
    forced=[r for r in history if r['time']<=release+1e-10]
    report.update({'elapsed_seconds':time.monotonic()-started,'history':history,
                   'summary':summarize(history,min(15.,release_duration),last),
                   'clamp_summary':summarize(forced,min(15.,clamp_duration),last),
                   'whole_run_numerical_checks':summarize(history,clamp_duration+release_duration,last)['numerical_checks']})
    (output/'analysis.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    payload=json.dumps({'config':asdict(sim.config),'frames':frames},separators=(',',':'),allow_nan=False)
    (output/'trajectory.json').write_text(payload)
    from .surface import viewer_template
    template=viewer_template()
    template=template.replace('One cell. Two possible identities.',f'Forced signal response · {case.replace("_"," ")}.')
    template=template.replace('A freely evolving 3D aggregate of deformable cells.',
                              f'Imposed chemical pattern until t={release:g}, then autonomous signaling. This is not spontaneous symmetry breaking.')
    (output/'viewer.html').write_text(template.replace('__SIMULATION_DATA__',payload))
    sim.checkpoint(output/'final_state.npz')
    (output/'lineage.json').write_text(json.dumps(sim.lineage,indent=2)+'\n')
    return report


def compare(directories, output):
    output=Path(output)
    if output.exists():raise FileExistsError('choose a fresh comparison directory')
    reports=[json.loads((Path(p)/'analysis.json').read_text()) for p in directories]
    by_case={r['mode']:r for r in reports}
    if len(by_case)!=len(reports) or set(by_case)!=set(CASES):
        raise ValueError('provide exactly one run for each of the five patch cases')
    baseline=reports[0]
    common=('clamp_duration','release_duration','release_time','sample_interval','axis','amplitude','width','strong_tension','initial_cell_ids')
    for r in reports:
        if (r['provenance']['sha256']!=baseline['provenance']['sha256']
                or any(r['protocol'][k]!=baseline['protocol'][k] for k in common)
                or any(r['parameters'][k]!=baseline['parameters'][k] for k in baseline['parameters'] if k not in ('feedback','fate_tension'))):
            raise ValueError('patch comparisons require matched checkpoint and non-intervention parameters')
    pairs=[]
    for forced,control in [('patch_default','uniform_default'),('patch_strong','uniform_strong'),('patch_default','patch_no_feedback')]:
        a=by_case[forced];b=by_case[control]
        if a['parameters']['fate_tension']!=b['parameters']['fate_tension']:
            raise ValueError('paired patch and control must match fate-tension strength')
        if [r['time'] for r in a['history']]!=[r['time'] for r in b['history']]:
            raise ValueError('comparison samples must match')
        endpoints={}
        for name,t in [('end_of_clamp',a['protocol']['release_time']),('end_of_release',a['history'][-1]['time'])]:
            left=next(r for r in a['history'] if np.isclose(r['time'],t))
            right=next(r for r in b['history'] if np.isclose(r['time'],t))
            ca=np.asarray(left['covariance']);cb=np.asarray(right['covariance'])
            endpoints[name]={'time':t,'axis_ratio_difference':left['axis_ratio']-right['axis_ratio'],
                             'imposed_axis_ratio_difference':left['imposed_axis_ratio']-right['imposed_axis_ratio'],
                             'covariance_relative_difference':float(np.linalg.norm(ca-cb)/np.linalg.norm(cb)),
                             'axis_angle_difference_degrees':axial_angle(left['principal_axis'],right['principal_axis'])}
        pairs.append({'case':forced,'control':control,'endpoints':endpoints,
                      'both_numerical_screens_pass':all(a['whole_run_numerical_checks'].values()) and all(b['whole_run_numerical_checks'].values())})
    result={'schema':1,'comparisons':pairs,
            'runs':[{'directory':str(p),'case':r['mode'],'summary':r['summary'],
                     'clamp_summary':r['clamp_summary'],'whole_run_numerical_checks':r['whole_run_numerical_checks']} for p,r in zip(directories,reports)],
            'interpretation':'Externally forced chemical pattern followed by release. Differences measure response of existing mechanics, not spontaneous axis selection or validated morphogenesis.'}
    output.mkdir(parents=True)
    (output/'analysis.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    plot(reports,output)
    plot_final_shapes(directories,output,title='After release · shared camera, axes and fate colors')
    clamp_index=round(baseline['protocol']['clamp_duration']/baseline['protocol']['sample_interval'])
    plot_final_shapes(directories,output,frame_index=clamp_index,filename='clamp_shapes.png',
                      title='End of imposed chemical clamp · shared camera, axes and fate colors')
    return result


def plot(reports,output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(3,2,figsize=(11,10),constrained_layout=True)
    for r in reports:
        h=r['history'];t=[m['time'] for m in h];name=r['mode'].replace('_',' ')
        for ax,key,scale in [(axes[0,0],'axis_ratio',1),(axes[0,1],'imposed_axis_ratio',1),
                             (axes[1,0],'activator_std',1),(axes[1,1],'max_cell_volume_error',100),
                             (axes[2,0],'boundary_occupancy',1),(axes[2,1],'min_radius_grid_cells',1)]:
            ax.plot(t,[scale*m[key] for m in h],label=name)
    for ax,label in zip(axes.ravel(),['Shape axis ratio','Ratio along imposed axis / transverse','Activator standard deviation','Largest cell volume error (%)','Boundary occupancy','Smallest radius / grid spacing']):
        ax.set(xlabel='Simulation time',ylabel=label)
        ax.axvline(reports[0]['protocol']['release_time'],color='gray',ls=':',label='Release')
    axes[2,0].axhline(.01,color='gray',ls='--',lw=1)
    axes[2,1].axhline(4.,color='gray',ls='--',lw=1)
    axes[0,0].legend(fontsize=7)
    fig.suptitle('Forced chemical-patch response → autonomous release · no prescribed shape')
    fig.savefig(output/'response.png',dpi=170);plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    advance=sub.add_parser('run')
    advance.add_argument('--output',type=Path,required=True);advance.add_argument('--checkpoint',type=Path,required=True)
    advance.add_argument('--case',choices=CASES,default='patch_default')
    advance.add_argument('--axis',nargs=3,type=float,default=(0.,0.,1.))
    for name,value in [('clamp-duration',30.),('release-duration',30.),('sample-interval',.6),('amplitude',.8),('width',.35),('strong-tension',.75)]:
        advance.add_argument('--'+name,type=float,default=value)
    analyze=sub.add_parser('compare');analyze.add_argument('--runs',nargs='+',required=True);analyze.add_argument('--output',type=Path,required=True)
    args=vars(parser.parse_args());command=args.pop('command')
    if command=='run':run(**args)
    else:compare(args.pop('runs'),args['output'])


if __name__=='__main__':main()
