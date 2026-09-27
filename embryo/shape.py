"""Finite-window persistence and causal controls for coupled embryo shape.

This is a pilot on the existing deformable-cell/contact-graph model, not a
continuum-limit or ensemble claim. Geometry and axes are never prescribed.
"""

import argparse
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.ndimage import label
from scipy.sparse.csgraph import connected_components

from .model import Config, Simulation, occupancy
from .signaling import normalized_graph
from .graph_analysis import cycle
from .signaling import stability


MODES = ('full', 'no_feedback', 'no_polarity_tension')
CRITERIA = {'axis_ratio_min': 1.2, 'axis_gap_min': .05, 'axis_rotation_max_degrees': 15.,
            'late_ratio_relative_range_max': .05, 'feedback_axis_ratio_excess_min': .05,
            'cell_volume_error_max': .05, 'boundary_occupancy_max': .01,
            'clipped_fraction_max': .01, 'radius_grid_cells_min': 4.,
            'activator_std_min': .1}


def axial_angle(a, b):
    """An unoriented axis: a and -a have zero angular difference."""
    if a is None or b is None:
        return None
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    if a.shape != (3,) or b.shape != (3,) or not np.isfinite([a,b]).all():
        raise ValueError('axes must be finite three-vectors')
    norm = np.linalg.norm(a)*np.linalg.norm(b)
    if norm <= 1e-14:
        return None
    return float(np.degrees(np.arccos(np.clip(abs(a@b)/norm, 0., 1.))))


def shape_tensor(weights, coordinates, gap_min=.05):
    """Diffuse-union inertia, with no arbitrary axis at a degenerate spectrum."""
    weights=np.asarray(weights, dtype=float).ravel()
    points=np.asarray(coordinates, dtype=float).reshape(3,-1)
    if (len(weights)!=points.shape[1] or not np.isfinite(weights).all()
            or not np.isfinite(points).all() or np.any(weights<0) or weights.sum()<=0):
        raise ValueError('nonnegative nonempty weights must match finite coordinates')
    center=points@weights/weights.sum()
    offset=points-center[:,None]
    covariance=(offset*weights)@offset.T/weights.sum()
    eigen,vectors=np.linalg.eigh(covariance)
    gap=float((eigen[-1]-eigen[-2])/max(eigen[-1],1e-30))
    axis=vectors[:,-1].tolist() if gap>=gap_min else None
    return {'centroid':center.tolist(), 'covariance':covariance.tolist(),
            'eigenvalues':eigen.tolist(), 'principal_axis':axis, 'long_axis_gap':gap,
            'tensor_axis_ratio':float(np.sqrt(eigen[-1]/max(eigen[0],1e-30)))}


def observe(simulation):
    sim=simulation
    h=occupancy(sim.phi); union=np.minimum(h.sum(axis=0),1)
    row=sim.metrics()
    row.update(shape_tensor(union,sim.xyz,CRITERIA['axis_gap_min']))
    first_axis=sim.lineage[0].get('division_axis')
    row['angle_to_first_cleavage_degrees']=axial_angle(row['principal_axis'],first_axis)
    graph=sim.signaling_graph()
    row['contact_components']=int(connected_components(graph.weights,directed=False,return_labels=False))
    row['max_cell_components_phi_0_5']=int(max(label(p>.5)[1] for p in sim.phi))
    row['union_components_occupancy_0_1']=int(label(union>.1)[1])
    row['union_components_occupancy_0_5']=int(label(union>.5)[1])
    volumes=sim.volumes();centers=sim.centers()
    signal_mean=float(volumes@sim.activator/volumes.sum())
    dipole=np.sum((volumes*(sim.activator-signal_mean))[:,None]*(centers-np.array(row['centroid'])),axis=0)/volumes.sum()
    row['signal_dipole_magnitude']=float(np.linalg.norm(dipole))
    row['signal_shape_angle_degrees']=axial_angle(row['principal_axis'],dipole)
    return row


def summarize(history, late_duration, last_cleavage):
    end=history[-1]['time'];start=end-late_duration
    late=[row for row in history if row['time']>=start-1e-10]
    complete=history[0]['time']<=start+1e-10 and len(late)>=3
    reference=late[0]['principal_axis']
    angles=[axial_angle(reference,r['principal_axis']) for r in late]
    all_identifiable=all(r['principal_axis'] is not None for r in late)
    ratios=np.array([r['axis_ratio'] for r in late])
    ratio_range=float(np.ptp(ratios)/ratios.mean())
    rotation=max((a for a in angles if a is not None),default=None)
    checks={
        'late_window_recorded':bool(complete),
        'cleavage_finished_before_late_window':bool(last_cleavage is not None and last_cleavage<start and all(r['dividing_cells']==0 for r in late)),
        'persistent_anisotropy':bool(np.min(ratios)>=CRITERIA['axis_ratio_min']),
        'identifiable_long_axis':all_identifiable,
        'axis_memory':bool(all_identifiable and rotation is not None and rotation<=CRITERIA['axis_rotation_max_degrees']),
        'late_shape_ratio_stable':bool(ratio_range<=CRITERIA['late_ratio_relative_range_max']),
        'late_signaling_contrast':bool(min(r['activator_std'] for r in late)>=CRITERIA['activator_std_min']),
    }
    quality={
        'cell_volumes_within_5_percent':all(r['max_cell_volume_error']<=CRITERIA['cell_volume_error_max'] for r in late),
        'clear_of_domain_boundary':all(r['boundary_occupancy']<=CRITERIA['boundary_occupancy_max'] for r in late),
        'little_phase_field_clipping':all(r['clipped_fraction']<=CRITERIA['clipped_fraction_max'] for r in late),
        'minimum_radius_resolved_by_4_grid_cells':all(r['min_radius_grid_cells']>=CRITERIA['radius_grid_cells_min'] for r in late),
        'each_cell_connected_at_phi_0_5':all(r['max_cell_components_phi_0_5']==1 for r in late),
        'contact_graph_connected':all(r['contact_components']==1 for r in late),
    }
    return {'late_window':[start,end], 'late_samples':len(late), 'minimum_late_axis_ratio':float(ratios.min()),
            'mean_late_axis_ratio':float(ratios.mean()),'late_ratio_relative_range':ratio_range,
            'maximum_late_axis_rotation_degrees':rotation,
            'final_angle_to_first_cleavage_degrees':history[-1]['angle_to_first_cleavage_degrees'],
            'checks':checks,'numerical_checks':quality,
            'finite_window_shape_persistence':all(value for key,value in checks.items() if key!='late_signaling_contrast'),
            'numerical_screen_passed':all(quality.values()),
            'interpretation':'Finite-window single-trajectory screen; grid/time/domain convergence, rotation tests, and seed ensembles remain required.'}


def apply_mode(sim, mode):
    if mode not in MODES: raise ValueError('unknown shape control')
    changes={'feedback':True}
    if mode=='no_feedback': changes['feedback']=False
    if mode=='no_polarity_tension': changes['polarity_tension']=0.
    sim.config=replace(sim.config,**changes)
    sim.config.validate()


def run(output, mode='full', checkpoint=None, duration=30., late_duration=15.,
        seed=7, grid=40, dt=.015, sample_interval=.6):
    output=Path(output)
    if output.exists(): raise FileExistsError('choose a fresh output directory')
    for name,value in [('duration',duration),('late_duration',late_duration),('sample_interval',sample_interval),('dt',dt)]:
        if isinstance(value,bool) or not np.isfinite(value) or value<=0: raise ValueError(f'{name} must be positive and finite')
    if late_duration>duration or sample_interval>late_duration/2: raise ValueError('record at least three late-window samples')
    if checkpoint is None:
        sim=Simulation(Config(seed=seed,grid=grid,dt=dt))
        provenance={'origin':'zygote','seed':seed}
    else:
        checkpoint=Path(checkpoint)
        sim=Simulation.restore(checkpoint)
        if sim.divisions or len(sim.phi)!=sim.config.max_cells:
            raise ValueError('branch checkpoint must have finished all cleavages at the cell cap')
        provenance={'origin':'mature_checkpoint','path':str(checkpoint),'sha256':hashlib.sha256(checkpoint.read_bytes()).hexdigest(),
                    'time':sim.time,'seed':sim.config.seed,'original_config':asdict(sim.config)}
    # Branching preserves the original clock and step; no interpolation or
    # phase-field resampling is used. Grid/dt flags apply only to zygote starts.
    actual_dt=sim.config.dt
    steps=round(duration/actual_dt);save_every=round(sample_interval/actual_dt)
    if steps<1 or save_every<1 or not np.isclose(steps*actual_dt,duration) or not np.isclose(save_every*actual_dt,sample_interval):
        raise ValueError('duration and sample_interval must be integer multiples of the simulation dt')
    if sim.config.fate_noise or sim.config.partition_noise or sim.config.exposure_bias or sim.config.neighbor_inhibition:
        raise ValueError('shape screen requires legacy fate noise and independent geometry/fate biases to be zero')
    apply_mode(sim,mode)
    sim.config=replace(sim.config,steps=sim.step_number+steps,save_every=save_every)
    output.mkdir(parents=True)
    report={'schema':1,'mode':mode,'provenance':provenance,'parameters':asdict(sim.config),
            'duration':duration,'late_duration':late_duration,'criteria':CRITERIA,
            'scope':f'3D diffuse-cell mechanics with {sim.config.signal_transport} GM signaling. No prescribed shape or imposed chemical axis.',
            'preflight':{'actual_initial_graph':sim.graph_snapshot(),
                         'operator':sim.config.signal_transport}}
    (output/'config.json').write_text(json.dumps(asdict(sim.config),indent=2)+'\n')
    (output/'protocol.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    history=[];frames=[];started=time.monotonic()
    for step in range(steps+1):
        if step%save_every==0 or step==steps:
            row=observe(sim);history.append(row)
            frames.append({'metrics':row,'cells':sim.surfaces(max_points=450),'graph':sim.graph_snapshot()})
            if step%(save_every*10)==0 or step==steps:
                print(f'{output.name}: t={sim.time:.2f}, cells={len(sim.phi)}, ratio={row["axis_ratio"]:.4f}, sigma(a)={row["activator_std"]:.3g}, elapsed={time.monotonic()-started:.1f}s',flush=True)
                (output/'history.json').write_text(json.dumps(history,indent=2,allow_nan=False)+'\n')
        if step<steps:sim.step()
    last=max((r['division'] for r in sim.lineage if r['division'] is not None),default=None)
    report.update({'last_cleavage':last,'summary':summarize(history,late_duration,last),'elapsed_seconds':time.monotonic()-started,'history':history})
    (output/'analysis.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    payload=json.dumps({'config':asdict(sim.config),'frames':frames},separators=(',',':'),allow_nan=False)
    (output/'trajectory.json').write_text(payload)
    from .surface import viewer_template
    template=viewer_template()
    origin_label='Mature-state continuation' if checkpoint is not None else 'Development from one zygote'
    template=template.replace('One cell. Two possible identities.',f'Shape persistence · {mode.replace("_", " ")}.')
    template=template.replace('A freely evolving 3D aggregate of deformable cells.',
        f'{origin_label}. Actual simulated cell surfaces; see the comparison report for causal and numerical checks.')
    (output/'viewer.html').write_text(template.replace('__SIMULATION_DATA__',payload))
    sim.checkpoint(output/'final_state.npz')
    (output/'lineage.json').write_text(json.dumps(sim.lineage,indent=2)+'\n')
    print(json.dumps(report['summary'],indent=2),flush=True)
    return report


def extend(source, output, until=180., late_duration=30., sample_interval=.6):
    """Continue an existing branch, retaining its original intervention ancestry.

    Source files remain untouched. Joined histories/playbacks cover the old and
    new windows; continuation hashes document that t=90 controls have distinct
    states descended from their common t=60 checkpoint.
    """
    source=Path(source);output=Path(output)
    previous=json.loads((source/'analysis.json').read_text())
    checkpoint=source/'final_state.npz'
    sim=Simulation.restore(checkpoint)
    if (not np.isfinite(until) or until<=sim.time
            or not np.isclose(previous['history'][-1]['time'],sim.time)
            or previous['parameters']!=asdict(sim.config)):
        raise ValueError('extension requires a consistent completed source and a later end time')
    old_payload=json.loads((source/'trajectory.json').read_text())
    if not np.isclose(old_payload['frames'][-1]['metrics']['time'],sim.time):
        raise ValueError('source trajectory must end at the checkpoint time')
    result=run(output,mode=previous['mode'],checkpoint=checkpoint,
               duration=until-sim.time,late_duration=late_duration,sample_interval=sample_interval)
    continuation={'source':str(source),'checkpoint_sha256':result['provenance']['sha256'],
                  'source_analysis_sha256':hashlib.sha256((source/'analysis.json').read_bytes()).hexdigest(),
                  'start_time':sim.time,'end_time':until,'elapsed_seconds':result['elapsed_seconds']}
    result['continuations']=previous.get('continuations',[])+[continuation]
    result['provenance']=previous['provenance']
    result['history']=previous['history']+result['history'][1:]
    result['duration']=until-result['history'][0]['time']
    result['summary']=summarize(result['history'],late_duration,result['last_cleavage'])
    result['elapsed_seconds_scope']='most recent continuation only'
    (output/'analysis.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    (output/'history.json').write_text(json.dumps(result['history'],indent=2,allow_nan=False)+'\n')
    protocol={k:v for k,v in result.items() if k not in ('history','summary')}
    (output/'protocol.json').write_text(json.dumps(protocol,indent=2,allow_nan=False)+'\n')
    payload=json.loads((output/'trajectory.json').read_text())
    old_json=json.dumps(payload,separators=(',',':'),allow_nan=False)
    payload['frames']=old_payload['frames']+payload['frames'][1:]
    new_json=json.dumps(payload,separators=(',',':'),allow_nan=False)
    (output/'trajectory.json').write_text(new_json)
    viewer=output/'viewer.html'
    viewer.write_text(viewer.read_text().replace(old_json,new_json))
    return result


def compare(directories, output):
    output=Path(output)
    if output.exists(): raise FileExistsError('choose a fresh comparison directory')
    reports=[json.loads((Path(p)/'analysis.json').read_text()) for p in directories]
    full=[r for r in reports if r['mode']=='full' and r['provenance']['origin']=='mature_checkpoint']
    if len(full)!=1: raise ValueError('provide exactly one mature full-feedback continuation')
    baseline=full[0]
    branches=[r for r in reports if r['provenance'].get('sha256')==baseline['provenance']['sha256'] and r['mode']!='full']
    developmental=[r for r in reports if r['provenance']['origin']=='zygote' and r['mode']=='no_feedback']
    results=[]
    for control in branches+developmental:
        matched=control['provenance']['origin']=='mature_checkpoint'
        branch_ignored={'feedback','polarity_tension','steps','save_every'}
        if matched and any(control['parameters'][key]!=value for key,value in baseline['parameters'].items() if key not in branch_ignored):
            raise ValueError('matched branches must retain non-intervention parameters')
        # The only intentional difference for the developmental comparison is
        # mechanical feedback. A shared seed alone is not a matched trajectory.
        ignored={'feedback', 'steps', 'save_every'}
        if not matched and any(control['parameters'][key]!=value for key,value in baseline['parameters'].items() if key not in ignored):
            raise ValueError('developmental control must match non-intervention parameters')
        if baseline['summary']['late_window']!=control['summary']['late_window']:
            raise ValueError('matched branches must share the late observation window')
        lower=baseline['summary']['late_window'][0]
        reference={round(r['time'],9):r for r in control['history']}
        pairs=[(r,reference[round(r['time'],9)]) for r in baseline['history'] if r['time']>=lower-1e-10 and round(r['time'],9) in reference]
        full_late=[r for r in baseline['history'] if r['time']>=lower-1e-10]
        if len(pairs)<3 or len(pairs)!=len(full_late):
            raise ValueError('controls must cover every full-feedback late sample')
        excess=np.array([a['axis_ratio']-b['axis_ratio'] for a,b in pairs])
        results.append({'control':control['mode'],'origin':control['provenance']['origin'],'matched_checkpoint':matched,'late_samples':len(pairs),
                        'mean_axis_ratio_excess':float(excess.mean()),'minimum_axis_ratio_excess':float(excess.min()),
                        'excess_above_0_05_throughout_late_window':bool(np.min(excess)>=CRITERIA['feedback_axis_ratio_excess_min'])})
    result={'schema':1,'runs':[{'directory':str(p),'mode':r['mode'],'origin':r['provenance']['origin'],'summary':r['summary']} for p,r in zip(directories,reports)],
            'matched_branch_comparisons':[r for r in results if r['matched_checkpoint']],
            'developmental_comparisons':[r for r in results if not r['matched_checkpoint']],
            'interpretation':'Removing feedback from a mature state tests maintenance, not its role during development. The zygote control is a separate developmental history. No ensemble or numerical convergence claim.'}
    output.mkdir(parents=True)
    (output/'analysis.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    plot_comparison(reports,output)
    if all((Path(p)/'trajectory.json').exists() for p in directories):
        plot_final_shapes(directories,output)
    return result


def plot_comparison(reports,output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(3,2,figsize=(11,10),constrained_layout=True)
    for r in reports:
        h=r['history'];times=[m['time'] for m in h]
        name=r['mode'].replace('_',' ')+(' · from zygote' if r['provenance']['origin']=='zygote' else ' · branch')
        axes[0,0].plot(times,[m['axis_ratio'] for m in h],label=name)
        axes[0,1].plot(times,[m['angle_to_first_cleavage_degrees'] for m in h],label=name)
        axes[1,0].plot(times,[m['activator_std'] for m in h],label=name)
        axes[1,1].plot(times,[100*m['max_cell_volume_error'] for m in h],label=name)
        axes[2,0].plot(times,[m['boundary_occupancy'] for m in h],label=name)
        axes[2,1].plot(times,[m['min_radius_grid_cells'] for m in h],label=name)
    axes[0,0].axhline(CRITERIA['axis_ratio_min'],color='gray',ls=':',lw=1)
    axes[0,0].set(ylabel='Shape axis ratio',title='Anisotropy and causal controls')
    axes[0,0].legend(fontsize=7)
    axes[0,1].set(ylabel='Unoriented angle (degrees)',title='Alignment with first cleavage axis')
    axes[1,0].set(ylabel='Activator standard deviation',title='Signaling contrast')
    axes[1,1].set(ylabel='Largest cell volume error (%)',title='Numerical diagnostic')
    axes[2,0].set(ylabel='Boundary occupancy',title='Domain-size screen')
    axes[2,0].axhline(CRITERIA['boundary_occupancy_max'],color='gray',ls=':',lw=1)
    axes[2,1].set(ylabel='Smallest radius / grid spacing',title='Cell-resolution screen')
    axes[2,1].axhline(CRITERIA['radius_grid_cells_min'],color='gray',ls=':',lw=1)
    for ax in axes.ravel(): ax.set_xlabel('Simulation time')
    fig.suptitle('Coupled shape-persistence pilot · finite trajectories, no ensemble claim')
    fig.savefig(output/'persistence.png',dpi=170);plt.close(fig)


def plot_final_shapes(directories, output, *, frame_index=-1, filename="final_shapes.png",
                      title="Final 3D aggregates · matching camera and spatial limits"):
    """Same camera, axes, and fate scale across actual final surface samples."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    fig=plt.figure(figsize=(11, 4.5*((len(directories)+1)//2)), constrained_layout=True)
    norm=Normalize(-1,1); cmap=plt.get_cmap('coolwarm')
    axes=[]
    for k,directory in enumerate(directories):
        directory=Path(directory)
        payload=json.loads((directory/'trajectory.json').read_text())
        report=json.loads((directory/'analysis.json').read_text())
        frame=payload['frames'][frame_index]; row=frame['metrics']
        lineage=json.loads((directory/'lineage.json').read_text())
        ax=fig.add_subplot((len(directories)+1)//2,2,k+1,projection='3d');axes.append(ax)
        for cell in frame['cells']:
            p=np.asarray(cell['points'])
            if len(p): ax.scatter(*p.T,s=2,c=[cmap(norm(np.tanh(cell['fate'])))],alpha=.55,linewidths=0)
        center=np.asarray(row['centroid']);extent=payload['config']['extent']
        overlays=[(row['principal_axis'],'#166d5c','Measured long axis'),
                  (lineage[0].get('division_axis'),'#bb922d','First cleavage axis')]
        imposed=report.get('protocol',{}).get('axis')
        if imposed is not None:
            overlays.append((imposed,'#9c4165','Imposed chemical axis'))
        for axis,color,label_text in overlays:
            if axis is not None:
                endpoints=center[:,None]+np.asarray(axis)[:,None]*np.array([-.65,.65])[None,:]
                ax.plot(*endpoints,color=color,lw=2,label=label_text)
        ax.set(xlim=(-extent,extent),ylim=(-extent,extent),zlim=(-extent,extent))
        ax.set_box_aspect((1,1,1));ax.view_init(elev=24,azim=35);ax.set_axis_off()
        origin='from zygote' if report['provenance']['origin']=='zygote' else 'mature-state branch'
        ax.set_title(f"{report['mode'].replace('_',' ')} · {origin}\nt={row['time']:g}, axis ratio={row['axis_ratio']:.4f}",fontsize=10)
        if k==0:ax.legend(fontsize=7,loc='lower left')
    fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),ax=axes,shrink=.55,label='tanh(fate variable), shared scale')
    fig.suptitle(title)
    fig.savefig(output/filename,dpi=150);plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    advance=sub.add_parser('run');advance.add_argument('--output',type=Path,required=True)
    advance.add_argument('--mode',choices=MODES,default='full')
    advance.add_argument('--checkpoint',type=Path)
    for name,value,kind in [('duration',30.,float),('late-duration',15.,float),('seed',7,int),('grid',40,int),('dt',.015,float),('sample-interval',.6,float)]:
        advance.add_argument('--'+name,type=kind,default=value)
    continuation=sub.add_parser('extend')
    continuation.add_argument('--source',type=Path,required=True)
    continuation.add_argument('--output',type=Path,required=True)
    continuation.add_argument('--until',type=float,default=180.)
    continuation.add_argument('--late-duration',type=float,default=30.)
    continuation.add_argument('--sample-interval',type=float,default=.6)
    analyze=sub.add_parser('compare');analyze.add_argument('--runs',nargs='+',required=True);analyze.add_argument('--output',type=Path,required=True)
    args=vars(parser.parse_args());command=args.pop('command')
    try:
        if command=='run':run(**args)
        elif command=='extend':extend(**args)
        else:compare(args['runs'],args['output'])
    except (ValueError,FileExistsError) as error:parser.error(str(error))


if __name__=='__main__':main()
