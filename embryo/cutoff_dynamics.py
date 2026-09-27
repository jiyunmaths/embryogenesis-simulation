"""Matched signaling-cutoff interventions with polarity filtering held fixed."""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from .model import Simulation
from .resolution import write_json, _steps


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def prepare(checkpoint, output, until=90., interval=.6):
    checkpoint, output = Path(checkpoint).resolve(), Path(output)
    if output.exists():
        raise FileExistsError('choose a fresh output directory')
    sim = Simulation.restore(checkpoint)
    if sim.config.signal_transport != 'conservative' or sim.divisions or len(sim.phi) != sim.config.max_cells:
        raise ValueError('requires a mature conservative checkpoint with no pending or future divisions')
    _steps(until-sim.time, sim.config.dt)
    sample_steps = _steps(interval, sim.config.dt)
    if round((until-sim.time)/sim.config.dt) % sample_steps:
        raise ValueError('duration must be divisible by observation interval')
    if sim.config.graph_contact_cutoff != .02:
        raise ValueError('baseline checkpoint must use cutoff 0.02')
    effective = sim.config.graph_contact_cutoff if sim.config.polarity_contact_cutoff == -1 else sim.config.polarity_contact_cutoff
    if effective != .02:
        raise ValueError('baseline polarity cutoff must be 0.02')
    protocol = {'checkpoint':str(checkpoint), 'checkpoint_sha256':sha(checkpoint),
                'model_sha256':sha(Path(__file__).with_name('model.py')),
                'source_config':asdict(sim.config), 'start':sim.time, 'until':until,
                'sample_steps':sample_steps, 'signal_cutoffs':[.02,.01,.04], 'polarity_cutoff':.02,
                'criteria': {'signal_absolute_rms_max':.01, 'contrast_absolute_difference_max':.005,
                             'axis_ratio_relative_max':.01, 'continuous_fate_absolute_max':.05,
                             'boundary_occupancy_max':.01, 'cell_volume_error_max':.05,
                             'minimum_radius_grid_cells':4.},
                'scope':'Single-seed interventions on a mature checkpoint, fixed polarity cutoff, full downstream feedback. Sampled per-cell signal RMS uses equilibrium concentration 1 as scale; contrast uses unweighted cell standard deviation. Exact labels must match baseline at every observation. This tests identity retention, not initial differentiation. Every-step volume/clipping/radius checks; boundary checked at observations. Known coarse-radius limitation is reported separately from sensitivity.'}
    output.mkdir(parents=True)
    write_json(output/'protocol.json',protocol)
    write_json(output/'status.json',{'state':'prepared','protocol_sha256':sha(output/'protocol.json')})
    return protocol


def observation(sim):
    return {'metrics':sim.metrics(), 'ids':sim.ids.tolist(),
            **{key:getattr(sim,key).tolist() for key in ('activator','inhibitor','fate','polarity')},
            'labels':np.where(sim.fate>sim.config.fate_threshold,1,np.where(sim.fate < -sim.config.fate_threshold,-1,0)).tolist()}


def compare(output):
    output = Path(output)
    protocol = json.loads((output/'protocol.json').read_text())
    reports = {str(c):json.loads((output/f'cutoff-{c:g}'/'analysis.json').read_text()) for c in protocol['signal_cutoffs']}
    base = reports['0.02']['history']
    rows = []
    for cutoff in protocol['signal_cutoffs'][1:]:
        history = reports[str(cutoff)]['history']
        if len(history)!=len(base) or any(a['ids']!=b['ids'] or a['metrics']['time']!=b['metrics']['time'] for a,b in zip(history,base)):
            raise ValueError('observations must match identities and times')
        row = {'cutoff':cutoff}
        for species in ('activator','inhibitor'):
            row[species+'_max_rms_difference'] = max(float(np.sqrt(np.mean((np.array(a[species])-b[species])**2))) for a,b in zip(history,base))
            row[species+'_max_contrast_difference'] = max(abs(float(np.std(a[species])-np.std(b[species]))) for a,b in zip(history,base))
        row['max_axis_ratio_relative_difference'] = max(abs(a['metrics']['axis_ratio']/b['metrics']['axis_ratio']-1) for a,b in zip(history,base))
        row['max_continuous_fate_difference'] = max(float(np.max(abs(np.array(a['fate'])-b['fate']))) for a,b in zip(history,base))
        row['max_label_disagreements'] = max(int(np.count_nonzero(np.array(a['labels'])!=b['labels'])) for a,b in zip(history,base))
        rows.append(row)
    limits=protocol['criteria']
    checks = {'signals_within_0_01':all(max(r[k+'_max_rms_difference'] for k in ('activator','inhibitor'))<limits['signal_absolute_rms_max'] for r in rows),
              'contrast_within_0_005':all(max(r[k+'_max_contrast_difference'] for k in ('activator','inhibitor'))<limits['contrast_absolute_difference_max'] for r in rows),
              'axis_ratio_within_1_percent':all(r['max_axis_ratio_relative_difference']<limits['axis_ratio_relative_max'] for r in rows),
              'continuous_fate_within_0_05':all(r['max_continuous_fate_difference']<limits['continuous_fate_absolute_max'] for r in rows),
              'all_sampled_identity_labels_match':all(r['max_label_disagreements']==0 for r in rows)}
    quality = {'sampled_boundaries_clear':all(max(x['metrics']['boundary_occupancy'] for x in r['history'])<limits['boundary_occupancy_max'] for r in reports.values()),
               'every_step_volumes_within_5_percent':all(r['max_volume_error']<limits['cell_volume_error_max'] for r in reports.values()),
               'no_phase_field_clipping':all(r['max_clipped_fraction']==0 for r in reports.values()),
               'cells_resolved':all(r['min_radius_grid_cells']>=limits['minimum_radius_grid_cells'] for r in reports.values())}
    unchanged = sha(protocol['checkpoint'])==protocol['checkpoint_sha256']
    result = {'protocol':protocol,'pairs':rows,'sensitivity_checks':checks,'quality_checks':quality,
              'source_unchanged':unchanged,'sensitivity_passed':all(checks.values()),
              'all_checks_pass':all(checks.values()) and all(quality.values()) and unchanged}
    write_json(output/'comparison.json',result)
    lines=['# Dynamic signaling-cutoff sensitivity','',f'Sensitivity checks pass: **{result["sensitivity_passed"]}**',f'All sensitivity and quality checks pass: **{result["all_checks_pass"]}**','',protocol['scope'],'']
    lines += [f'- {key}: {"PASS" if value else "FAIL"}' for key,value in {**checks,**quality,'source_unchanged':unchanged}.items()]
    lines += ['', 'Detailed per-branch maxima are in comparison.json.']
    (output/'RESULTS.md').write_text('\n'.join(lines)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(10,7),constrained_layout=True)
    for cutoff,report in reports.items():
        times=[x['metrics']['time'] for x in report['history']]
        for ax,key in zip(axes.flat,('activator_std','inhibitor_std','axis_ratio','fate_a')):
            ax.plot(times,[x['metrics'][key] for x in report['history']],label=cutoff)
            ax.set(xlabel='Model time',ylabel=key.replace('_',' '))
    axes[0,0].legend(title='Signaling cutoff')
    fig.suptitle('Matched mature-state branches; polarity cutoff fixed at 0.02')
    fig.savefig(output/'comparison.png',dpi=140)
    plt.close(fig)
    return result


def run(output):
    output=Path(output)
    p=json.loads((output/'protocol.json').read_text())
    state=json.loads((output/'status.json').read_text())
    if state['state']!='prepared' or sha(output/'protocol.json')!=state['protocol_sha256']:
        raise ValueError('requires unchanged prepared protocol')
    if sha(p['checkpoint'])!=p['checkpoint_sha256'] or sha(Path(__file__).with_name('model.py'))!=p['model_sha256']:
        raise ValueError('source checkpoint or model changed since preparation')
    completed=[]
    try:
        for cutoff in p['signal_cutoffs']:
            path=output/f'cutoff-{cutoff:g}'
            path.mkdir()
            sim=Simulation.restore(p['checkpoint'])
            sim.config.graph_contact_cutoff=cutoff
            sim.config.polarity_contact_cutoff=p['polarity_cutoff']
            sim.config.steps=round(p['until']/sim.config.dt)
            sim.config.validate()
            history=[observation(sim)]
            frames=[{'metrics':history[0]['metrics'],'cells':sim.surfaces(max_points=200),'graph':sim.graph_snapshot()}]
            max_volume=history[0]['metrics']['max_cell_volume_error']
            min_radius=history[0]['metrics']['min_radius_grid_cells']
            clipped=sim.clipped_fraction
            start=time.monotonic()
            write_json(output/'status.json',{**state,'state':'running','cutoff':cutoff,'time':sim.time,'completed_cutoffs':completed})
            while sim.step_number<sim.config.steps:
                sim.step()
                volumes=sim.volumes()
                max_volume=max(max_volume,float(np.max(abs(volumes/sim.target-1))))
                min_radius=min(min_radius,float(np.min((3*volumes/(4*np.pi))**(1/3))/sim.dx))
                clipped=max(clipped,sim.clipped_fraction)
                if (sim.step_number-round(p['start']/sim.config.dt)) % p['sample_steps']==0:
                    history.append(observation(sim))
                    write_json(path/'history.json',history)
                    write_json(output/'status.json',{**state,'state':'running','cutoff':cutoff,'time':sim.time,'completed_cutoffs':completed})
                    if (len(history)-1)%10==0:
                        sim.checkpoint(path/'progress_state.npz')
                        print(f'cutoff={cutoff:g}, t={sim.time:g}, elapsed={time.monotonic()-start:.1f}s',flush=True)
            sim.checkpoint(path/'final_state.npz')
            frames.append({'metrics':history[-1]['metrics'],'cells':sim.surfaces(max_points=200),'graph':sim.graph_snapshot()})
            payload=json.dumps({'config':asdict(sim.config),'frames':frames},allow_nan=False)
            (path/'trajectory.json').write_text(payload)
            from .surface import viewer_template
            template=viewer_template().replace('A freely evolving 3D aggregate of deformable cells.','Matched cutoff intervention: two endpoint surfaces; intermediate metrics in history.json.')
            (path/'viewer.html').write_text(template.replace('__SIMULATION_DATA__',payload))
            write_json(path/'analysis.json',{'cutoff':cutoff,'config':asdict(sim.config),'history':history,
                       'elapsed_seconds':time.monotonic()-start,'max_volume_error':max_volume,
                       'min_radius_grid_cells':min_radius,'max_clipped_fraction':clipped})
            completed.append(cutoff)
        result=compare(output)
        write_json(output/'status.json',{**state,'state':'completed','completed_cutoffs':completed,
                   'sensitivity_passed':result['sensitivity_passed'],'all_checks_pass':result['all_checks_pass']})
        return result
    except Exception as error:
        write_json(output/'status.json',{**state,'state':'failed','completed_cutoffs':completed,'error':str(error)})
        raise


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['prepare','run'])
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--checkpoint',type=Path,default=Path('outputs/domain-conservative/grid-56/final_state.npz'))
    parser.add_argument('--until',type=float,default=90.)
    args=parser.parse_args()
    prepare(args.checkpoint,args.output,args.until) if args.action=='prepare' else run(args.output)


if __name__=='__main__':main()
