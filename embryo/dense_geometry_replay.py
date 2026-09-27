"""Capture every-step geometry from exact continuation, then assess chemical replay."""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from .geometry_replay import ARMS, trajectory
from .joint_fate import discrepancies
from .moving_causal import JointMovingSimulation
from .resolution import _steps, write_json


def run(source, output):
    source=Path(source).resolve();output=Path(output)
    if output.exists():raise FileExistsError('choose a fresh output directory')
    source_protocol=json.loads((source/'protocol.json').read_text())
    for path,digest in source_protocol['source_sha256'].items():
        if hashlib.sha256(Path(path).read_bytes()).hexdigest()!=digest:
            raise ValueError('moving source code changed; cannot claim exact continuation')
    if json.loads((source/'status.json').read_text())['state']!='completed':raise ValueError('source incomplete')
    sim=JointMovingSimulation.restore(source/'full/state-18.npz')
    p={'source':str(source),'dt':sim.config.dt,'start':sim.time,'duration':60.,'arms':list(ARMS),
       'capture_interval':sim.config.dt,'snapshot_strides':[1,2], 'replay_dts':[sim.config.dt,sim.config.dt/2],
       'criteria':{'signal_rms_max':.01,'fate_difference_max':.05},
       'input_sha256':{str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted((source/'full').glob('state-*.npz'))+[source/'full/history.json']},
       'source_sha256':{**source_protocol['source_sha256'],**{str(Path(__file__).with_name(n).resolve()):hashlib.sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in ['dense_geometry_replay.py','geometry_replay.py']}},
       'scope':'Exact continuation from full t18 checkpoint, every-step conductance/volume capture, bitwise source checkpoint comparisons, followed by four chemistry-only replays. No altered equations or retuning.'}
    output.mkdir(parents=True);write_json(output/'protocol.json',p);write_json(output/'status.json',{'state':'capturing','time':sim.time})
    try:
        steps=_steps(p['duration'],sim.config.dt);every=_steps(6.,sim.config.dt)
        times=[];conductance=[];volumes=[];audits=[];initial=sim.joint_state.copy();started=time.monotonic()
        for i in range(steps+1):
            graph=sim.signaling_graph();times.append(sim.time);conductance.append(graph.weights.copy());volumes.append(sim.volumes())
            if i%every==0:
                reference=JointMovingSimulation.restore(source/'full'/f'state-{sim.time:g}.npz')
                audit={'time':sim.time,'phase_bitwise_equal':bool(np.array_equal(sim.phi,reference.phi)),
                       'joint_bitwise_equal':bool(np.array_equal(sim.joint_state,reference.joint_state)),
                       'polarity_bitwise_equal':bool(np.array_equal(sim.polarity,reference.polarity))}
                audits.append(audit)
                if not all(v for k,v in audit.items() if k!='time'):raise ValueError(f'source checkpoint mismatch: {audit}')
                np.savez_compressed(output/'geometry.npz',times=times,conductance=conductance,volumes=volumes,initial=initial,ids=sim.ids)
                write_json(output/'capture-audit.json',audits)
                write_json(output/'status.json',{'state':'capturing','time':sim.time,'elapsed_seconds':time.monotonic()-started})
                print(f'Captured t={sim.time:g}; checkpoint identical; elapsed={time.monotonic()-started:.1f}s',flush=True)
            if i<steps:sim.step()
        times=np.array(times);conductance=np.array(conductance);volumes=np.array(volumes)
        history=json.loads((source/'full/history.json').read_text());observations=np.array([r['time'] for r in history])
        recorded=np.array([[np.array(r['activator'])-1,np.array(r['inhibitor'])-1,r['fate']] for r in history])
        write_json(output/'status.json',{'state':'replaying'})
        rows=[];runs={};summary={}
        for arm in ARMS:
            for stride in p['snapshot_strides']:
                for dt in p['replay_dts']:
                    value=trajectory(times[::stride],conductance[::stride],volumes[::stride],initial,sim.config,arm,dt,observations)
                    runs[arm,stride,dt]=value
                    np.savez_compressed(output/f'{arm}-stride-{stride}-dt-{dt}.npz',times=observations,state=value)
            coarse=runs[arm,1,p['replay_dts'][0]];fine=runs[arm,1,p['replay_dts'][1]]
            rows.append({'arm':arm,'step_halving':discrepancies(coarse,fine),
                         'snapshot_decimation':discrepancies(fine,runs[arm,2,p['replay_dts'][1]])})
            late=observations>=63-1e-10;contrast=fine[:,0].std(axis=1)
            labels=np.where(fine[-1,2]>.55,1,np.where(fine[-1,2]<-.55,-1,0))
            summary[arm]={'late_contrast_min':float(contrast[late].min()),'late_contrast_max':float(contrast[late].max()),
                          'final_a':int(sum(labels==1)),'final_b':int(sum(labels==-1)),'final_uncommitted':int(sum(labels==0))}
            write_json(output/'progress.json',{'rows':rows,'summary':summary});print(arm,rows[-1],flush=True)
        fidelity={str(dt):discrepancies(runs['replay',1,dt],recorded) for dt in p['replay_dts']}
        c=p['criteria'];okay=lambda d:d['max_signal_rms']<c['signal_rms_max'] and d['max_fate_difference']<c['fate_difference_max'] and d['final_labels_match']
        checks={'capture_checkpoint_identity':all(all(v for k,v in a.items() if k!='time') for a in audits),
                'all_step_halving_checks':all(okay(r['step_halving']) for r in rows),
                'all_snapshot_decimation_checks':all(okay(r['snapshot_decimation']) for r in rows),
                'both_replays_match_recorded_full':all(okay(d) for d in fidelity.values())}
        result={'protocol':p,'checks':checks,'attribution_ready':all(checks.values()),'rows':rows,'fidelity':fidelity,'summary':summary,
                'scope':'One recorded geometry/seed; diagnostic interventions may violate moving-volume conservation. Halving chemistry dt does not refine the mechanical source trajectory.'}
        write_json(output/'comparison.json',result)
        lines=['# Dense geometry replay','',f'Attribution prerequisites pass: **{result["attribution_ready"]}**','',
               '| Arm | Late contrast min–max | Final A / B / uncommitted |','|---|---:|---:|']
        for arm,r in summary.items():lines.append(f'| {arm} | {r["late_contrast_min"]:.6g}–{r["late_contrast_max"]:.6g} | {r["final_a"]} / {r["final_b"]} / {r["final_uncommitted"]} |')
        lines+=['',*[f'- {k}: {v}' for k,v in checks.items()], '',result['scope']]
        (output/'RESULTS.md').write_text('\n'.join(lines)+'\n');write_json(output/'status.json',{'state':'completed','attribution_ready':result['attribution_ready']})
        return result
    except Exception as error:
        write_json(output/'status.json',{'state':'failed','error':str(error)});raise


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--source',type=Path,default=Path('outputs/moving-causal'));parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();run(args.source,args.output)
