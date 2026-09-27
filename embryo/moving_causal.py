"""Matched moving-geometry feedback controls using verified joint fate stages."""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from .model import Simulation
from .joint_fate import advance
from .causal_signaling import initial_signals, preflight
from .shape import observe
from .resolution import write_json, _steps


ARMS=('full','no_feedback','no_self_activation','no_signal_to_fate')


class JointMovingSimulation(Simulation):
    def step(self, *, prescribed_signals=None):
        if prescribed_signals is not None:raise ValueError('joint experiment does not accept signal clamps')
        if self.divisions or len(self.phi)!=self.config.max_cells:raise ValueError('joint moving pilot requires no further division')
        graph=self.signaling_graph()
        kinetics='full' if self.arm=='no_feedback' else self.arm
        self.joint_state=advance(self.joint_state,graph.delta,self.config,kinetics,self.config.dt)
        self.activator=1+self.joint_state[0];self.inhibitor=1+self.joint_state[1];self.fate=self.joint_state[2].copy()
        # This pathway reuses mechanics and polarity ordering, while skipping
        # the base reaction/fate update and applying dilution exactly once below.
        super().step(prescribed_signals=(self.activator,self.inhibitor))

    def update_fate(self,contacts,exposure):
        pass # Fate was advanced at the same RK stages as both regulators.

    def mechanical_step(self):
        before=self.volumes()
        super().mechanical_step()
        after=self.volumes();ratio=before/after
        for k in (0,1):
            self.joint_state[k]=self.joint_state[k]*ratio+(before-after)/after
        self.activator=1+self.joint_state[0];self.inhibitor=1+self.joint_state[1]
        if np.any(self.activator<=0) or np.any(self.inhibitor<=0):raise FloatingPointError('nonpositive diluted regulator')

    def graph_snapshot(self,contacts=None):
        result=super().graph_snapshot(contacts)
        graph=self.signaling_graph(contacts)
        arm='full' if self.arm=='no_feedback' else self.arm
        p=preflight(graph.eigenvalues,self.config.signal_beta,self.config.signal_da,self.config.signal_dh,arm)
        growth=np.array(p['growth_rates']);nonzero=graph.eigenvalues>1e-10
        result.update(intervention=self.arm,local_growth_rate=p['local_max_growth'],local_stable=p['local_stable'],
                      growth_rates=p['growth_rates'],unstable_modes=p['unstable_spatial_modes'],
                      diffusion_driven_instability=bool(p['local_stable'] and p['unstable_spatial_modes']),
                      maximum_spatial_growth=float(growth[nonzero].max()) if nonzero.any() else None)
        if arm=='no_self_activation':result['continuous_lambda_band']=None
        return result

    def metrics(self):
        result=super().metrics()
        result['unstable_graph_modes']=len(self.graph_snapshot()['unstable_modes'])
        result['causal_arm']=self.arm
        return result

    def checkpoint(self,path):
        super().checkpoint(path)
        with np.load(path,allow_pickle=False) as saved:payload={key:saved[key].copy() for key in saved.files}
        payload.update(joint_state=self.joint_state,causal_arm=np.array(self.arm))
        np.savez_compressed(path,**payload)

    @classmethod
    def restore(cls,path):
        sim=super().restore(path)
        with np.load(path,allow_pickle=False) as saved:
            sim.arm=str(saved['causal_arm']) if 'causal_arm' in saved else 'full'
            sim.joint_state=saved['joint_state'].copy() if 'joint_state' in saved else np.stack([sim.activator-1,sim.inhibitor-1,sim.fate])
        return sim


def initialize(checkpoint,arm,seed):
    if arm not in ARMS:raise ValueError('unknown moving intervention')
    sim=JointMovingSimulation.restore(checkpoint);sim.arm=arm
    sim.config.feedback=arm!='no_feedback'
    a,h=initial_signals(sim.volumes(),[seed],.001)
    sim.activator=a[0];sim.inhibitor=h[0];sim.fate=np.zeros(len(sim.phi));sim.polarity[:]=0
    sim.joint_state=np.stack([sim.activator-1,sim.inhibitor-1,sim.fate])
    return sim


def prepare(checkpoint,validation,output,duration=60.,interval=.6):
    checkpoint=Path(checkpoint).resolve();validation=Path(validation).resolve();output=Path(output)
    if output.exists():raise FileExistsError('choose a fresh output directory')
    gate=json.loads((validation/'comparison.json').read_text())
    if not gate['passed']:raise ValueError('joint fate reference validation must pass first')
    if hashlib.sha256(checkpoint.read_bytes()).hexdigest()!=gate['protocol']['checkpoint_sha256']:raise ValueError('moving source must match validated geometry')
    sim=Simulation.restore(checkpoint)
    if sim.divisions or len(sim.phi)!=sim.config.max_cells or sim.config.fate_noise or sim.config.neighbor_inhibition or sim.config.exposure_bias:
        raise ValueError('requires mature source without additional fate forcing')
    if sim.config.signal_transport!='conservative' or sim.metrics()['min_radius_grid_cells']<4:raise ValueError('requires resolved conservative source')
    _steps(duration,sim.config.dt);_steps(interval,sim.config.dt)
    if round(duration/sim.config.dt)%round(interval/sim.config.dt):raise ValueError('observation times must align')
    p={'checkpoint':str(checkpoint),'checkpoint_sha256':hashlib.sha256(checkpoint.read_bytes()).hexdigest(),
       'validation':str(validation),'validation_sha256':hashlib.sha256((validation/'comparison.json').read_bytes()).hexdigest(),
       'arms':list(ARMS),'seed':7,'start':sim.time,'duration':duration,'interval':interval,'late_duration':15.,
       'source_config':asdict(sim.config),'criteria':{'shape_excess_min':.05,'persistent_contrast_min':.1,'volume_max':.05,'boundary_max':.01,'radius_min':4.},
       'source_sha256':{str(Path(__file__).with_name(n).resolve()):hashlib.sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in ('model.py','joint_fate.py','moving_causal.py')},
       'scope':'Single matched mature-geometry perturbation pilot. Reset signals to paired unbiased deviations, fate to zero, polarity to zero; retain geometry/lineage. Joint signal/fate RK2, standard polarity/mechanics, exact amount dilution. Geometry already anisotropic: test feedback-specific excess, not spontaneous embryogenesis. No cleavage. Core live solver unchanged.'}
    output.mkdir(parents=True);write_json(output/'protocol.json',p);write_json(output/'status.json',{'state':'prepared','protocol_sha256':hashlib.sha256((output/'protocol.json').read_bytes()).hexdigest()})
    return p


def compare(output,p):
    reports={arm:json.loads((output/arm/'analysis.json').read_text()) for arm in p['arms']}
    late={arm:[r for r in report['history'] if r['time']>=p['start']+p['duration']-p['late_duration']-1e-10] for arm,report in reports.items()}
    ratios={arm:float(np.mean([r['axis_ratio'] for r in rows])) for arm,rows in late.items()}
    excess={arm:ratios['full']-ratios[arm] for arm in p['arms'] if arm!='full'}
    effects={'feedback_specific_shape_excess':excess['no_feedback']>=p['criteria']['shape_excess_min'],
             'self_activation_specific_shape_excess':excess['no_self_activation']>=p['criteria']['shape_excess_min'],
             'full_persistent_signal_contrast':min(r['activator_std'] for r in late['full'])>=p['criteria']['persistent_contrast_min']}
    quality={arm:{'volume':r['max_volume_error']<p['criteria']['volume_max'],
                  'radius':r['min_radius']>=p['criteria']['radius_min'],'no_clipping':r['max_clipping']==0,
                  'boundary':max(h['boundary_occupancy'] for h in r['history'])<p['criteria']['boundary_max']} for arm,r in reports.items()}
    result={'protocol':p,'late_axis_ratios':ratios,'shape_excess_over_controls':excess,'causal_effect_checks':effects,'quality_checks':quality,
            'all_quality_pass':all(all(q.values()) for q in quality.values()),'interpretation':'Mechanism checks may fail scientifically even when numerics pass. One geometry/seed, no developmental or ensemble causality claim.'}
    write_json(output/'comparison.json',result)
    lines=['# Moving-geometry feedback controls','',result['interpretation'],'','| Arm | Late mean axis ratio | Full minus control |','|---|---:|---:|']
    lines += [f'| {arm} | {ratios[arm]:.6f} | {excess.get(arm,0):.6f} |' for arm in p['arms']]
    lines += ['',*[f'- {k}: {v}' for k,v in effects.items()],f'- All quality checks pass: {result["all_quality_pass"]}']
    (output/'RESULTS.md').write_text('\n'.join(lines)+'\n')
    return result


def run(output):
    output=Path(output);p=json.loads((output/'protocol.json').read_text());state=json.loads((output/'status.json').read_text())
    if state['state']!='prepared' or hashlib.sha256((output/'protocol.json').read_bytes()).hexdigest()!=state['protocol_sha256']:raise ValueError('requires unchanged prepared protocol')
    for path,digest in {**p['source_sha256'],p['checkpoint']:p['checkpoint_sha256'],str(Path(p['validation'])/'comparison.json'):p['validation_sha256']}.items():
        if hashlib.sha256(Path(path).read_bytes()).hexdigest()!=digest:raise ValueError('prepared input changed')
    completed=[]
    try:
        for arm in p['arms']:
            path=output/arm;path.mkdir();sim=initialize(p['checkpoint'],arm,p['seed'])
            steps=_steps(p['duration'],sim.config.dt);every=_steps(p['interval'],sim.config.dt)
            sim.config.steps=sim.step_number+steps;sim.config.save_every=every
            history=[];frames=[];max_volume=0.;min_radius=float('inf');max_clipping=0.;started=time.monotonic()
            for i in range(steps+1):
                volume=sim.volumes();max_volume=max(max_volume,float(np.max(abs(volume/sim.target-1))))
                min_radius=min(min_radius,float(np.min((3*volume/(4*np.pi))**(1/3))/sim.dx));max_clipping=max(max_clipping,sim.clipped_fraction)
                if i%every==0:
                    row=observe(sim);row.update(activator=sim.activator.tolist(),inhibitor=sim.inhibitor.tolist(),fate=sim.fate.tolist())
                    history.append(row);write_json(path/'history.json',history)
                    write_json(output/'status.json',{**state,'state':'running','arm':arm,'time':sim.time,'completed_arms':completed})
                    if i%(10*every)==0 or i==steps:
                        sim.checkpoint(path/f'state-{sim.time:g}.npz');frames.append({'metrics':row,'cells':sim.surfaces(max_points=200),'graph':sim.graph_snapshot()})
                        print(f'{arm}: t={sim.time:g}, elapsed={time.monotonic()-started:.1f}s',flush=True)
                if i<steps:sim.step()
            sim.checkpoint(path/'final_state.npz')
            write_json(path/'analysis.json',{'history':history,'max_volume_error':max_volume,'min_radius':min_radius,'max_clipping':max_clipping,'elapsed_seconds':time.monotonic()-started})
            from .surface import viewer_template
            payload=json.dumps({'config':asdict(sim.config),'frames':frames},allow_nan=False)
            (path/'viewer.html').write_text(viewer_template().replace('__SIMULATION_DATA__',payload));completed.append(arm)
        result=compare(output,p);write_json(output/'status.json',{**state,'state':'completed','completed_arms':completed,'all_quality_pass':result['all_quality_pass']})
        return result
    except Exception as error:
        write_json(output/'status.json',{**state,'state':'failed','completed_arms':completed,'error':str(error)});raise


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=['prepare','run']);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--checkpoint',type=Path,default=Path('outputs/development-refinement/space-72/state-18.npz'));parser.add_argument('--validation',type=Path,default=Path('outputs/joint-fate-validation'))
    args=parser.parse_args();prepare(args.checkpoint,args.validation,args.output) if args.action=='prepare' else run(args.output)
