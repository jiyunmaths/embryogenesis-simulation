"""Short matched mechanical-response screen, not a long-time patterning sweep."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from .attribute_development import AttributeSimulation
from .feedback_spectrum import summarize


class ResponseSimulation(AttributeSimulation):
    law='linear'
    def material_coefficients(self):
        if self.law=='linear':
            tension,attraction=super().material_coefficients()
        else:
            c=self.config;r=np.tanh(self.activator-1)
            tension=c.surface_tension*np.exp(c.fate_tension*r)
            attraction=c.adhesion*np.exp(c.fate_adhesion*r[:,None]*r[None,:])
            np.fill_diagonal(attraction,0)
        mask=~np.eye(len(tension),dtype=bool)
        if np.any(tension<=0) or np.any(attraction[mask]<0):
            raise ValueError('Nonpositive tension or negative adhesion coefficient')
        return tension,attraction


def worker(args):
    root,arm=args;root=Path(root);p=json.loads((root/'protocol.json').read_text())
    out=root/arm['name'];out.mkdir()
    sim=ResponseSimulation.restore(p['checkpoint']);sim.law=arm['law']
    sim.attribute_mode='direct';sim.config.feedback=True
    sim.config.fate_tension=arm['gamma'];sim.config.fate_adhesion=1.4*arm['gamma']
    sim.config.polarity_tension=.35 if arm['polarity'] else 0.
    # Polarity dynamics retained in every arm; only its mechanical action differs.
    sim.config.dt=arm.get('dt',.0075)
    sim.step_number=int(round(sim.time/sim.config.dt))
    start=sim.time;steps=round(p['duration']/sim.config.dt)
    records=[];max_volume=0.;max_clip=0.;min_radius=float('inf');started=time.monotonic()
    for step in range(steps+1):
        v=sim.volumes();max_volume=max(max_volume,float(np.max(abs(v/sim.target-1))))
        min_radius=min(min_radius,float(np.min((3*v/(4*np.pi))**(1/3))/sim.dx));max_clip=max(max_clip,sim.clipped_fraction)
        if step%round(.15/sim.config.dt)==0:
            row,_=summarize(sim);row['elapsed']=sim.time-start
            row['log_activator_sd']=float(np.std(np.log(sim.activator)))
            row['minimum_tension']=float(sim.material_coefficients()[0].min())
            records.append(row)
            (out/'history.json').write_text(json.dumps(records,indent=2)+'\n')
            (out/'status.json').write_text(json.dumps(dict(step=step,steps=steps))+'\n')
        if step<steps:sim.step()
    metrics=sim.metrics()
    quality=dict(max_volume_error=max_volume,max_clipping=max_clip,min_radius=min_radius,
                 boundary=float(metrics['boundary_occupancy']))
    quality['pass']=bool(max_volume<.05 and max_clip==0 and min_radius>=4 and quality['boundary']<.01)
    final=records[-1];result=dict(arm=arm,quality=quality,initial=records[0],final=final,
        growth_change=final['maximum_spatial_growth']-records[0]['maximum_spatial_growth'],
        total_conductance_change=final['total_conductance']-records[0]['total_conductance'],
        elapsed_seconds=time.monotonic()-started)
    # Explicit law metadata; experimental state must not be restored as an ordinary linear AttributeSimulation.
    np.savez_compressed(out/'response_state.npz',phi=sim.phi,activator=sim.activator,inhibitor=sim.inhibitor,
                        polarity=sim.polarity,volumes=sim.volumes(),ids=sim.ids,law=np.array(sim.law))
    (out/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def run(checkpoint,output):
    arms=[dict(name='baseline',law='linear',gamma=0.,polarity=False),
          dict(name='polarity_only',law='linear',gamma=0.,polarity=True),
          dict(name='historical_full',law='linear',gamma=.25,polarity=True)]
    arms += [dict(name=f'linear_{g:g}',law='linear',gamma=g,polarity=False) for g in [.1,.25,.5,1.]]
    arms += [dict(name=f'exponential_{g:g}',law='exponential',gamma=g,polarity=False) for g in [.1,.25,.5,1.,1.5,2.]]
    p=dict(checkpoint=str(checkpoint.resolve()),duration=.6,arms=arms,
           scope='Common no-feedback t=90 state; freely evolving chemistry/mechanics/polarity; only polarity mechanical action is ablated. Short-response screen cannot establish long-time resonance or developmental pattern selection. Exponential extension changes the constitutive law.',
           checkpoint_sha256=hashlib.sha256(checkpoint.read_bytes()).hexdigest(),
           source_sha256={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in [Path(__file__),Path(__file__).with_name('attribute_development.py'),Path(__file__).with_name('model.py')]})
    output.mkdir(parents=True,exist_ok=False);(output/'protocol.json').write_text(json.dumps(p,indent=2)+'\n')
    with ProcessPoolExecutor(max_workers=2) as pool:
        rows=list(pool.map(worker,[(str(output),arm) for arm in arms]))
    (output/'results.json').write_text(json.dumps(rows,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--checkpoint',type=Path,default=Path('outputs/attribute-development/no_feedback/state-90.npz'))
    p.add_argument('--output',type=Path,default=Path('outputs/feedback-response'))
    args=p.parse_args();run(args.checkpoint,args.output)
