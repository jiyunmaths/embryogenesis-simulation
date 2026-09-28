"""Instantaneous mechanical force interventions on one common chemical state."""
import argparse
import json
from pathlib import Path
import numpy as np
from .attribute_development import AttributeSimulation
from .feedback_spectrum import summarize


def run(checkpoint,output):
    rows=[]
    for dt in [.0075,.00375]:
        for name,cg,ca,cp in [('baseline',0,0,0),('tension',.25,0,0),('adhesion',0,.35,0),('activity',.25,.35,0),('polarity',0,0,.35),('full',.25,.35,.35)]:
            sim=AttributeSimulation.restore(checkpoint);sim.config.feedback=True
            sim.config.fate_tension=cg;sim.config.fate_adhesion=ca;sim.config.polarity_tension=cp;sim.config.dt=dt
            before,graph=summarize(sim);v0=sim.volumes();a=sim.activator.copy();c0=sim.centers()
            sim.mechanical_step() # No chemistry, polarity update, or dilution in this tangent probe.
            after,_=summarize(sim)
            rows.append(dict(name=name,dt=dt,growth_derivative=(after['maximum_spatial_growth']-before['maximum_spatial_growth'])/dt,
                conductance_derivative=(after['total_conductance']-before['total_conductance'])/dt,
                activity_volume_rate_correlation=float(np.corrcoef(a,(sim.volumes()/v0-1)/dt)[0,1]),
                max_center_speed=float(np.linalg.norm(sim.centers()-c0,axis=1).max()/dt)))
    for row in rows:
        base=next(x for x in rows if x['dt']==row['dt'] and x['name']=='baseline')
        row['excess_growth_derivative']=row['growth_derivative']-base['growth_derivative']
    (output/'force_split.json').write_text(json.dumps(rows,indent=2)+'\n')
    return rows


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--checkpoint',type=Path,default=Path('outputs/attribute-development/no_feedback/state-90.npz'))
    p.add_argument('--output',type=Path,default=Path('outputs/feedback-spectrum'))
    args=p.parse_args();run(args.checkpoint,args.output)
