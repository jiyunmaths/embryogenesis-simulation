"""Matched tension-only and adhesion-only short continuations."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from .feedback_response import ResponseSimulation
from .feedback_spectrum import summarize


def run(checkpoint,output):
    output.mkdir(parents=True,exist_ok=False)
    (output/'protocol.json').write_text(json.dumps(dict(checkpoint=str(checkpoint),duration=.6,dt=.0075,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        checkpoint_sha256=hashlib.sha256(checkpoint.read_bytes()).hexdigest(),
        scope='Same initial state as feedback-response; tension-only versus adhesion-only; full signal/dilution evolution, polarity dynamics retained but no polarity tension.'),indent=2)+'\n')
    results=[]
    for name,cg,ca in [('tension_only',.25,0.),('adhesion_only',0.,.35)]:
        sim=ResponseSimulation.restore(checkpoint);sim.config.feedback=True;sim.attribute_mode='direct'
        sim.config.fate_tension=cg;sim.config.fate_adhesion=ca;sim.config.polarity_tension=0.
        h=[];maxerr=0.;clip=0.;minradius=float('inf')
        for step in range(81):
            v=sim.volumes();maxerr=max(maxerr,float(np.max(abs(v/sim.target-1))))
            minradius=min(minradius,float(np.min((3*v/(4*np.pi))**(1/3))/sim.dx));clip=max(clip,sim.clipped_fraction)
            if step%20==0:
                row,_=summarize(sim);h.append(row)
                (output/(name+'_history.json')).write_text(json.dumps(h,indent=2)+'\n')
            if step<80:sim.step()
        boundary=sim.metrics()['boundary_occupancy']
        results.append(dict(name=name,history=h,growth_change=h[-1]['maximum_spatial_growth']-h[0]['maximum_spatial_growth'],
            max_volume_error=maxerr,max_clipping=clip,min_radius=minradius,boundary=boundary,
            quality_pass=bool(maxerr<.05 and clip==0 and minradius>=4 and boundary<.01)))
        (output/'results.json').write_text(json.dumps(results,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--checkpoint',type=Path,default=Path('outputs/attribute-development/no_feedback/state-90.npz'))
    p.add_argument('--output',type=Path,default=Path('outputs/feedback-components'))
    args=p.parse_args();run(args.checkpoint,args.output)
