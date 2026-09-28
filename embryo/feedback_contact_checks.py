"""Endpoint contact-closure diagnostics for the feedback spectrum analysis."""
import argparse
import json
from pathlib import Path
import numpy as np
from .attribute_development import AttributeSimulation
from .signaling import stability
from .transport import contact_transport,transport_graph


def run(source,output):
    rows=[]
    for mode in ['direct','no_feedback']:
        sim=AttributeSimulation.restore(source/mode/'state-90.npz')
        contacts,_=sim.contacts();v=sim.volumes();centers=sim.centers();c=sim.config
        for cutoff in [0.,.005,.01,.02,.04,.08]:
            graph=transport_graph(contact_transport(contacts,v,centers,c.interface_width,cutoff))
            report=stability(graph,c.signal_beta,c.signal_da,c.signal_dh)
            rows.append(dict(mode=mode,cutoff=cutoff,edges=report['edges'],max_growth=report['maximum_spatial_growth'],
                             max_eigenvalue=float(graph.eigenvalues[-1]),unstable_modes=report['unstable_modes']))
    (output/'cutoff_checks.json').write_text(json.dumps(rows,indent=2)+'\n')
    history=json.loads((source/'direct/history.json').read_text())
    amplitudes=[]
    for time in [6,12,18,36,60,90]:
        frame=min(history,key=lambda x:abs(x['metrics']['time']-time))
        a=np.exp([x['attributes']['log_activator'] for x in frame['cells']])
        polarity=np.array([x['attributes']['polarity_magnitude'] for x in frame['cells']])
        amplitudes.append(dict(time=frame['metrics']['time'],
            max_activity_tension_fraction=float(np.max(.25*abs(np.tanh(a-1)))),
            max_polarity_directional_modulation_bound=float(np.max(.35*polarity))))
    (output/'modulation_amplitudes.json').write_text(json.dumps(amplitudes,indent=2)+'\n')
    return rows


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,default=Path('outputs/attribute-development'))
    p.add_argument('--output',type=Path,default=Path('outputs/feedback-spectrum'))
    args=p.parse_args();run(args.source,args.output)
