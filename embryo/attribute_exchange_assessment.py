"""Post-run stability and ordering diagnostics for the chemical exchange assay."""
import argparse
import json
from pathlib import Path
import numpy as np


def assess(output):
    result=json.loads((output/'results.json').read_text())
    config=json.loads((output/'protocol.json').read_text())['config']
    data=np.load(output/'trajectories.npz')
    initial=data['initial'];final=data['exchanges'][:,-1];delta=data['delta'];n=len(delta)
    diagnostic={}
    for outcome in ['destination','reorganized']:
        selected=[k for k,row in enumerate(result['pairs']) if row['outcome']==outcome]
        diagnostic[outcome]=dict(count=len(selected),reversed_activator_ordering=sum(
            bool((initial[0,i]-initial[0,j])*(final[k,0,i]-final[k,0,j])<0)
            for k in selected for i,j in [result['pairs'][k]['indices']]),
            closer_to_transferred_than_destination=sum(
                result['pairs'][k]['transferred_ratio']<result['pairs'][k]['destination_ratio']
                for k in selected))
    growth=[]
    for a,b in np.concatenate([final,data['seeded_resets'][:,-1]]):
        jac=np.block([[np.diag(2*a/b-1)+config['signal_da']*delta,np.diag(-a*a/b**2)],
            [np.diag(2*config['signal_beta']*a),-config['signal_beta']*np.eye(n)+config['signal_dh']*delta]])
        growth.append(float(np.linalg.eigvals(jac).real.max()))
    diagnostic['largest_real_jacobian_eigenvalues']=growth
    diagnostic['all_tested_final_states_locally_stable']=bool(max(growth)<0)
    diagnostic['seeded_resets_within_1e_5_of_original']=sum(
        row['destination_distance']<1e-5 for row in result['seeded_uniform_resets'])
    diagnostic['scope']='Post-run descriptive diagnostics; not predeclared outcome criteria. Final-state Jacobians ordered as 120 exchanges then 20 seeded resets.'
    (output/'assessment.json').write_text(json.dumps(diagnostic,indent=2)+'\n')
    return diagnostic


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path('outputs/attribute-exchange'))
    assess(parser.parse_args().output)
