"""Compare selected matched feedback response at two mechanical timesteps."""
import json
from pathlib import Path
import numpy as np


def compare(coarse=Path('outputs/feedback-response/historical_full'),fine=Path('outputs/feedback-response-refinement/historical_full_half_dt')):
    a=json.loads((coarse/'history.json').read_text());b=json.loads((fine/'history.json').read_text())
    if len(a)!=len(b) or not np.allclose([x['time'] for x in a],[x['time'] for x in b]):
        raise ValueError('Unmatched sample times')
    growth=max(abs(x['maximum_spatial_growth']-y['maximum_spatial_growth']) for x,y in zip(a,b))
    c=np.load(coarse/'response_state.npz');f=np.load(fine/'response_state.npz')
    chemical=max(float(np.abs(np.log(c[key]/f[key])).max()) for key in ['activator','inhibitor'])
    result=dict(max_sampled_growth_rate_difference=growth,max_endpoint_chemical_log_difference=chemical,
                growth_limit=1e-4,chemical_limit=1e-4,passed=bool(growth<1e-4 and chemical<1e-4),
                scope='Selected historical full-coupling response only; does not certify every coupling or spatial resolution.')
    (fine.parent/'comparison.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__=='__main__':
    print(json.dumps(compare(),indent=2))
