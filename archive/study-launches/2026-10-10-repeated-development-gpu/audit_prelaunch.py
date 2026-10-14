"""Independent initial-field/control/pin audit before trajectories are launched."""
from pathlib import Path
import hashlib,json
import numpy as np
from scipy.integrate import quad
from scipy.special import expit

ROOT=Path('output/studies/repeated-development-gpu-qualification-v2')
HERE=Path(__file__).parent


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()


if __name__=='__main__':
    p=json.loads((ROOT/'protocol.json').read_text());pins={**p['source_sha256'],**p['input_sha256']}
    preservation=json.loads(Path(p['preserved_manifest']).read_text())
    pins.update(preservation['preserved_sha256'])
    for f,h in pins.items():assert sha(f)==h,f
    states=[];rows=[];targets=[]
    target=quad(lambda r:4*np.pi*r*r*(lambda x:x*x*(3-2*x))(expit(np.sqrt(2)*(.8-r)/.085)),
        0,.8+30*.085,epsabs=1e-12,epsrel=1e-12)[0]
    for j in p['jobs']:
        c=j['config'];assert c['signal_transport']=='conservative' and not c['differentiation'] and c['division_orientation']=='shape'
        assert c['max_cells']==16 and c['division_interval']==2. and c['signal_da']==.02 and c['signal_dh']==.55
        assert c['feedback']==(j['branch']!='off')
        assert c['fate_tension']==(.25 if j['branch'] in ('tension_adhesion','full') else 0.)
        assert c['fate_adhesion']==(.35 if j['branch'] in ('tension_adhesion','full') else 0.)
        assert c['polarity_tension']==(.35 if j['branch'] in ('polarity','full') else 0.)
        with np.load(j['source'],allow_pickle=False) as z:
            m=json.loads(str(z['metadata']));assert m['time']==0 and m['step_number']==0 and m['divisions']=={} and m['config']==c
            n=c['grid'];dx=2*c['extent']/n;axis=(np.arange(n)+.5)*dx-c['extent']
            radius=np.sqrt(sum(x*x for x in np.meshgrid(axis,axis,axis,indexing='ij')))
            np.testing.assert_array_equal(z['phi'][0],expit(np.sqrt(2)*(.8-radius)/.085).astype(np.float32))
            np.testing.assert_array_equal(z['target'],[target]);targets.append(float(z['target'][0]))
            for key,value in [('activator',[1.03]),('inhibitor',[.98]),('ids',[0]),('parents',[-1]),('fate',[0.])]:np.testing.assert_array_equal(z[key],value)
            assert np.all(z['polarity']==0) and np.all(z['phase_carry']==0) and np.all(z['chemical_correction']==0)
            assert z['phi'].dtype==np.float32 and z['phase_carry'].dtype==z['chemical_correction'].dtype==np.float64
            states.append([m[k] for k in ('rng','fate_rng','signal_rng')]+[z['due'].tolist()])
            rows.append(dict(case=j['key'],analytic_exact=True,initial_attributes_exact=True))
    assert len(set(targets))==1 and all(s==states[0] for s in states)
    calibration=json.loads((ROOT/'calibration.json').read_text());cm=json.loads((ROOT/'calibration-manifest.json').read_text())
    assert cm['calibration_sha256']==sha(ROOT/'calibration.json')
    assert calibration['protocol_sha256']==cm['protocol_sha256']==sha(ROOT/'protocol.json')
    assert len(calibration['rows'])==36 and calibration['passed'] and max(r['relative_l2'] for r in calibration['rows'])<=.0025
    for f,h in pins.items():assert sha(f)==h,f
    result=dict(passed=True,protocol_sha256=sha(ROOT/'protocol.json'),cases=16,new_scientific_histories=0,
        source_pins=len(p['source_sha256']),input_pins=len(p['input_sha256']),preservation_pins=len(preservation['preserved_sha256']),
        calibration_max=max(r['relative_l2'] for r in calibration['rows']),matched_initial_random_streams=True,rows=rows,
        checked_sha256=pins,calibration_sha256=sha(ROOT/'calibration.json'),script_sha256=sha(__file__))
    (HERE/'prelaunch-audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('rows','checked_sha256')},indent=2))
