"""Literal analytic-source and control audit, without candidate solver imports."""
from pathlib import Path
import hashlib,json
import numpy as np
from scipy.special import expit
from scipy.integrate import quad

ROOT=Path('output/studies/development-gpu-cleavage-qualification-v2')
HERE=Path(__file__).parent
def sha(f):return hashlib.sha256(Path(f).read_bytes()).hexdigest()


def main():
    p=json.loads((ROOT/'protocol.json').read_text());watched={};rows=[];targets=[];random_states=[]
    for f,h in {**p['source_sha256'],**p['input_sha256']}.items():assert sha(f)==h,f;watched[f]=h
    for j in p['jobs']:
        c=j['config'];assert c['grid']==72 and c['max_cells']==2 and c['signal_transport']=='conservative' and not c['differentiation']
        assert c['feedback']==(j['branch']!='off')
        assert c['fate_tension']==(.25 if j['branch'] in ('tension_adhesion','full') else 0.)
        assert c['fate_adhesion']==(.35 if j['branch'] in ('tension_adhesion','full') else 0.)
        assert c['polarity_tension']==(.35 if j['branch'] in ('polarity','full') else 0.)
        with np.load(j['source'],allow_pickle=False) as z:
            meta=json.loads(str(z['metadata']));assert meta['time']==0 and meta['step_number']==0 and meta['divisions']=={} and meta['config']==c
            assert z['phi'].dtype==np.float32 and z['phi'].shape==(1,72,72,72)
            assert np.all(z['phase_carry']==0) and np.all(z['chemical_correction']==0)
            assert z['phase_carry'].dtype==z['chemical_correction'].dtype==np.float64
            for key,value in [('activator',[1.03]),('inhibitor',[.98]),('ids',[0]),('parents',[-1]),('fate',[0.])]:np.testing.assert_array_equal(z[key],value)
            np.testing.assert_array_equal(z['polarity'],[[0.,0.,0.]] if j['axis']=='axial' else [[.08,.03,-.02]])
            dx=2*c['extent']/c['grid'];axis=(np.arange(c['grid'])+.5)*dx-c['extent'];xyz=np.stack(np.meshgrid(axis,axis,axis,indexing='ij'))
            radius=np.sqrt(np.sum(xyz**2,axis=0));np.testing.assert_array_equal(z['phi'][0],expit(np.sqrt(2)*(.4-radius)/c['interface_width']).astype(np.float32))
            w=c['interface_width']
            def h(q):return q*q*(3-2*q)
            target=quad(lambda r:4*np.pi*r*r*h(expit(np.sqrt(2)*(.4-r)/w)),0,.4+30*w,epsabs=1e-12,epsrel=1e-12)[0]
            np.testing.assert_array_equal(z['target'],[target]);targets.append(target)
            assert str(z['chemical_method'])==p['method'] and int(z['development_schema'])==1
            random_states.append([meta[k] for k in ('rng','fate_rng','signal_rng')]);rows.append(dict(case=j['key'],initial_exact=True,target=target))
    assert len(set(targets))==1 and all(s==random_states[0] for s in random_states)
    for f,h in watched.items():assert sha(f)==h,f
    result=dict(passed=True,cases=16,new_scientific_histories=0,source_pins=len(p['source_sha256']),input_pins=len(p['input_sha256']),
        protocol_sha256=sha(ROOT/'protocol.json'),records=rows,matched_random_states=True,checked_sha256=watched,script_sha256=sha(__file__))
    (HERE/'input-audit.json').write_text(json.dumps(result,indent=2)+'\n');print('All initial arrays, controls, RNG states and pins pass.')


if __name__=='__main__':main()
