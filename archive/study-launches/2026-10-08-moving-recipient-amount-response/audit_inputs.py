"""Literal preparation audit; no project helper imports or stepping."""
import json,hashlib
from pathlib import Path
import numpy as np

ROOT=Path('output/studies/moving-recipient-amount-response')
WATCHED={}
def sha(file):return hashlib.sha256(Path(file).read_bytes()).hexdigest()
def pin(file,expected=None):
    name=str(Path(file).resolve());h=sha(name)
    assert expected is None or h==expected,name
    assert name not in WATCHED or WATCHED[name]==h,name
    WATCHED[name]=h;return name
def read(file):return json.loads(Path(pin(file)).read_text())
def data(file,keys):
    with np.load(pin(file),allow_pickle=False) as z:return {k:z[k].copy() for k in keys}

def main():
    p=read(ROOT/'protocol.json');parent=read(Path(p['parent'])/'protocol.json')
    for collection in ('source_sha256','input_sha256'):
        for f,h in p[collection].items():pin(f,h)
    old={c['key']:c for c in parent['contexts']};assert (len(p['contexts']),len(p['pairs']),len(p['jobs']))==(12,6,96)
    assert p['primary_reference']==p['initial_physics']=='S' and p['histories']==[7,8,9] and p['new_histories']==0
    assert p['method']=='compensated SSP-RK2 and error-free dilution, phase-carry CUDA v1'
    assert p['levels']==[dict(level='coarse',dt=.00375),dict(level='fine',dt=.001875)]
    assert sum(j['stage']=='wait' for j in p['jobs'])==24 and sum(j['stage']=='response' for j in p['jobs'])==72
    states={};physical={};closures=[]
    for c in p['contexts']:
        assert c['matching_reference']=='S' and c['physics_label']=='S' and c['chemical_origin'] in ('S','R')
        assert all(c[k]==old[c['key']][k] for k in old[c['key']])
        z=data(c['source'],('delta','volumes','initial','ids'));states[c['key']]=z
        assert z['initial'].dtype==z['volumes'].dtype==z['delta'].dtype==np.float64 and z['initial'].shape==(2,16) and np.all(z['initial']>0) and np.isfinite(z['initial']).all()
        d=z['delta'];v=z['volumes'];assert np.max(abs(d.sum(1)))<1e-11 and np.max(abs(v@d))<1e-12 and np.max(abs(v[:,None]*d-(v[:,None]*d).T))<1e-12
        x=data(c['physical_source'],('phi','phase_carry','rounding','metadata','ids','fate','attribute_mode','precision_arm'));meta=json.loads(str(x['metadata']))
        assert meta['time']==510. and x['phi'].dtype==np.float32 and x['phi'].shape==(16,72,72,72) and x['phase_carry'].dtype==np.float64 and x['phase_carry'].shape==x['phi'].shape
        assert np.array_equal(x['ids'],z['ids']) and np.isfinite(x['phase_carry']).all() and np.all(x['fate']==0) and str(x['attribute_mode'])=='direct' and str(x['precision_arm'])=='phase_carry'
        assert not meta['divisions'] and meta['config']['signal_transport']=='conservative'
        assert tuple(meta['config'][k] for k in ('signal_beta','signal_da','signal_dh','fate_tension','fate_adhesion','polarity_tension'))==(2.,.02,.55,.25,.35,.35)
        physical.setdefault(c['seed'],c['physical_source']);assert physical[c['seed']]==c['physical_source']
    for pair in p['pairs']:
        a,b=[states[pair['contexts'][C]] for C in ('S','R')];i=a['ids'].tolist().index(pair['cell']);others=np.arange(16)!=i
        for k in ('delta','volumes','ids'):np.testing.assert_array_equal(a[k],b[k])
        np.testing.assert_array_equal(a['initial'][:,i],b['initial'][:,i]);v=a['volumes']
        surrounding=float(np.max(abs((b['initial'][:,others]@v[others])/(a['initial'][:,others]@v[others])-1)))
        whole=float(np.max(abs((b['initial']@v)/(a['initial']@v)-1)));assert max(surrounding,whole)<=2e-13
        closures.append(dict(history=pair['seed'],cell=pair['cell'],surrounding_relative_error=surrounding,whole_relative_error=whole))
    for j in p['jobs']:
        assert j['duration']==60. and j['dt']==dict((x['level'],x['dt']) for x in p['levels'])[j['level']]
        assert j['start']==(510. if j['stage']=='wait' else 570.)
        assert (j['target'] is None)==(j['dose'] is None)
        for duration in (.15,.6,3.,6.,60.,j['start']):assert abs(duration/j['dt']-round(duration/j['dt']))<1e-8
        if j['stage']=='wait':assert j['source']==next(c['source'] for c in p['contexts'] if c['key']==j['context'])
        else:assert not Path(j['source']).exists() # future actual fine endpoint, never a copied frozen result
    guards=read('archive/study-document-snapshots/2026-10-08-before-moving-recipient-amount-response/manifest.json')['unchanged_guards']
    for f,h in guards.items():pin(f,h)
    for f,h in WATCHED.items():assert sha(f)==h
    result=dict(passed=True,protocol_sha256=sha(ROOT/'protocol.json'),histories=3,new_histories=0,contexts=12,nested_recipients=6,wait_paths=24,response_paths=72,physical_sources_referenced=3,primary_reference='S',initial_physics='S',matched_budget_checks=closures,maximum_initial_amount_relative_error=max(max(x['surrounding_relative_error'],x['whole_relative_error']) for x in closures),historical_guards_unchanged=len(guards),inputs_unchanged=True,inputs_sha256=WATCHED,script_sha256=sha(__file__))
    Path(__file__).with_name('preparation-audit.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ('passed','contexts','wait_paths','response_paths','maximum_initial_amount_relative_error','historical_guards_unchanged')},indent=2))

if __name__=='__main__':main()
