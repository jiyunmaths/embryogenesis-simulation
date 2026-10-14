"""Literal prepared-array/source reconstruction; no project helper imports."""
import hashlib
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[3]
STUDY=ROOT/'output/studies/phase-carry-recipient-delayed-response'


def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


p=read(STUDY/'protocol.json');old=read(Path(p['parent'])/'protocol.json')
watched={**p['source_sha256'],**p['input_sha256'],str(STUDY/'protocol.json'):sha(STUDY/'protocol.json')}
for name,expected in watched.items():assert sha(name)==expected,name
assert p['start']==570. and p['wait']==60. and p['original_matching_time']==510.
assert (p['beta'],p['da'],p['db'])==(2.,.02,.55)
with np.load(Path(p['parent'])/'production/fine/paths.npz') as z:
    times=z['times'].copy();parent={str(k):v.copy() for k,v in zip(z['trial_keys'],z['trajectory'].transpose(1,0,2,3))}
    assert np.array_equal(times,np.arange(401)*.15)
other={}
for label in ('coarse','diagnostic'):
    with np.load(Path(p['parent'])/'production'/label/'paths.npz') as z:
        assert np.array_equal(z['times'],times)
        other[label]={str(k):v[-1].copy() for k,v in zip(z['trial_keys'],z['trajectory'].transpose(1,0,2,3))}
contexts={};sources=[];settled=0
for c in p['contexts']:
    original=next(x for x in old['contexts'] if x['key']==c['parent_context'])
    control=next(t for t in old['trials'] if t['context']==original['key'] and t['target'] is None)
    with np.load(original['source']) as z:previous={k:z[k].copy() for k in ('delta','volumes','ids')}
    with np.load(c['source']) as z:data={k:z[k].copy() for k in ('delta','volumes','ids','initial')}
    for key in ('delta','volumes','ids'):assert np.array_equal(data[key],previous[key])
    u=data['initial'];delta=data['delta'];v=data['volumes'];path=parent[control['key']]
    assert np.array_equal(u,path[-1]) and u.dtype==np.float64 and np.all(u>0)
    assert len(set(data['ids']))==len(v) and u.shape==(2,len(v))
    assert np.max(abs(delta.sum(axis=1)))<1e-10
    assert np.max(abs((v[:,None]*delta)-(v[:,None]*delta).T))<1e-10
    assert np.max(abs(v@delta))<1e-10
    assert np.array_equal(u@v,np.asarray(c['challenge_global_species_amounts']))
    assert np.array_equal(u*v,np.asarray(c['challenge_per_cell_species_amounts']))
    mask=times>=36.
    drift=float(np.max(abs(np.log(path[mask]/u[None]))))
    a,b=u;rates=np.stack((a*a/b-a+.02*(delta@a),2.*(a*a-b)+.55*(delta@b)))/u
    rate=float(np.max(abs(rates)))
    record=c['settling']
    assert abs(drift-record['preceding_max_raw_log_drift'])<=1e-14
    assert abs(rate-record['challenge_max_absolute_log_rate'])<=1e-14
    flag=bool(drift<=.001 and rate<=.0001)
    assert flag==record['descriptive_settling_screen'];settled+=flag
    ends={label:other[label][control['key']] for label in other}
    with np.load(Path(p['parent'])/'references'/control['key']/'paths.npz') as z:
        assert np.array_equal(z['times'],times)
        ends.update({k:z[k][-1].copy() for k in ('dop853','radau')})
    errors={k:float(np.max(abs(np.log(u/x)))) for k,x in ends.items()}
    assert all(x<=1e-8 for x in errors.values())
    assert errors==record['endpoint_log_discrepancies']
    contexts[c['key']]=data
    sources.append(dict(context=c['key'],parent_control=control['key'],endpoint_errors=errors,
        preceding_max_log_drift=drift,challenge_max_log_rate=rate,descriptive_settling_screen=flag))

pairs={};physical=[]
for f in p['factorials']:
    for P in ('S','R'):
        names=[f['contexts'][C+'_'+P] for C in ('S','R')]
        a,b=[contexts[k] for k in names]
        for key in ('delta','volumes','ids'):assert np.array_equal(a[key],b[key])
        i=a['ids'].tolist().index(f['cell'])
        values=[float(a['initial'][0,i]),float(b['initial'][0,i])]
        for key in names:pairs[key,f['cell']]=values
        physical.append(dict(seed=f['seed'],cell=f['cell'],matching_reference=f['matching_reference'],physics=P,
            target_species_ratio=(b['initial'][:,i]/a['initial'][:,i]).tolist(),
            surrounding_global_species_amount_ratio=((b['initial']@b['volumes'])/(a['initial']@a['volumes'])).tolist()))

cases=[];unique=set();maximum_amount_error=0.
for t in p['trials']:
    c=contexts[t['context']];u=c['initial'];expected=u.copy()
    with np.load(t['source']) as z:actual=z['initial'].copy()
    if t['target'] is None:
        assert t['factor']==1. and t['initial_log_amplitude']==0.
    else:
        i=c['ids'].tolist().index(t['target']);v=c['volumes'][i]
        minimum=min(pairs[t['context'],t['target']])
        factor=.9 if t['dose']=='fractional' else float(1.-.1*minimum/float(u[0,i]))
        expected[0,i]*=factor
        requested=.1*float(u[0,i]*v) if t['dose']=='fractional' else float(.1*minimum*v)
        measured=float((u[0,i]-expected[0,i])*v)
        error=abs(measured/requested-1);maximum_amount_error=max(maximum_amount_error,error)
        assert factor==t['factor'] and float(abs(np.log(factor)))==t['initial_log_amplitude']
        assert requested==t['requested_removed_amount'] and measured==t['removed_activator_amount']
        assert t['initial_target_volume']==float(v) and t['initial_target_activator_amount']==float(u[0,i]*v)
        assert error==t['pulse_amount_relative_error'] and error<=2e-13
    assert np.array_equal(actual,expected)
    signature=(t['context'],*((x.dtype.str,x.shape,x.tobytes()) for x in (c['delta'],c['volumes'],c['ids'],actual)))
    unique.add(signature);cases.append(dict(trial=t['key'],source_sha256=sha(t['source']),exact_preparation=True))
assert (len(contexts),len(cases),sum(t['target'] is None for t in p['trials']))==(45,141,45)
assert len(unique)==p['unique_initial_value_problems']
for name,expected in watched.items():assert sha(name)==expected,name
result=dict(passed=True,protocol_sha256=sha(STUDY/'protocol.json'),contexts=45,logical_cases=141,
    controls=45,pulses=96,unique_initial_value_problems=len(unique),source_checks=sources,
    source_gate_count=45,descriptively_settled_backgrounds=settled,paired_context_checks=physical,cases=cases,
    maximum_pulse_relative_amount_error=maximum_amount_error,
    inputs_unchanged=True,inputs_sha256=watched,script_sha256=sha(__file__),
    scope='Actual unpulsed t570 states preserved exactly, no rematching; both delayed doses checked; whole-network amounts and target states unmatched.')
(Path(__file__).with_name('input-audit.json')).write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ('inputs_sha256','source_checks','paired_context_checks','cases')},indent=2))
