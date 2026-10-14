"""Reconstruct all actual-start gates before the moving science starts."""
import json,hashlib
from pathlib import Path
import numpy as np

ROOT=Path('output/studies/moving-recipient-amount-response')
WATCHED={}
def sha(f):return hashlib.sha256(Path(f).read_bytes()).hexdigest()
def pin(f,h=None):
    name=str(Path(f).resolve());now=sha(name);assert h is None or now==h,name
    assert name not in WATCHED or WATCHED[name]==now,name;WATCHED[name]=now;return name
def read(f):return json.loads(Path(pin(f)).read_text())
def arrays(f,keys):
    with np.load(pin(f),allow_pickle=False) as z:return {k:z[k].copy() for k in keys}
def errors(a,b):
    d,e=np.array(a['delta']),np.array(b['delta'])
    return dict(chemical_log_max=float(abs(np.log(np.array(a['chemistry'])/b['chemistry'])).max()),polarity_abs_max=float(abs(np.array(a['polarity'])-b['polarity']).max()),relative_axis_max=abs(a['axis_ratio']/b['axis_ratio']-1),relative_volume_max=float(abs(np.array(a['volumes'])/b['volumes']-1).max()),relative_transport_max=float(np.linalg.norm(d-e)/max(np.linalg.norm(e),1e-300)))

def main():
    p=read(ROOT/'protocol.json');ph=sha(ROOT/'protocol.json');q=read(ROOT/'qualification.json');assert q['passed'] and q['jobs']==24 and q['method']==p['method'] and q['protocol_sha256']==ph
    for key in ('source_sha256','input_sha256'):
        for f,h in p[key].items():pin(f,h)
    times=np.arange(5)*.15;rows=[];actual_starts={}
    for job in [j for j in p['jobs'] if j['stage']=='wait']:
        folder=ROOT/job['key']/'prefix';r=read(folder/'result.json');assert r['passed'] and r['job']==job and r['protocol_sha256']==ph and r['source_sha256']==sha(job['source']) and r['method']==p['method']
        for f,h in r['evidence_sha256'].items():pin(f,h)
        gpu=read(folder/'gpu-history.json');cpu=read(folder/'cpu-history.json');start=arrays(job['source'],('initial','volumes','ids','delta'))
        assert len(gpu)==len(cpu)==5
        for h in (gpu,cpu):
            assert np.allclose([x['elapsed'] for x in h],times,rtol=0,atol=1e-10)
            assert np.array_equal(h[0]['chemistry'],start['initial'])
            for x,t in zip(h,times):
                assert np.array_equal(x['ids'],start['ids']) and abs(x['time']-job['start']-t)<1e-9
                u=np.array(x['chemistry']);v=np.array(x['volumes']);d=np.array(x['delta'])
                assert np.isfinite(u).all() and np.all(u>0) and np.all(v>0) and x['boundary_occupancy']<.01
                assert np.max(abs(d.sum(1)))<1e-11 and np.max(abs(v@d))<1e-12 and np.max(abs(v[:,None]*d-(v[:,None]*d).T))<1e-12
        measured={key:max(errors(a,b)[key] for a,b in zip(gpu,cpu)) for key in r['errors']}
        for key,value in measured.items():assert value==r['errors'][key] and value<=p['prefix_limits'][key]
        a=arrays(folder/'gpu-phi.npz',('phi',))['phi'];b=arrays(folder/'cpu-phi.npz',('phi',))['phi'];phi=float(abs(a-b).max())
        assert phi==r['final_phi_abs_max'] and phi<=p['prefix_limits']['final_phi_abs_max']
        assert r['independent_chemical_steps']==8 and r['cpu_chemical_log_max']<=1e-11 and r['amount_accounting_max']<=2e-14 and r['restart_exact_steps']==4
        cp=arrays(folder/'restart.npz',('metadata','activator','inhibitor','chemical_correction','chemical_method','phi','phase_carry','rounding','precision_arm'))
        u=np.array([cp['activator'],cp['inhibitor']]);e=cp['chemical_correction'];assert e.dtype==np.float64 and e.shape==u.shape and np.isfinite(e).all() and np.all(abs(e)<=2*np.spacing(u)) and str(cp['chemical_method'])==p['method']
        assert str(cp['precision_arm'])=='phase_carry' and cp['phase_carry'].dtype==np.float64 and cp['phase_carry'].shape==cp['phi'].shape and np.isfinite(cp['phase_carry']).all()
        meta=json.loads(str(cp['metadata']));assert abs(meta['time']-job['start']-8*job['dt'])<1e-9
        # Confirm that actual GPU geometry reproduces the capacity/operator used
        # to construct the accepted initial chemical intervention.
        volume_error=float(abs(np.array(gpu[0]['volumes'])/start['volumes']-1).max());operator_error=float(np.linalg.norm(np.array(gpu[0]['delta'])-start['delta'])/np.linalg.norm(start['delta']))
        assert volume_error<=2e-13 and operator_error<=1e-12
        actual_starts[job['context'],job['level']]=gpu[0]
        rows.append(dict(job=job['key'],errors=measured,final_phi_abs_max=phi,chemical_step_error=r['cpu_chemical_log_max'],dilution_error=r['amount_accounting_max'],initial_capacity_error=volume_error,initial_operator_error=operator_error,checkpoint_method_and_corrections_verified=True,runtime_exact_restart_steps=r['restart_exact_steps'],passed=True))
    for pair in p['pairs']:
        for level in ('coarse','fine'):
            a,b=[actual_starts[pair['contexts'][C],level] for C in ('S','R')];v=np.array(a['volumes']);i=a['ids'].index(pair['cell']);others=np.arange(16)!=i;x,y=np.array(a['chemistry']),np.array(b['chemistry'])
            assert np.array_equal(a['volumes'],b['volumes']) and np.array_equal(a['delta'],b['delta']) and np.array_equal(x[:,i],y[:,i])
            assert np.max(abs((x@v)/(y@v)-1))<=2e-13 and np.max(abs((x[:,others]@v[others])/(y[:,others]@v[others])-1))<=2e-13
    guards=read('archive/study-document-snapshots/2026-10-08-before-moving-recipient-amount-response/manifest.json')['unchanged_guards']
    for f,h in guards.items():pin(f,h)
    for f,h in WATCHED.items():assert sha(f)==h
    result=dict(passed=True,method=p['method'],protocol_sha256=ph,actual_start_gates=24,independent_histories=3,new_histories=0,rows=rows,inputs_unchanged=True,inputs_sha256=WATCHED,script_sha256=sha(__file__),scope='Separate actual-start saved-array gate reconstruction. Exact restart and eight-step chemistry remain source/runtime evidence; no full-horizon scientific result or native equivalence.')
    Path(__file__).with_name('qualification-audit.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(dict(passed=True,actual_starts=24,max_chemical_error=max(x['chemical_step_error'] for x in rows),max_dilution_error=max(x['dilution_error'] for x in rows)),indent=2))

if __name__=='__main__':main()
