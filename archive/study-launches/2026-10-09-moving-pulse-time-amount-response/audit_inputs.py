"""Independent literal input reconstruction; no embryo imports or stepping."""
import hashlib,json
from pathlib import Path
import numpy as np
ROOT=Path('output/studies/moving-pulse-time-amount-response')
WATCHED={}
def sha(f):return hashlib.sha256(Path(f).read_bytes()).hexdigest()
def pin(f,h=None):
    name=str(Path(f).resolve());now=sha(name);assert h is None or now==h,name
    assert name not in WATCHED or WATCHED[name]==now,name;WATCHED[name]=now;return name
def read(f):return json.loads(Path(pin(f)).read_text())
def arrays(f,keys):
    with np.load(pin(f),allow_pickle=False) as z:return {k:z[k].copy() for k in keys}
def main():
    p=read(ROOT/'protocol.json');parent=read(Path(p['parent'])/'protocol.json');local=read(Path(p['local_parent'])/'protocol.json');mp=read(Path(p['physical_parent'])/'protocol.json')
    reviews={}
    for kind,file in [('challenge','docs/pulse_time_amount_response_assessment.json'),('moving','docs/moving_recipient_amount_response_assessment.json')]:
        r=read(file);a=read(pin(r['independent_review'],r['independent_review_sha256']));assert r['full_study_accepted'] and a['passed'] and a['accepted_at_tested_scope'];reviews[kind]=a
    for key in ('source_sha256','input_sha256'):
        for f,h in p[key].items():pin(f,h)
    assert (len(p['contexts']),len(p['pairs']),len(p['jobs']))==(12,6,48)
    assert p['histories']==[7,8,9] and p['independent_histories']==3 and p['new_histories']==0
    assert p['primary_reference']==p['budget_reference']==p['initial_physics']=='S'
    assert p['start']==570. and p['wait']==0. and p['response_duration']==60. and p['pilot_duration']==6.
    assert p['method']=='compensated SSP-RK2 and error-free dilution, phase-carry CUDA v1'
    assert p['levels']==[dict(level='coarse',dt=.00375),dict(level='fine',dt=.001875)]
    states={};starts={};closures=[];signatures=set();physical_hashes={}
    for c in p['contexts']:
        old=next(x for x in parent['contexts'] if x['key']==c['key']);assert all(c[k]==old[k] for k in old if k!='source') and c['accepted_context_source']==old['source']
        assert (c['matching_reference'],c['budget_reference'],c['physical_origin'])==('S','S','S') and c['recipient']==c['targets'][0]
        z=arrays(c['source'],('initial','volumes','delta','ids'));before=arrays(c['accepted_context_source'],tuple(z));states[c['key']]=z
        for k in z:np.testing.assert_array_equal(z[k],before[k])
        pin(c['accepted_context_source'],reviews['challenge']['evidence_sha256'][str(Path(c['accepted_context_source']).resolve())])
        d,v,u,ids=[z[k] for k in ('delta','volumes','initial','ids')];assert u.dtype==v.dtype==d.dtype==np.float64 and u.shape==(2,16) and np.all(u>0) and np.isfinite(u).all()
        assert np.max(abs(d.sum(1)))<1e-11 and np.max(abs(v@d))<1e-12 and np.max(abs(v[:,None]*d-(v[:,None]*d).T))<1e-12
        pc=next(x for x in local['contexts'] if x['key']==c['parent_context']);assert c['physical_parent_context']==pc['physical_parent_context']
        job=next(j for j in mp['jobs'] if (j['context'],j['stage'],j['level'])==(pc['physical_parent_context'],'wait','fine'))
        source=Path(p['physical_parent'])/job['key']/'latest_state.npz';assert c['physical_source']==str(source)
        for name in ('history.json','result.json','latest_state.npz'):pin(source.with_name(name),reviews['moving']['evidence_sha256'][str(source.with_name(name).resolve())])
        h=read(source.with_name('history.json'));r=read(source.with_name('result.json'));assert r['quality_pass'] and r['job']==job and h[-1]['time']==570. and h[-1]['elapsed']==60.
        for k in ('delta','volumes','ids'):np.testing.assert_array_equal(z[k],h[-1][k])
        names=('phi','phase_carry','rounding','polarity','target','ids','parents','due','fate','metadata','attribute_mode','precision_arm','chemical_correction','chemical_method')
        x=arrays(source,names);meta=json.loads(str(x['metadata']));assert meta['time']==570. and not meta['divisions']
        assert x['phi'].dtype==np.float32 and x['phi'].shape==(16,72,72,72) and x['phase_carry'].dtype==np.float64 and x['phase_carry'].shape==x['phi'].shape and np.isfinite(x['phase_carry']).all()
        assert x['rounding'].dtype==np.int32 and x['rounding'].shape==(2,) and np.all(x['rounding']>=0) and str(x['precision_arm'])=='phase_carry' and str(x['attribute_mode'])=='direct' and np.all(x['fate']==0)
        assert str(x['chemical_method'])==p['method'];np.testing.assert_array_equal(ids,x['ids'])
        cfg=meta['config'];assert cfg['signal_transport']=='conservative' and cfg['feedback'] and not cfg['differentiation']
        assert tuple(cfg[k] for k in ('signal_beta','signal_da','signal_dh','fate_tension','fate_adhesion','polarity_tension'))==(2.,.02,.55,.25,.35,.35)
        for level,dt in [('coarse',.00375),('fine',.001875)]:
            retimed=dict(cfg,dt=dt,steps=round(cfg['steps']*cfg['dt']/dt),save_every=round(cfg['save_every']*cfg['dt']/dt));assert retimed==p['accepted_configs'][c['key']+'_'+level]
        ph=hashlib.sha256()
        for name in names[:-2]:
            a=x[name];ph.update(name.encode());ph.update(a.dtype.str.encode());ph.update(str(a.shape).encode());ph.update(a.tobytes())
        physical_hashes[c['key']]=ph.hexdigest()
    for j in p['jobs']:
        c=next(c for c in p['contexts'] if c['key']==j['context']);z=states[c['key']];i=z['ids'].tolist().index(c['recipient'])
        assert j['stage']=='response' and j['start']==570. and j['duration']==60. and j['dt']==dict((x['level'],x['dt']) for x in p['levels'])[j['level']]
        assert (j['target'] is None)==(j['dose'] is None) and j['dose'] in (None,'fractional')
        s=arrays(j['source'],('initial','correction','delta','volumes','ids'));record=read(Path(j['source']).with_suffix('.json'));assert record['physical_source']==c['physical_source'] and record['paths_sha256']==sha(j['source']) and record['physical_source_sha256']==sha(c['physical_source'])
        expected=z['initial'].copy()
        if j['dose']:expected[0,i]*=.9
        np.testing.assert_array_equal(s['initial'],expected);assert s['correction'].dtype==np.float64 and s['correction'].shape==expected.shape and np.all(s['correction']==0)
        for k in ('delta','volumes','ids'):np.testing.assert_array_equal(s[k],z[k])
        trial=next(t for t in parent['trials'] if t['context']==c['key'] and t['target']==j['target']);assert record['parent_trial']==trial['key'] and record['record']['factor']==trial['factor']
        pin(trial['source'],reviews['challenge']['evidence_sha256'][str(Path(trial['source']).resolve())]);original=arrays(trial['source'],('initial',));np.testing.assert_array_equal(s['initial'],original['initial'])
        signatures.add((physical_hashes[c['key']],s['initial'].tobytes()));starts[c['key'],j['dose']]=s
        for duration in (.15,.6,3.,6.,60.,570.):assert abs(duration/j['dt']-round(duration/j['dt']))<1e-8
    for pair in p['pairs']:
        c=[next(c for c in p['contexts'] if c['key']==pair['contexts'][C]) for C in ('S','R')];assert c[0]['physical_source']==c[1]['physical_source']
        for dose in (None,'fractional'):
            a,b=[starts[pair['contexts'][C],dose] for C in ('S','R')];v=a['volumes'];i=a['ids'].tolist().index(pair['cell']);others=np.arange(16)!=i
            np.testing.assert_array_equal(a['initial'][:,i],b['initial'][:,i])
            errors=[float(np.max(abs(x/y-1))) for x,y in [(a['initial']@v,b['initial']@v),(a['initial'][:,others]@v[others],b['initial'][:,others]@v[others])]];assert max(errors)<=2e-13
            closures.append(dict(history=pair['seed'],cell=pair['cell'],dose=dose,global_relative_error=errors[0],surrounding_relative_error=errors[1]))
    assert len(signatures)==p['distinct_initial_value_problems']
    for f,h in WATCHED.items():assert sha(f)==h,f
    result=dict(passed=True,protocol_sha256=sha(ROOT/'protocol.json'),contexts=12,nested_recipients=6,total_paths=48,logical_starts=24,distinct_initial_value_problems=len(signatures),independent_histories=3,new_histories=0,matched_budget_checks=closures,maximum_initial_amount_relative_error=max(max(x['global_relative_error'],x['surrounding_relative_error']) for x in closures),physical_sources_referenced=len(set(c['physical_source'] for c in p['contexts'])),inputs_unchanged=True,inputs_sha256=WATCHED,script_sha256=sha(__file__),historical_noninput_export_exception=p['parent_historical_export_note'])
    Path(__file__).with_name('preparation-audit.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps({k:result[k] for k in ('passed','total_paths','distinct_initial_value_problems','maximum_initial_amount_relative_error')},indent=2))
if __name__=='__main__':main()
