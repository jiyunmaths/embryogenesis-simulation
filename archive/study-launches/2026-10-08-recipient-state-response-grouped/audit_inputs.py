"""Literal recipient-matching audit; no project imports or integration."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

WATCHED={}
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def pin(p,h=None):
    p=str(Path(p).resolve());actual=digest(p)
    assert h is None or actual==h,p
    assert p not in WATCHED or WATCHED[p]==actual,p
    WATCHED[p]=actual;return p
def read(p):return json.loads(Path(pin(p)).read_text())
def arrays(p):
    with np.load(pin(p),allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}
def signature(z,u=None):
    parts=[z[k] for k in ('delta','volumes','ids')]+[z['initial'] if u is None else u]
    return tuple((a.dtype.str,a.shape,a.tobytes()) for a in parts)

def audit(root):
    root=Path(root).resolve();p=read(root/'protocol.json');old=read(Path(p['parent'])/'protocol.json')
    compact=read('docs/phase_carry_shared_context_response_assessment.json')
    assert compact['full_study_accepted'] and compact['protocol_sha256']==digest(Path(p['parent'])/'protocol.json')
    proof=read(compact['independent_review']);assert proof['passed'] and proof['inputs_unchanged']
    assert digest(compact['independent_review'])==compact['independent_review_sha256']
    for field in ('source_sha256','input_sha256'):
        for f,h in p[field].items():pin(f,h)
    assert (p['beta'],p['da'],p['db'],p['start'],p['duration'],p['interval'])==(2.,.02,.55,510.,60.,.15)
    assert [(l['level'],l['dt']) for l in p['levels']]==[('coarse',.00375),('fine',.001875),('diagnostic',.0009375)]
    assert p['reference_settings']==dict(rtol=1e-13,atol=1e-15,max_step=.15)
    assert (len(p['contexts']),len(p['trials']),len(p['factorials']))==(45,93,12)
    donors={c['key']:arrays(c['source']) for c in old['contexts']}
    for c in old['contexts']:pin(c['source'],proof['inputs_sha256'][c['source']])
    prepared={c['key']:arrays(c['source']) for c in p['contexts']}
    actual={signature(z):key for key,z in prepared.items()};assert len(actual)==45
    cases={};intended=set();reset_details=[];pulse_errors=[];physical_checks=[]
    for f in old['factorials']:
        seed,cell=f['seed'],f['cell']
        for L in ('S','R'):
            reference=donors[f['contexts'][L+'_'+L]];i=reference['ids'].tolist().index(cell)
            item=next(x for x in p['factorials'] if (x['seed'],x['cell'],x['matching_reference'])==(seed,cell,L))
            for C in ('S','R'):
                for P in ('S','R'):
                    donor=donors[f['contexts'][C+'_'+P]];assert np.array_equal(donor['ids'],reference['ids'])
                    expected={k:v.copy() for k,v in donor.items()};expected['initial'][:,i]=reference['initial'][:,i]
                    sig=signature(expected);intended.add(sig);key=actual[sig]
                    assert item['contexts'][C+'_'+P]==key
                    cases[f'{seed}_{cell}_{L}_{C}_{P}']=key
    assert intended==set(actual)
    context_meta={c['key']:c for c in p['contexts']}
    for key,z in prepared.items():
        c=context_meta[key];donor=donors[c['parent_context']];u=z['initial'];v=z['volumes'];d=z['delta'];ids=z['ids']
        assert u.dtype==v.dtype==d.dtype==np.float64 and u.shape==(2,16) and d.shape==(16,16) and len(set(ids.tolist()))==16
        assert np.isfinite(u).all() and np.all(u>0) and np.isfinite(v).all() and np.all(v>0)
        g=v[:,None]*d;off=d.copy();np.fill_diagonal(off,0.)
        assert np.all(off>=0) and abs(d.sum(axis=1)).max()<1e-11
        assert abs(v@d).max()<1e-12 and abs(g-g.T).max()<1e-12
        counts=[max(1,int(np.ceil(l['dt']*max(1+.02*float(-np.diag(d).min()),2+.55*float(-np.diag(d).min()))/.2))) for l in p['levels']]
        assert counts==[1,1,1];physical_checks.append(dict(context=key,substeps=counts))
        before=donor['initial']*v;after=u*v;record=c['matching_accounting']
        for name,value in [('pre_matching_per_cell_species_amounts',before),('post_matching_per_cell_species_amounts',after),('matching_per_cell_species_amount_change',after-before),('pre_matching_global_species_amounts',before.sum(axis=1)),('post_matching_global_species_amounts',after.sum(axis=1)),('matching_global_species_amount_change',(after-before).sum(axis=1))]:
            assert np.array_equal(record[name],value),(key,name)
        assert c['donor_amount_accounting']==old['assignments'][c['parent_context']]
        for match in c['matching_targets']:
            i=ids.tolist().index(match['cell'])
            assert np.array_equal(match['original_pair'],donor['initial'][:,i])
            assert np.array_equal(match['matched_pair'],u[:,i])
        reset_details.append(dict(context=key,matching_reference=c['matching_reference'],global_species_amount_change=(after-before).sum(axis=1).tolist(),native=bool(np.array_equal(u,donor['initial']))))
    starts={};controls={};native_cases=[];unique=set()
    for t in p['trials']:
        z=prepared[t['context']];u=arrays(t['source'])['initial'];expected=z['initial'].copy()
        if t['target'] is None:controls[t['context']]=t['key'];assert t['factor']==1.
        else:
            i=z['ids'].tolist().index(t['target']);expected[0,i]*=.9
            assert t['dose']=='fractional' and t['factor']==.9 and t['initial_log_amplitude']==abs(np.log(.9))
            requested=.1*float(z['initial'][0,i]*z['volumes'][i]);removed=float((z['initial'][0,i]-u[0,i])*z['volumes'][i])
            assert requested==t['requested_removed_amount'] and removed==t['removed_activator_amount']
            error=abs(removed/requested-1.);assert error<=2e-13;pulse_errors.append(error)
        assert np.array_equal(u,expected);starts[t['key']]=u;unique.add(signature(z,u))
        if 'native_parent_trial' in t:
            prior=next(x for x in old['trials'] if x['key']==t['native_parent_trial'])
            assert prior['context']==context_meta[t['context']]['parent_context'] and prior['target']==t['target'] and prior['dose']==t['dose']
            np.testing.assert_array_equal(arrays(prior['source'])['initial'],u);native_cases.append(t['key'])
    assert len(unique)==93 and len(controls)==45 and len(pulse_errors)==48 and len(native_cases)==45
    for f in p['factorials']:
        for P in ('S','R'):
            ts=[next(t for t in p['trials'] if t['context']==f['contexts'][C+'_'+P] and t['target']==f['cell']) for C in ('S','R')]
            z=prepared[ts[0]['context']];i=z['ids'].tolist().index(f['cell'])
            assert np.array_equal(starts[ts[0]['key']][:,i],starts[ts[1]['key']][:,i])
            assert ts[0]['requested_removed_amount']==ts[1]['requested_removed_amount']
    for f,h in list(WATCHED.items()):assert digest(f)==h,f
    return dict(passed=True,protocol_sha256=digest(root/'protocol.json'),contexts=45,logical_cases=93,unique_initial_value_problems=93,controls=45,pulses=48,native_cases=45,cases=cases,matching_accounting=reset_details,physical_checks=physical_checks,maximum_pulse_relative_amount_error=max(pulse_errors),inputs_unchanged=True,inputs_sha256=WATCHED,script_sha256=digest(__file__))

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--root',type=Path,default=Path('output/studies/phase-carry-recipient-state-response-grouped'));a.add_argument('--output',type=Path,default=Path('archive/study-launches/2026-10-08-recipient-state-response-grouped/input-audit.json'));args=a.parse_args()
    r=audit(args.root);args.output.write_text(json.dumps(r,indent=2)+'\n');print({k:r[k] for k in ('passed','contexts','logical_cases','unique_initial_value_problems','native_cases','maximum_pulse_relative_amount_error')})
