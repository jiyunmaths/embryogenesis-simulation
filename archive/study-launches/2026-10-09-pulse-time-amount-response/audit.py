"""Literal input and actual-start array audit; no project helpers or stepping."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np

ROOT=Path('output/studies/pulse-time-amount-response')
PINS={}


def sha(file):
    h=hashlib.sha256()
    with Path(file).open('rb') as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()


def pin(file,expected=None):
    name=str(Path(file).resolve());value=PINS.get(name)
    if value is None:value=sha(file);PINS[name]=value
    assert expected is None or value==expected,name
    return Path(name)


def read(file):return json.loads(pin(file).read_text())


def arrays(file,keys):
    with np.load(pin(file),allow_pickle=False) as z:return {k:z[k].copy() for k in keys}



def inputs():
    p=read(ROOT/'protocol.json')
    for collection in ('source_sha256','input_sha256'):
        for file,expected in p[collection].items():pin(file,expected)
    assert (len(p['contexts']),len(p['trials']),len(p['factorials']))==(168,336,42)
    assert p['start']==570. and p['duration']==60. and p['histories']==[7,8,9] and p['new_histories']==0
    assert p['criteria']['pulse_amount_relative_max']==p['amount_relative_max']==2e-13
    parent=read(Path(p['parent'])/'protocol.json');review=read('docs/pulse_time_recipient_response_assessment.json')
    assert review['full_study_accepted'] and review['protocol_sha256']==sha(Path(p['parent'])/'protocol.json')
    pin(review['independent_review'],review['independent_review_sha256']);audit=read(review['independent_review'])
    assert audit['passed'] and audit['accepted_at_tested_scope']
    assert p['parent_historical_export_note']==review['historical_guard_audit']
    saved={};states={};largest=0.;amount_error=0.;fraction_error=0.
    for c in parent['contexts']:
        pin(c['source'],audit['evidence_sha256'][str(Path(c['source']).resolve())])
        saved[c['key']]=arrays(c['source'],('initial','ids','volumes','delta'))
    signatures=set()
    for c in p['contexts']:
        z=arrays(c['source'],('initial','ids','volumes','delta'));states[c['key']]=z
        original=saved[c['parent_context']];expected=original['initial'].copy()
        for k in ('ids','volumes','delta'):np.testing.assert_array_equal(z[k],original[k])
        i=z['ids'].tolist().index(c['targets'][0]);others=z['ids']!=c['targets'][0];B=c['budget_reference']
        assert c['matching_reference']!='native' or B=='native'
        if B!='native':
            f=next(f for f in parent['factorials'] if (f['seed'],f['cell'],f['matching_reference'])==(c['seed'],c['targets'][0],c['matching_reference']))
            ref=saved[f['contexts'][B+'_'+c['physical_origin']]]
            before=expected[:,others]@z['volumes'][others];common=ref['initial'][:,others]@z['volumes'][others];scale=common/before
            expected[:,others]*=scale[:,None];after=expected[:,others]@z['volumes'][others]
            a=original['initial'][:,others]*z['volumes'][others]/before[:,None];b=expected[:,others]*z['volumes'][others]/after[:,None]
            ae=float(abs(after/common-1).max());fe=float(abs(b/a-1).max());amount_error=max(amount_error,ae);fraction_error=max(fraction_error,fe)
            assert max(ae,fe)<=2e-13
            ac=c['amount_matching']
            values=dict(target_pair=expected[:,i],surrounding_indices=np.flatnonzero(others),pre_surrounding_species_amounts=before,
                common_surrounding_species_amounts=common,scale_factors=scale,post_surrounding_species_amounts=after,
                pre_global_species_amounts=original['initial']@z['volumes'],post_global_species_amounts=expected@z['volumes'],
                pre_per_cell_species_amounts=original['initial']*z['volumes'],post_per_cell_species_amounts=expected*z['volumes'],
                per_cell_species_amount_change=(expected-original['initial'])*z['volumes'],
                pre_surrounding_allocation_fractions=a,post_surrounding_allocation_fractions=b,
                maximum_relative_amount_error=ae,maximum_relative_allocation_error=fe)
            for k,v in values.items():np.testing.assert_array_equal(ac[k],v)
        else:assert c['amount_matching'] is None
        np.testing.assert_array_equal(expected,z['initial']);np.testing.assert_array_equal(z['initial'][:,i],original['initial'][:,i])
        np.testing.assert_array_equal(c['post_reset_global_species_amounts'],expected@z['volumes'])
        v=z['volumes'];d=z['delta'];off=d.copy();np.fill_diagonal(off,0.)
        assert np.all(v>0) and len(set(z['ids'].tolist()))==16 and np.all(off>=0)
        assert abs(d.sum(1)).max()<1e-11 and abs(v@d).max()<1e-12 and abs(v[:,None]*d-(v[:,None]*d).T).max()<1e-12
        assert expected.dtype==np.float64 and np.isfinite(expected).all() and np.all(expected>0)
    trial_values={}
    for t in p['trials']:
        z=states[t['context']];values=arrays(t['source'],('initial',))['initial'];expected=z['initial'].copy()
        if t['target'] is not None:
            i=z['ids'].tolist().index(t['target']);expected[0,i]*=.9
            amount=(z['initial'][0,i]-expected[0,i])*z['volumes'][i];wanted=.1*z['initial'][0,i]*z['volumes'][i]
            assert amount==t['removed_activator_amount'] and abs(amount/wanted-1)<=2e-13 and t['factor']==.9
            assert t['initial_log_amplitude']==abs(np.log(.9));largest=max(largest,abs(amount/wanted-1))
        np.testing.assert_array_equal(values,expected);trial_values[t['key']]=values
        signatures.add(tuple((z[k].dtype.str,z[k].shape,z[k].tobytes()) for k in ('delta','volumes','ids'))+(values.tobytes(),))
    matched_pairs=0;budget_pairs=0;paired_budget_error=0.
    for f in p['factorials']:
        for P in ('S','R'):
            a,b=[states[f['contexts'][C+'_'+P]] for C in ('S','R')]
            for k in ('ids','volumes','delta'):np.testing.assert_array_equal(a[k],b[k])
            if f['matching_reference']=='native':continue
            i=a['ids'].tolist().index(f['cell']);np.testing.assert_array_equal(a['initial'][:,i],b['initial'][:,i]);matched_pairs+=1
            pair=[next(t for t in p['trials'] if t['context']==f['contexts'][C+'_'+P] and t['target'] is not None) for C in ('S','R')]
            assert pair[0]['removed_activator_amount']==pair[1]['removed_activator_amount']
            if f['budget_reference']!='native':
                others=a['ids']!=f['cell'];v=a['volumes'];budget_pairs+=1
                for av,bv in [(a['initial']@v,b['initial']@v),(a['initial'][:,others]@v[others],b['initial'][:,others]@v[others]),(trial_values[pair[0]['key']]@v,trial_values[pair[1]['key']]@v)]:
                    error=float(abs(bv/av-1).max());assert error<=2e-13;paired_budget_error=max(paired_budget_error,error)
    assert matched_pairs==72 and budget_pairs==48
    assert len(signatures)==p['unique_physical_initial_value_problems']
    return p,dict(contexts=168,logical_starts=336,paired_matching_checks=matched_pairs,paired_budget_checks=budget_pairs,
        distinct_physical_problems=len(signatures),max_pulse_amount_relative_error=largest,max_budget_relative_error=amount_error,
        max_allocation_fraction_relative_error=fraction_error,max_paired_budget_relative_error=paired_budget_error)


def main(action):
    p,details=inputs()
    if action=='prepare':
        result=dict(passed=True,protocol_sha256=sha(ROOT/'protocol.json'),independent_histories=3,new_histories=0,
            parent_noninput_export_exception_preserved=True,**details)
        file=Path(__file__).with_name('preparation-audit.json')
    else:
        pref=ROOT/'prefix';q=read(pref/'protocol.json');summary=read(pref/'summary.json')
        assert summary['passed'] and summary['protocol_sha256']==sha(pref/'protocol.json') and q['trials']==p['trials']
        rows=[]
        for level in p['levels']:
            name=level['level'];dt=level['dt'];prod=pref/'production'/name
            gpu=arrays(prod/'paths.npz',('trajectory','times','trial_keys'));cpu=arrays(pref/'cpu-prefix'/name/'paths.npz',('trajectory',))
            shadow=arrays(pref/'uninterrupted'/name/'paths.npz',('trajectory','state','correction'))
            cp=arrays(prod/'checkpoint.npz',('state','correction','history','metadata'))
            commit=read(prod/'checkpoint-commit.json');pin(prod/'checkpoint.npz',commit['checkpoint_sha256'])
            record=read(prod/'result.json');pin(prod/'paths.npz',record['paths_sha256'])
            assert record['method']==summary['method'] and record['protocol_sha256']==sha(pref/'protocol.json')
            assert cp['state'].dtype==cp['correction'].dtype==np.float64 and np.isfinite(cp['correction']).all()
            assert np.all(abs(cp['correction'])<=2*np.spacing(cp['state']))
            np.testing.assert_array_equal(gpu['trajectory'],shadow['trajectory'])
            np.testing.assert_array_equal(cp['history'],gpu['trajectory'])
            np.testing.assert_array_equal(cp['state'],shadow['state']);np.testing.assert_array_equal(cp['correction'],shadow['correction'])
            np.testing.assert_array_equal(gpu['times'],np.arange(5)*.15)
            np.testing.assert_array_equal(gpu['trial_keys'],[t['key'] for t in p['trials']])
            meta=json.loads(str(cp['metadata']));assert meta['step']==round(.6/dt) and meta['protocol_sha256']==sha(pref/'protocol.json')
            raw=float(abs(np.log(gpu['trajectory']/cpu['trajectory'])).max());normalized=[]
            for i,t in enumerate(p['trials']):
                initial=arrays(t['source'],('initial',))['initial'];np.testing.assert_array_equal(gpu['trajectory'][0,i],initial)
                if t['target'] is None:continue
                j=next(j for j,c in enumerate(p['trials']) if c['context']==t['context'] and c['target'] is None)
                a,b=gpu['trajectory'],cpu['trajectory']
                normalized.append(float(abs(np.log(a[:,i]/a[:,j])-np.log(b[:,i]/b[:,j])).max()/t['initial_log_amplitude']))
            old=next(r for r in summary['levels'] if r['level']==level)
            assert raw==old['raw_log_error'] and max(normalized)==old['normalized_response_error']
            assert raw<=1e-8 and max(normalized)<=.01
            rows.append(dict(level=name,logical_starts=336,raw_log_error=raw,normalized_response_error=max(normalized),restart_bitwise_identical=True,passed=True))
        result=dict(passed=True,protocol_sha256=sha(pref/'protocol.json'),summary_sha256=sha(pref/'summary.json'),
            actual_context_timestep_cases=1008,rows=rows,**details)
        file=pref/'independent-audit.json'
    for path,expected in PINS.items():assert sha(path)==expected,path
    result.update(evidence_sha256=PINS,script_sha256=sha(__file__))
    file.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('evidence_sha256',)},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('prepare','prefix'));main(parser.parse_args().action)
