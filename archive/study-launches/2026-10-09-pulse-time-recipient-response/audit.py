"""Literal input and actual-start array audit; no project helpers or stepping."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np

ROOT=Path('output/studies/pulse-time-recipient-response')
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
    assert (len(p['contexts']),len(p['trials']),len(p['factorials']))==(72,144,18)
    assert p['start']==570. and p['duration']==60. and p['histories']==[7,8,9] and p['new_histories']==0
    parent=read(Path(p['parent'])/'protocol.json');review=read('docs/moving_recipient_amount_response_assessment.json')
    assert review['full_study_accepted'] and review['protocol_sha256']==sha(Path(p['parent'])/'protocol.json')
    audit=read(review['independent_review']);pin(review['independent_review'],review['independent_review_sha256'])
    assert audit['passed'] and audit['accepted_at_tested_scope']
    assert p['parent_historical_export_note']==review['historical_guard_audit']
    endpoints={};states={};largest=0.
    for c in parent['contexts']:
        job=next(j for j in parent['jobs'] if (j['context'],j['stage'],j['level'])==(c['key'],'wait','fine'))
        folder=Path(p['parent'])/job['key']
        for name in ('history.json','result.json','latest_state.npz'):
            file=folder/name;pin(file,audit['evidence_sha256'][str(file.resolve())])
        h=read(folder/'history.json');cp=arrays(folder/'latest_state.npz',('activator','inhibitor','ids'))
        u=np.array([cp['activator'],cp['inhibitor']]);np.testing.assert_array_equal(u,h[-1]['chemistry'])
        assert h[-1]['time']==570. and h[-1]['elapsed']==60.
        endpoints[c['key']]=dict(initial=u,ids=cp['ids'],volumes=np.asarray(h[-1]['volumes']),delta=np.asarray(h[-1]['delta']))
    for c in p['contexts']:
        z=arrays(c['source'],('initial','ids','volumes','delta'));states[c['key']]=z
        chemical=endpoints[c['chemical_parent_context']];physical=endpoints[c['physical_parent_context']]
        np.testing.assert_array_equal(z['ids'],physical['ids']);np.testing.assert_array_equal(z['ids'],chemical['ids'])
        np.testing.assert_array_equal(z['volumes'],physical['volumes']);np.testing.assert_array_equal(z['delta'],physical['delta'])
        expected=chemical['initial'].copy();i=z['ids'].tolist().index(c['targets'][0]);M=c['matching_reference']
        if M!='native':
            pair=next(pair for pair in parent['pairs'] if (pair['seed'],pair['cell'])==(c['seed'],c['targets'][0]))
            expected[:,i]=endpoints[pair['contexts'][M]]['initial'][:,i]
        np.testing.assert_array_equal(z['initial'],expected)
        v=z['volumes'];d=z['delta'];assert expected.dtype==np.float64 and np.isfinite(expected).all() and np.all(expected>0)
        assert abs(d.sum(1)).max()<1e-11 and abs(v@d).max()<1e-12 and abs(v[:,None]*d-(v[:,None]*d).T).max()<1e-12
        original=chemical['initial']*chemical['volumes'];reassigned=chemical['initial']*v;post=expected*v
        np.testing.assert_array_equal(c['reassignment_per_cell_amount_change'],reassigned-original)
        np.testing.assert_array_equal(c['reassignment_global_amount_change'],(reassigned-original).sum(axis=1))
        np.testing.assert_array_equal(c['post_reset_global_species_amounts'],expected@v)
        if M!='native':
            account=c['matching_accounting']
            np.testing.assert_array_equal(account['matching_per_cell_species_amount_change'],post-reassigned)
            np.testing.assert_array_equal(account['matching_global_species_amount_change'],(post-reassigned).sum(axis=1))
    for t in p['trials']:
        z=states[t['context']];values=arrays(t['source'],('initial',))['initial'];expected=z['initial'].copy()
        if t['target'] is not None:
            i=z['ids'].tolist().index(t['target']);expected[0,i]*=.9
            amount=(z['initial'][0,i]-expected[0,i])*z['volumes'][i];wanted=.1*z['initial'][0,i]*z['volumes'][i]
            assert amount==t['removed_activator_amount'] and abs(amount/wanted-1)<=2e-13
            largest=max(largest,abs(amount/wanted-1));assert t['factor']==.9
        np.testing.assert_array_equal(values,expected)
    matched_pairs=0
    for f in p['factorials']:
        for P in ('S','R'):
            a,b=[states[f['contexts'][C+'_'+P]] for C in ('S','R')]
            np.testing.assert_array_equal(a['volumes'],b['volumes']);np.testing.assert_array_equal(a['delta'],b['delta'])
            if f['matching_reference']!='native':
                i=a['ids'].tolist().index(f['cell']);np.testing.assert_array_equal(a['initial'][:,i],b['initial'][:,i])
                pair=[next(t for t in p['trials'] if t['context']==f['contexts'][C+'_'+P] and t['target'] is not None) for C in ('S','R')]
                assert pair[0]['removed_activator_amount']==pair[1]['removed_activator_amount'];matched_pairs+=1
    return p,dict(contexts=72,logical_starts=144,paired_matching_checks=matched_pairs,max_pulse_amount_relative_error=largest)


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
            rows.append(dict(level=name,logical_starts=144,raw_log_error=raw,normalized_response_error=max(normalized),restart_bitwise_identical=True,passed=True))
        result=dict(passed=True,protocol_sha256=sha(pref/'protocol.json'),summary_sha256=sha(pref/'summary.json'),
            actual_context_timestep_cases=432,rows=rows,**details)
        file=pref/'independent-audit.json'
    for path,expected in PINS.items():assert sha(path)==expected,path
    result.update(evidence_sha256=PINS,script_sha256=sha(__file__))
    file.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('evidence_sha256',)},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('prepare','prefix'));main(parser.parse_args().action)
