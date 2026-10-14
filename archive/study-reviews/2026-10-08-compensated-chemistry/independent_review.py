"""Reconstruct backend qualification from saved arrays; no project imports."""
import importlib.util
import json
from pathlib import Path
import numpy as np

HELPER=Path('archive/study-reviews/2026-10-08-recipient-amount-delayed-response/independent_review.py')
spec=importlib.util.spec_from_file_location('array_review',HELPER)
r=importlib.util.module_from_spec(spec);spec.loader.exec_module(r)
ROOT=Path('output/studies/compensated-chemistry-validation')
OLD=Path('output/studies/phase-carry-recipient-amount-delayed-response')
METHOD='float64 SSP-RK2 increment with Kahan accumulation v1'


def correction(state,carry):
    assert carry.dtype==np.float64 and carry.shape==state.shape and np.isfinite(carry).all()
    assert np.all(abs(carry)<=2*np.spacing(state))


def layout(p):
    contexts={c['key']:r.read_npz(c['source']) for c in p['contexts']}
    groups={};mapping={}
    for t in p['trials']:
        z=contexts[t['context']];sig=tuple((z[k].dtype.str,z[k].shape,z[k].tobytes()) for k in ('delta','volumes','ids'))
        groups.setdefault(sig,[])
        mapping[t['key']]=(list(groups).index(sig),len(groups[sig]));groups[sig].append(t['key'])
    return contexts,groups,mapping


def production(folder,p,times):
    contexts,groups,mapping=layout(p);width=max(map(len,groups.values()));n=16
    paths={};rows=[]
    ph=r.digest(folder/'protocol.json')
    for level in p['levels']:
        label=level['level'];sub=folder/'production'/label
        report=r.read_json(sub/'result.json');z=r.read_npz(sub/'paths.npz');cp=r.read_npz(sub/'checkpoint.npz');commit=r.read_json(sub/'checkpoint-commit.json')
        assert report['protocol_sha256']==ph and report['method']==METHOD and report['level']==level
        assert report['paths_sha256']==r.digest(sub/'paths.npz') and report['checkpoint_sha256']==r.digest(sub/'checkpoint.npz')
        assert commit['checkpoint_sha256']==report['checkpoint_sha256'] and commit['method']==METHOD and commit['protocol_sha256']==ph
        assert np.array_equal(z['times'],times) and list(z['trial_keys'])==[t['key'] for t in p['trials']]
        r.valid(z['trajectory'],(len(times),len(p['trials']),2,16));r.valid(cp['state'],(2,16*len(groups),width))
        correction(cp['state'],cp['correction']);assert np.array_equal(z['trajectory'],cp['history'])
        meta=json.loads(str(cp['metadata']));assert meta['method']==METHOD and meta['protocol_sha256']==ph and meta['level']==level
        assert meta['step']==commit['step']==round(p['duration']/level['dt'])
        for i,t in enumerate(p['trials']):
            key=t['key'];gi,col=mapping[key];path=z['trajectory'][:,i]
            assert np.array_equal(path[0],r.read_npz(t['source'])['initial'])
            assert np.array_equal(path[-1],cp['state'][:,16*gi:16*(gi+1),col]);paths[label,key]=path
        for gi,group in enumerate(groups.values()):
            for col in range(len(group),width):
                assert np.array_equal(cp['state'][:,16*gi:16*(gi+1),col],cp['state'][:,16*gi:16*(gi+1),0])
                assert np.array_equal(cp['correction'][:,16*gi:16*(gi+1),col],cp['correction'][:,16*gi:16*(gi+1),0])
        assert report['actual_paths']==len(p['trials']) and report['padding_paths']==len(groups)*width-len(p['trials'])
        rows.append(dict(level=label,method_verified=True,correction_checkpoint_verified=True,positive_saved_states=True))
    return paths,rows


def main():
    r.pin(HELPER);r.pin(__file__)
    p=r.read_json(ROOT/'protocol.json');old=r.read_json(OLD/'protocol.json');summary=r.read_json(ROOT/'summary.json');status=r.read_json(ROOT/'status.json')
    for field in ('source_sha256','input_sha256'):
        for file,h in p[field].items():r.pin(file,h)
    assert p['solver']['method']==METHOD and p['criteria']==old['criteria'] and p['levels']==old['levels']
    selected={c['key'] for c in old['contexts'] if c['matching_reference']=='S' and c['physics_label']=='S'}
    assert p['trials']==[t for t in old['trials'] if t['context'] in selected]
    assert len(p['trials'])==36 and len(p['contexts'])==12
    pref=r.read_json(ROOT/'prefix/protocol.json');assert pref['trials']==old['trials'] and pref['contexts']==old['contexts']
    assert pref['duration']==.6 and pref['checkpoint_interval']==.3
    prepaths,prestages=production(ROOT/'prefix',pref,np.arange(5)*.15);preport=r.read_json(ROOT/'prefix/summary.json')
    prechecks=[]
    for level in p['levels']:
        label=level['level'];z=r.read_npz(ROOT/'prefix/uninterrupted'/label/'paths.npz');cp=r.read_npz(ROOT/'prefix/production'/label/'checkpoint.npz')
        gpu=np.stack([prepaths[label,t['key']] for t in pref['trials']],axis=1)
        assert np.array_equal(z['trajectory'],gpu) and np.array_equal(z['state'],cp['state']) and np.array_equal(z['correction'],cp['correction'])
        saved=r.read_json(ROOT/'prefix/production'/label/'result.json');assert saved['resumed_from_step']==round(.3/level['dt'])
        cpu=r.read_npz(ROOT/'prefix/cpu-prefix'/label/'paths.npz')['trajectory'];r.valid(cpu,gpu.shape)
        assert np.array_equal(cpu[0],gpu[0]);raw=r.log_error(gpu,cpu);norm=[]
        for i,t in enumerate(pref['trials']):
            if t['target'] is None:continue
            j=next(j for j,c in enumerate(pref['trials']) if c['context']==t['context'] and c['target'] is None)
            norm.append(float(abs(np.log(gpu[:,i]/gpu[:,j])-np.log(cpu[:,i]/cpu[:,j])).max()/abs(np.log(t['factor']))))
        previous=next(x for x in preport['levels'] if x['level']==level)
        r.compare_record(raw,previous['raw_log_error'],label+'/prefix/raw');r.compare_record(max(norm),previous['normalized_response_error'],label+'/prefix/response')
        passed=raw<=1e-8 and max(norm)<=.01;r.compare_record(passed,previous['passed'],label+'/prefix/pass')
        prechecks.append(dict(level=label,raw_log_error=raw,normalized_response_error=max(norm),gpu_restart_bitwise_identical=True,passed=bool(passed)))
    times=np.arange(401)*.15;paths,stages=production(ROOT,p,times);contexts,_,_=layout(p)
    noise={};accuracy=[];orders=[];rawpairs=[]
    for t in p['trials']:
        key=t['key'];sub=OLD/'references'/key;z=r.read_npz(sub/'paths.npz');report=r.read_json(sub/'result.json')
        assert report['trial']==t and report['passed'] and report['protocol_sha256']==r.digest(OLD/'protocol.json') and report['paths_sha256']==r.digest(sub/'paths.npz')
        assert np.array_equal(z['times'],times)
        for label in ('dop853','radau'):
            r.valid(z[label],(401,2,16));assert np.array_equal(z[label][0],paths['fine',key][0]);paths[label,key]=z[label]
        noise[key]=r.log_error(z['dop853'],z['radau']);errors=[]
        for level in p['levels']:
            label=level['level'];sub=ROOT/'cpu-verification'/label/key;z=r.read_npz(sub/'paths.npz');report=r.read_json(sub/'result.json')
            assert report['trial']==t and report['method']==METHOD and report['level']==level and report['protocol_sha256']==r.digest(ROOT/'protocol.json') and report['paths_sha256']==r.digest(sub/'paths.npz')
            assert np.array_equal(z['times'],times);r.valid(z['trajectory'],(401,2,16));assert np.array_equal(z['trajectory'][0],paths[label,key][0]);correction(z['trajectory'][-1],z['final_correction'])
            error=max(r.log_error(paths[label,key],paths[ref,key]) for ref in ('dop853','radau'));cpu=r.log_error(paths[label,key],z['trajectory']);errors.append(error)
            previous=next(x for x in summary['accuracy'] if (x['trial'],x['level'])==(key,label))
            row=dict(trial=key,level=label,raw_reference_log_error=error,cpu_torch_log_error=cpu,passed=bool(error<=.01 and cpu<=1e-8));accuracy.append(row)
            for name in ('raw_reference_log_error','cpu_torch_log_error','passed'):r.compare_record(row[name],previous[name],key+'/'+label+'/'+name)
        for a,b in [('coarse','fine'),('fine','diagnostic')]:
            error=r.log_error(paths[a,key],paths[b,key]);passed=error<=.01;rawpairs.append(dict(trial=key,first=a,second=b,raw_log_error=error,passed=bool(passed)))
            previous=next(x for x in summary['raw_pairs'] if (x['trial'],x['first'],x['second'])==(key,a,b))
            r.compare_record(error,previous['raw_log_error'],key+'/pair');r.compare_record(passed,previous['passed'],key+'/pair/pass')
        resolved=errors[-1]>100*noise[key];values=np.log2(np.array(errors[:-1])/errors[1:]).tolist() if resolved and all(x>0 for x in errors) else None
        row=dict(trial=key,errors=errors,reference_noise=noise[key],order_resolved=bool(resolved),orders=values,passed=bool(not resolved or values is not None and all(1.7<=x<=2.3 for x in values)));orders.append(row)
        previous=next(x for x in summary['orders'] if x['trial']==key)
        for name in ('reference_noise','order_resolved','passed'):r.compare_record(row[name],previous[name],key+'/order/'+name)
        for x,y in zip(errors,previous['errors']):assert x==y
        if values is None:assert previous['orders'] is None
        else:
            for x,y in zip(values,previous['orders']):r.compare_record(x,y,key+'/order')
    controls={t['context']:t['key'] for t in p['trials'] if t['target'] is None};changes={};metrics={}
    for label in ('dop853','radau','coarse','fine','diagnostic'):
        for t in p['trials']:
            if t['target'] is None:continue
            z=contexts[t['context']];i=z['ids'].tolist().index(t['target'])
            _,changes[label,t['key']],metrics[label,t['key']]=r.response(paths[label,t['key']],paths[label,controls[t['context']]],times,z['volumes'],i,t['factor'],24.)
    comparisons=[]
    for saved in summary['response_comparisons']:
        a,b=saved['first'],saved['second'];limit=1e-8 if a=='dop853' else .01;rows=[]
        assert saved['effects']==[]
        for t in p['trials']:
            if t['target'] is None:continue
            key=t['key'];x,y=metrics[a,key],metrics[b,key];error=float(abs(changes[a,key]-changes[b,key]).max()/abs(np.log(t['factor'])))
            auc=abs(y['target_activator_log_auc_per_log_pulse']/x['target_activator_log_auc_per_log_pulse']-1.);recovery={}
            for name in ('target_recovery_time','network_recovery_time'):
                u,v=x[name],y[name];recovery[name]=0. if u is None and v is None else None if u is None or v is None else abs(u-v)
            row=dict(trial=key,normalized_response_error=error,relative_auc_error=float(auc),recovery_time_errors=recovery,passed=bool(error<=limit and auc<=.02 and all(v is not None and v<=.3+1e-9 for v in recovery.values())));rows.append(row)
            previous=next(z for z in saved['pulse_rows'] if z['trial']==key)
            for name in ('normalized_response_error','relative_auc_error','passed'):r.compare_record(row[name],previous[name],key+'/response/'+name)
            for name,value in recovery.items():r.compare_record(value,previous['recovery_time_errors'][name],key+'/'+name)
        passed=all(x['passed'] for x in rows);r.compare_record(passed,saved['passed'],a+'/'+b)
        comparisons.append(dict(first=a,second=b,pulse_rows=rows,passed=passed))
    failed=next(x for x in r.read_json(OLD/'summary.json')['orders'] if not x['passed'])['trial'];row=next(x for x in orders if x['trial']==failed)
    drift=r.log_error(paths['fine',failed],paths['fine',failed][0]);diagnosis=not row['order_resolved'] and row['errors'][-1]<=100*row['reference_noise'] and drift>0
    r.compare_record(drift,summary['diagnosed_control']['stored_max_log_relaxation'],'diagnosed relaxation');r.compare_record(diagnosis,summary['diagnosed_control']['passed'],'diagnosis')
    for file,h in summary['evidence_sha256'].items():r.pin(file,h)
    guards=r.read_json('archive/study-document-snapshots/2026-10-08-before-compensated-chemistry/manifest.json')['unchanged_guards']
    for file,h in guards.items():r.pin(file,h)
    passed=bool(diagnosis and all(x<=1e-8 for x in noise.values()) and all(x['passed'] for x in accuracy+orders+rawpairs+comparisons+prechecks))
    r.compare_record(passed,summary['passed'],'overall');r.compare_record(passed,status['passed'],'status');assert status['state']==('completed' if passed else 'completed_with_unresolved_checks')
    for file,h in r.WATCHED.items():assert r.digest(file)==h
    result=dict(passed=passed,method=METHOD,protocol_sha256=r.digest(ROOT/'protocol.json'),counts=summary['counts'],accuracy=accuracy,orders=orders,raw_pairs=rawpairs,response_comparisons=comparisons,prefix=prechecks,checkpoint_checks=stages+prestages,diagnosed_control=dict(order_record=row,stored_max_log_relaxation=drift,passed=bool(diagnosis)),inputs_unchanged=True,inputs_sha256=r.WATCHED,unchanged_guards=len(guards),script_sha256=r.digest(__file__))
    out=Path(__file__).with_suffix('.json');out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(passed=passed,counts=result['counts'],diagnosed_control=result['diagnosed_control'],raw_error_max=max(x['raw_reference_log_error'] for x in accuracy),cpu_gpu_error_max=max(x['cpu_torch_log_error'] for x in accuracy)),indent=2))

if __name__=='__main__':main()
