"""Separate array reconstruction: no project helper imports or stepping."""
import argparse
import hashlib
import json
from pathlib import Path
from collections import defaultdict
import numpy as np
from scipy.integrate import trapezoid

WATCHED = {}
RECORDED_DIFFERENCES = []


def digest(file):
    h=hashlib.sha256()
    with Path(file).open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''): h.update(block)
    return h.hexdigest()


def pin(file, expected=None):
    name=str(Path(file).resolve()); actual=digest(name)
    if expected is not None and actual != expected: raise AssertionError('Input hash mismatch: '+name)
    if name in WATCHED and WATCHED[name] != actual: raise AssertionError('Input changed during review: '+name)
    WATCHED[name]=actual
    return name


def read_json(file):
    name=pin(file)
    return json.loads(Path(name).read_text())


def read_npz(file, keys=None):
    name=pin(file)
    with np.load(name,allow_pickle=False) as z:
        return {k:z[k].copy() for k in (z.files if keys is None else keys)}


def compare_record(value, recorded, label):
    if value is None or recorded is None:
        if value is not None or recorded is not None: raise AssertionError('Recorded null mismatch: '+label)
    elif isinstance(value,(bool,np.bool_)):
        if bool(value) != bool(recorded): raise AssertionError('Recorded decision mismatch: '+label)
    else:
        error=abs(float(value)-float(recorded))
        RECORDED_DIFFERENCES.append((error,label))
        if not np.isclose(value,recorded,rtol=1e-9,atol=2e-12): raise AssertionError(f'Recorded value mismatch: {label}: {value} vs {recorded}')


def valid(path, shape=None):
    if path.dtype != np.float64 or (shape is not None and path.shape != shape) or not np.isfinite(path).all() or np.any(path<=0):
        raise AssertionError('Invalid saved concentration trajectory')


def log_error(a,b): return float(np.max(np.abs(np.log(a/b))))


def recovery(times,deviation,threshold,window):
    future=np.logical_and.accumulate((deviation<=threshold)[::-1])[::-1]
    indices=np.flatnonzero(future & (times<=times[-1]-window))
    return float(times[indices[0]]) if len(indices) else None


def response(path,control,times,volumes,index,factor,window):
    amplitude=float(abs(np.log(factor))); change=np.log(path/control)
    wave=change[:,:,index]/amplitude
    target=np.max(np.abs(change[:,:,index]),axis=1)
    local=np.abs(change[:,0,index])
    other=np.arange(len(volumes))!=index
    neighbors=np.sqrt(np.sum(change[:,:,other]**2*volumes[other],axis=(1,2))/(2*sum(volumes[other])))
    separate=np.log(path)-np.log(control)
    network=np.sqrt(np.sum(separate**2*volumes,axis=(1,2))/(2*sum(volumes)))
    metrics=dict(target_activator_peak_gain=float(local.max()/amplitude),
        target_inhibitor_peak_gain=float(np.abs(change[:,1,index]).max()/amplitude),
        target_activator_log_auc_per_log_pulse=float(trapezoid(local,times)/amplitude),
        other_cells_peak_log_rms_per_log_pulse=float(neighbors.max()/amplitude),
        target_recovery_time=recovery(times,target,.1*amplitude,window),
        network_recovery_time=recovery(times,network,.1*network[0],window),
        final_log_rms_from_control=float(network[-1]),
        late_max_log_rms_from_control=float(network[times>=times[-1]-window].max()))
    return wave,change,metrics


def rms(wave,times):
    return float(np.sqrt(trapezoid(np.mean(np.asarray(wave)**2,axis=1),times)/(times[-1]-times[0])))



def audit(root):
    ROOT=Path.cwd();STUDY=Path(root).resolve()
    p=read_json(STUDY/'protocol.json');old=read_json(Path(p['parent'])/'protocol.json')
    watched={**p['source_sha256'],**p['input_sha256'],str(STUDY/'protocol.json'):digest(STUDY/'protocol.json')}
    for name,expected in watched.items():assert digest(name)==expected,name
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
        unique.add(signature);cases.append(dict(trial=t['key'],source_sha256=digest(t['source']),exact_preparation=True))
    assert (len(contexts),len(cases),sum(t['target'] is None for t in p['trials']))==(45,141,45)
    assert len(unique)==p['unique_initial_value_problems']
    for name,expected in watched.items():assert digest(name)==expected,name
    result=dict(passed=True,protocol_sha256=digest(STUDY/'protocol.json'),contexts=45,logical_cases=141,
        controls=45,pulses=96,unique_initial_value_problems=len(unique),source_checks=sources,
        source_gate_count=45,descriptively_settled_backgrounds=settled,paired_context_checks=physical,cases=cases,
        maximum_pulse_relative_amount_error=maximum_amount_error,
        inputs_unchanged=True,inputs_sha256=watched,script_sha256=digest(__file__),
        scope='Actual unpulsed t570 states preserved exactly, no rematching; both delayed doses checked; whole-network amounts and target states unmatched.')
    for name,expected in watched.items():pin(name,expected)
    return result

def main(root, output):
    root=Path(root).resolve();output=Path(output)
    prepared=audit(root)
    p=read_json(root/'protocol.json');ph=digest(root/'protocol.json')
    status=read_json(root/'status.json');summary=read_json(root/'summary.json')
    assert status['state'] in ('completed','completed_with_unresolved_checks')
    assert summary['protocol_sha256']==ph
    expected=dict(raw_log_max=.01,normalized_response_max=.01,relative_auc_max=.02,recovery_time_max=.3,effect_error_max=.01,reference_log_max=1e-8,reference_normalized_max=1e-8,cpu_torch_log_max=1e-8,order_min=1.7,order_max=2.3,reference_noise_factor=100.,pulse_amount_relative_max=2e-13)
    assert p['criteria']==expected
    times=np.arange(401)*.15;trials=p['trials'];contexts={c['key']:read_npz(c['source']) for c in p['contexts']}
    starts={t['key']:read_npz(t['source'])['initial'] for t in trials}
    controls={t['context']:t['key'] for t in trials if t['target'] is None}
    paths={};noise={};stages=[]
    for t in trials:
        folder=root/'references'/t['key'];r=read_json(folder/'result.json');z=read_npz(folder/'paths.npz')
        assert r['protocol_sha256']==ph and r['trial']==t and r['paths_sha256']==digest(folder/'paths.npz') and np.array_equal(z['times'],times)
        for label in ('dop853','radau'):
            valid(z[label],(401,2,16));assert np.array_equal(z[label][0],starts[t['key']]);paths[label,t['key']]=z[label]
        noise[t['key']]=log_error(z['dop853'],z['radau'])
        compare_record(noise[t['key']],r['raw_solver_log_error'],t['key']+'/reference')
    groups={};mapping={}
    for t in trials:
        z=contexts[t['context']];sig=tuple((z[k].dtype.str,z[k].shape,z[k].tobytes()) for k in ('delta','volumes','ids'))
        if sig not in groups:groups[sig]=[]
        mapping[t['key']]=(list(groups).index(sig),len(groups[sig]));groups[sig].append(t['key'])
    assert len(groups)==9 and max(map(len,groups.values()))==23
    for level in p['levels']:
        label=level['level'];folder=root/'production'/label;r=read_json(folder/'result.json');z=read_npz(folder/'paths.npz');cp=read_npz(folder/'checkpoint.npz');commit=read_json(folder/'checkpoint-commit.json')
        assert r['protocol_sha256']==ph and r['level']==level and r['paths_sha256']==digest(folder/'paths.npz') and r['checkpoint_sha256']==digest(folder/'checkpoint.npz')
        assert commit['protocol_sha256']==ph and commit['checkpoint_sha256']==r['checkpoint_sha256']
        assert list(z['trial_keys'])==[t['key'] for t in trials] and np.array_equal(z['times'],times)
        valid(z['trajectory'],(401,141,2,16));valid(cp['state'],(2,144,23));assert np.array_equal(z['trajectory'],cp['history'])
        meta=json.loads(str(cp['metadata']));assert meta['step']==round(60/level['dt']) and commit['step']==meta['step'] and meta['protocol_sha256']==ph and meta['level']==level
        for i,t in enumerate(trials):
            key=t['key'];gi,col=mapping[key];path=z['trajectory'][:,i]
            assert np.array_equal(path[0],starts[key]) and np.array_equal(path[-1],cp['state'][:,16*gi:16*(gi+1),col]);paths[label,key]=path
        padding_error=max([0.]+[log_error(cp['state'][:,16*gi:16*(gi+1),0],cp['state'][:,16*gi:16*(gi+1),col]) for gi,group in enumerate(groups.values()) for col in range(len(group),23)])
        assert r['actual_paths']==141 and r['padding_paths']==66
        stages.append(dict(level=label,logical_paths=141,padding=66,padding_error=padding_error,positive_saved_states=True,checkpoint_consistent=True,wall_seconds=r['wall_seconds']))
        for t in trials:
            folder=root/'cpu-verification'/label/t['key'];r=read_json(folder/'result.json');z=read_npz(folder/'paths.npz')
            assert r['protocol_sha256']==ph and r['trial']==t and r['level']==level and r['paths_sha256']==digest(folder/'paths.npz') and np.array_equal(z['times'],times)
            valid(z['trajectory'],(401,2,16));assert np.array_equal(z['trajectory'][0],starts[t['key']]);paths['cpu-'+label,t['key']]=z['trajectory']
    labels=('dop853','radau','coarse','fine','diagnostic');waves={};changes={};metrics={};contrasts={};norms={};diagonals={};identities=[]
    for label in labels:
        for t in trials:
            if t['target'] is None:continue
            z=contexts[t['context']];i=z['ids'].tolist().index(t['target']);key=t['key']
            waves[label,key],changes[label,key],metrics[label,key]=response(paths[label,key],paths[label,controls[t['context']]],times,z['volumes'],i,t['factor'],24.)
        for f in p['factorials']:
            seed,cell,L=f['seed'],f['cell'],f['matching_reference']
            for dose in ('fractional','amount'):
                four={}
                for name,context in f['contexts'].items():
                    t=next(t for t in trials if (t['context'],t['target'],t['dose'])==(context,cell,dose));four[name]=waves[label,t['key']]
                cc={P:four['R_'+P]-four['S_'+P] for P in ('S','R')};pp={C:four[C+'_R']-four[C+'_S'] for C in ('S','R')}
                interaction=four['R_R']-four['S_R']-four['R_S']+four['S_S'];diagonal=four['R_R']-four['S_S']
                kinds=[('chemical','sham_physics',cc['S']),('chemical','redistributed_physics',cc['R']),('physical','sham_chemistry',pp['S']),('physical','redistributed_chemistry',pp['R']),('interaction','difference_of_differences',interaction)]
                for kind,condition,wave in kinds:
                    key=f'{seed}_{cell}_{L}_{dose}_{kind}_{condition}';contrasts[label,key]=wave;norms[label,key]=rms(wave,times)
                diagonals[label,seed,cell,L,dose]=rms(diagonal,times)
                identities.extend(float(abs(x).max()) for x in [interaction-(cc['R']-cc['S']),interaction-(pp['R']-pp['S']),diagonal-cc['S']-pp['R'],diagonal-cc['R']-pp['S']])
        recorded=summary['response_summaries'][label]
        for t in trials:
            if t['target'] is None:continue
            for name,value in metrics[label,t['key']].items():compare_record(value,recorded['metrics'][t['key']][name],label+'/'+t['key']+'/'+name)
        assert len(recorded['effects'])==120
        for row in recorded['effects']:
            value=norms[label,row['key']];compare_record(value,row['normalized_response_shift_rms'],label+'/'+row['key'])
            if row['kind']=='chemical':compare_record(value>.01,row['response_shift_detected'],label+'/'+row['key']+'/decision')
            else:assert 'response_shift_detected' not in row
        for row in recorded['native_diagonal']:compare_record(diagonals[label,row['seed'],row['cell'],row['matching_reference'],row['dose']],row['normalized_response_shift_rms'],label+'/diagonal')
    specs=[('dop853','radau'),('coarse','fine'),('fine','diagnostic')]+[(l['level'],ref) for l in p['levels'] for ref in ('dop853','radau')]
    comparisons=[]
    for first,second in specs:
        limit=1e-8 if (first,second)==('dop853','radau') else .01;pulses=[];cr=[]
        saved=next(r for r in summary['response_comparisons'] if (r['first'],r['second'])==(first,second))
        for t in trials:
            if t['target'] is None:continue
            key=t['key'];a,b=metrics[first,key],metrics[second,key]
            error=float(abs(changes[first,key]-changes[second,key]).max()/abs(np.log(t['factor'])))
            area=abs(b['target_activator_log_auc_per_log_pulse']/a['target_activator_log_auc_per_log_pulse']-1.)
            recovery_errors={}
            for name in ('target_recovery_time','network_recovery_time'):
                x,y=a[name],b[name];recovery_errors[name]=0. if x is None and y is None else None if x is None or y is None else abs(x-y)
            passed=error<=limit and area<=.02 and all(x is not None and x<=.3+1e-9 for x in recovery_errors.values())
            row=dict(trial=key,normalized_response_error=error,relative_auc_error=float(area),recovery_time_errors=recovery_errors,passed=bool(passed));pulses.append(row)
            previous=next(r for r in saved['pulse_rows'] if r['trial']==key)
            for name in ('normalized_response_error','relative_auc_error','passed'):compare_record(row[name],previous[name],first+'/'+second+'/'+key+'/'+name)
            for name,x in recovery_errors.items():compare_record(x,previous['recovery_time_errors'][name],first+'/'+second+'/'+key+'/'+name)
        for row in summary['response_summaries'][first]['effects']:
            key=row['key'];error=float(abs(contrasts[first,key]-contrasts[second,key]).max());normerror=abs(norms[first,key]-norms[second,key]);same=None if row['kind']!='chemical' else (norms[first,key]>.01)==(norms[second,key]>.01)
            passed=error<=limit and normerror<=limit and same is not False
            cr.append(dict(key=key,direct_waveform_error=error,rms_error=normerror,same_scientific_decision=same,passed=bool(passed)))
            previous=next(r for r in saved['factorial_contrasts'] if r['key']==key)
            for name,x in cr[-1].items():
                if name!='key':compare_record(x,previous[name],first+'/'+second+'/'+key+'/'+name)
        passed=all(r['passed'] for r in pulses+cr);compare_record(passed,saved['passed'],first+'/'+second+'/pass')
        comparisons.append(dict(first=first,second=second,pulse_rows=pulses,contrast_rows=cr,passed=bool(passed)))
    accuracy=[];orders=[];rawpairs=[]
    for t in trials:
        key=t['key'];errors=[]
        for level in p['levels']:
            name=level['level'];error=max(log_error(paths[name,key],paths[ref,key]) for ref in ('dop853','radau'));cpu=log_error(paths[name,key],paths['cpu-'+name,key]);errors.append(error)
            row=dict(trial=key,level=name,raw_reference_log_error=error,cpu_torch_log_error=cpu,passed=bool(error<=.01 and cpu<=1e-8));accuracy.append(row)
            previous=next(r for r in summary['production_accuracy'] if (r['trial'],r['level'])==(key,name))
            for k,x in row.items():
                if k not in ('trial','level'):compare_record(x,previous[k],key+'/'+name+'/'+k)
        for first,second in [('coarse','fine'),('fine','diagnostic')]:
            error=log_error(paths[first,key],paths[second,key]);row=dict(trial=key,first=first,second=second,raw_log_error=error,passed=bool(error<=.01));rawpairs.append(row)
            previous=next(r for r in summary['raw_timestep_comparisons'] if (r['trial'],r['first'],r['second'])==(key,first,second))
            for k in ('raw_log_error','passed'):compare_record(row[k],previous[k],key+'/rawpair/'+k)
        resolved=errors[-1]>100*noise[key];order=np.log2(np.array(errors[:-1])/errors[1:]).tolist() if resolved and all(x>0 for x in errors) else None
        passed=not resolved or order is not None and all(1.7<=x<=2.3 for x in order);row=dict(trial=key,errors=errors,reference_noise=noise[key],order_resolved=bool(resolved),orders=order,passed=bool(passed));orders.append(row)
        previous=next(r for r in summary['orders'] if r['trial']==key)
        for k in ('reference_noise','order_resolved','passed'):compare_record(row[k],previous[k],key+'/order/'+k)
        if order is not None:
            for a,b in zip(order,previous['orders']):compare_record(a,b,key+'/order')
        else:assert previous['orders'] is None
    duplicates=[]
    for c in p['contexts']:
        for cell in c['targets']:
            ts=[next(t for t in trials if (t['context'],t['target'],t['dose'])==(c['key'],cell,d)) for d in ('fractional','amount')]
            if not np.array_equal(starts[ts[0]['key']],starts[ts[1]['key']]):continue
            for label in labels:
                error=log_error(paths[label,ts[0]['key']],paths[label,ts[1]['key']])
                duplicates.append(dict(context=c['key'],cell=cell,label=label,raw_log_error=error,passed=bool(error<=1e-8)))
    assert len(duplicates)==120
    archive={'times':times};scientific=[];late=times>=36;baseline=[];immediate=[]
    parent_compact=read_json('docs/phase_carry_recipient_state_response_assessment.json')
    for f in p['factorials']:
        seed,cell,L=f['seed'],f['cell'],f['matching_reference']
        for P,condition in [('S','sham_physics'),('R','redistributed_physics')]:
            cs=[f['contexts'][C+'_'+P] for C in ('S','R')];v=contexts[cs[0]]['volumes'];ids=contexts[cs[0]]['ids'];i=ids.tolist().index(cell)
            a,b=[paths['fine',controls[c]] for c in cs]
            d=np.log(b/a);target=d[:,:,i];network=np.sqrt(np.sum(d**2*v,axis=(1,2))/(2*v.sum()))
            desc=dict(seed=seed,cell=cell,matching_reference=L,physics=P,initial_species_ratios=(b[0,:,i]/a[0,:,i]).tolist(),final_species_ratios=(b[-1,:,i]/a[-1,:,i]).tolist(),late_target_log_rms=rms(target[late],times[late]),late_network_log_rms=float(np.sqrt(trapezoid(network[late]**2,times[late])/24.)))
            baseline.append(desc);saved=next(r for r in summary['baseline_descriptions']['paired_chemical_states'] if (r['seed'],r['cell'],r['matching_reference'],r['physics'])==(seed,cell,L,'sham' if P=='S' else 'redistributed'))
            for akey,bkey in [('late_target_log_rms','late_two_species_state_log_rms'),('late_network_log_rms','late_network_two_species_state_log_rms')]:compare_record(desc[akey],saved[bkey],str((seed,cell,L,P))+'/baseline')
            assert np.array_equal(desc['initial_species_ratios'],saved['initial_species_ratios']) and np.array_equal(desc['final_species_ratios'],saved['final_species_ratios'])
            for dose in ('fractional','amount'):
                key=f'{seed}_{cell}_{L}_{dose}_chemical_{condition}'
                row=dict(seed=seed,cell=cell,matching_reference=L,physics=P,dose=dose,full_effects={label:norms[label,key] for label in labels},above_screen=bool(norms['fine',key]>.01),fine_pulse_metrics={})
                for C,c in zip(('S','R'),cs):
                    t=next(t for t in trials if (t['context'],t['target'],t['dose'])==(c,cell,dose))
                    row['fine_pulse_metrics'][C]=metrics['fine',t['key']];archive[f'{seed}_{cell}_{L}_{P}_{dose}_response_{C}']=waves['fine',t['key']]
                archive[key]=contrasts['fine',key];scientific.append(row)
                if dose=='fractional':
                    old=next(x for x in parent_compact['scientific_comparisons'] if (x['seed'],x['cell'],x['matching_reference'],x['physics'])==(seed,cell,L,P))
                    record=dict(seed=seed,cell=cell,matching_reference=L,physics=condition,immediate_rms=old['full_effects']['fine'],delayed_rms=row['full_effects']['fine'],immediate_above_screen=old['above_screen'],delayed_above_screen=row['above_screen'])
                    saved=next(x for x in summary['immediate_vs_delayed'] if (x['seed'],x['cell'],x['matching_reference'],x['physics'])==(seed,cell,L,condition))
                    for name in ('immediate_rms','delayed_rms','immediate_above_screen','delayed_above_screen'):compare_record(record[name],saved[name],key+'/immediate/'+name)
                    immediate.append(record)
        for dose in ('fractional','amount'):
            for kind,condition in [('physical','sham_chemistry'),('physical','redistributed_chemistry'),('interaction','difference_of_differences')]:
                key=f'{seed}_{cell}_{L}_{dose}_{kind}_{condition}';archive[key]=contrasts['fine',key]
    outcomes=[]
    for seed in p['histories']:
        for L in ('S','R'):
            for dose in ('fractional','amount'):
                for P in ('S','R'):
                    rows=[r for r in scientific if (r['seed'],r['matching_reference'],r['dose'],r['physics'])==(seed,L,dose,P)]
                    count=sum(r['above_screen'] for r in rows);outcomes.append(dict(seed=seed,matching_reference=L,dose=dose,physics=P,above_screen=count,recipients=2))
                    saved=next(r for r in summary['history_outcomes'] if (r['seed'],r['matching_reference'],r['dose'],r['physics'])==(seed,L,dose,'sham_physics' if P=='S' else 'redistributed_physics'))
                    compare_record(count,saved['above_screen'],'history_outcome')
    reference_contrasts=[]
    for seed in p['histories']:
        for cell in p['selections'][str(seed)]['selected_ids']:
            fs={f['matching_reference']:f for f in p['factorials'] if (f['seed'],f['cell'])==(seed,cell)}
            for name in ('S_S','R_S','S_R','R_R'):
                for dose in ('fractional','amount'):
                    ts=[next(t for t in trials if (t['context'],t['target'],t['dose'])==(fs[L]['contexts'][name],cell,dose)) for L in ('S','R')]
                    difference=waves['fine',ts[0]['key']]-waves['fine',ts[1]['key']]
                    row=dict(seed=seed,cell=cell,chemical_physics=name,dose=dose,rms=rms(difference,times),maximum_normalized_difference=float(np.max(abs(difference))))
                    saved=next(x for x in summary['matching_reference_response_contrasts'] if (x['seed'],x['cell'],x['chemical_physics'],x['dose'])==(seed,cell,name,dose))
                    for key in ('rms','maximum_normalized_difference'):compare_record(row[key],saved[key],str((seed,cell,name,dose))+'/reference')
                    reference_contrasts.append(row)
    control_drift=[]
    for context,key in controls.items():
        path=paths['fine',key];v=contexts[context]['volumes']
        d=np.log(path/path[0]);curve=np.sqrt(np.sum(np.mean(d**2,axis=1)*v,axis=1)/v.sum())
        saved=next(x for x in summary['baseline_descriptions']['controls'] if x['context']==context)
        for metric,value in [('final_network_two_species_log_drift_rms',curve[-1]),('maximum_network_two_species_log_drift_rms',curve.max())]:compare_record(value,saved[metric],context+'/'+metric)
        control_drift.append(dict(context=context,maximum_raw_log_drift=float(np.max(abs(d))),maximum_network_log_rms_drift=float(curve.max())))
    for f,h in summary['evidence_sha256'].items():pin(f,h)
    q=read_json(root/'reference-qualification.json');assert q['protocol_sha256']==ph
    for f,h in q['evidence_sha256'].items():pin(f,h)
    overall=all(e<=1e-8 for e in noise.values()) and all(r['passed'] for r in accuracy+rawpairs+orders+comparisons+duplicates)
    compare_record(overall,summary['passed'],'overall');compare_record(overall,status['passed'],'overallstatus')
    for f,h in list(WATCHED.items()):assert digest(f)==h,f
    resolved=[x for r in orders if r['orders'] is not None for x in r['orders']]
    descriptive=[r for r in summary['response_summaries']['fine']['effects'] if r['kind']!='chemical']
    maxima=dict(raw_reference=max(noise.values()),cpu_torch=max(r['cpu_torch_log_error'] for r in accuracy),production_reference=max(r['raw_reference_log_error'] for r in accuracy),reference_normalized=max(r['normalized_response_error'] for r in comparisons[0]['pulse_rows']),normalized_response=max(r['normalized_response_error'] for c in comparisons[1:] for r in c['pulse_rows']),reference_signed_contrast=max(r['direct_waveform_error'] for r in comparisons[0]['contrast_rows']),signed_contrast=max(r['direct_waveform_error'] for c in comparisons[1:] for r in c['contrast_rows']),contrast_rms=max(r['rms_error'] for c in comparisons[1:] for r in c['contrast_rows']),relative_auc=max(r['relative_auc_error'] for c in comparisons for r in c['pulse_rows']),order_min=min(resolved) if resolved else None,order_max=max(resolved) if resolved else None,vector_identity=max(identities),duplicate_log=max(r['raw_log_error'] for r in duplicates))
    gate_counts=dict(reference_raw_pass=sum(e<=1e-8 for e in noise.values()),accuracy_pass=sum(r['passed'] for r in accuracy),order_pass=sum(r['passed'] for r in orders),resolved_orders=sum(r['order_resolved'] for r in orders),raw_pairs_pass=sum(r['passed'] for r in rawpairs),pulse_rows_pass=sum(r['passed'] for c in comparisons for r in c['pulse_rows']),contrast_rows_pass=sum(r['passed'] for c in comparisons for r in c['contrast_rows']))
    result=dict(passed=bool(overall),protocol_sha256=ph,counts=dict(backgrounds=45,logical_starts=141,unique_initial_value_problems=117,production_paths=423,cpu_paths=423,reference_trajectories=282,independent_histories=3,new_histories=0),input_audit=prepared,stage_checks=stages,reference_errors=noise,production_accuracy=accuracy,orders=orders,raw_timestep_comparisons=rawpairs,response_comparisons=comparisons,duplicate_continuity=duplicates,scientific_comparisons=scientific,descriptive_contrasts=descriptive,baseline_descriptions=baseline,history_outcomes=outcomes,immediate_vs_delayed=immediate,matching_reference_response_contrasts=reference_contrasts,control_drift=control_drift,inputs_unchanged=True,inputs_sha256=WATCHED,recorded_metrics_compared=len(RECORDED_DIFFERENCES),maximum_recorded_metric_difference=max(RECORDED_DIFFERENCES)[0],maxima=maxima,gate_counts=gate_counts,script_sha256=digest(__file__),scope='Separate delayed saved-array reconstruction, no project helper imports or integration. Actual fine unpulsed t570 states, no rematching; unequal local challenge chemistry and unmatched surrounding totals. Three existing histories. Unsaved Euler stages remain runtime/source evidence.')
    np.savez_compressed(output.with_suffix('.npz'),**archive);result['waveform_sha256']=digest(output.with_suffix('.npz'))
    output.with_suffix('.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ('passed','gate_counts','maxima','history_outcomes')},indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--root',type=Path,default=Path('output/studies/phase-carry-recipient-delayed-response'));parser.add_argument('--output',type=Path,default=Path('/tmp/recipient_delayed_independent_review'));args=parser.parse_args();main(args.root,args.output)
