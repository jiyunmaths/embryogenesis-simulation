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



read=read_json
arrays=read_npz
def signature(z,u=None):
    parts=[z[k] for k in ('delta','volumes','ids')]+[z['initial'] if u is None else u]
    return tuple((a.dtype.str,a.shape,a.tobytes()) for a in parts)

def audit(root):
    import importlib.util
    file=Path('archive/study-launches/2026-10-08-recipient-amount-response/audit_inputs.py')
    spec=importlib.util.spec_from_file_location('literal_amount_audit',file)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    result=module.audit(root);WATCHED.update(result['inputs_sha256']);pin(file)
    return result


def amount_curves(p,contexts,paths,controls,times,summary):
    rows=[]
    for f in p['factorials']:
        for P in ('S','R'):
            cs=[f['contexts'][C+'_'+P] for C in ('S','R')];v=contexts[cs[0]]['volumes'];ids=contexts[cs[0]]['ids'];others=ids!=f['cell']
            a,b=[paths['fine',controls[c]] for c in cs]
            ga,gb=a@v,b@v;sa,sb=a[:,:,others]@v[others],b[:,:,others]@v[others]
            gr=gb/ga;sr=sb/sa;error=max(float(np.max(abs(gr[0]-1))),float(np.max(abs(sr[0]-1))))
            if p['phase']=='immediate':assert error<=2e-13
            row=dict(seed=f['seed'],cell=f['cell'],matching_reference=f['matching_reference'],physics=P,times=times.tolist(),global_species_amounts_S=ga.tolist(),global_species_amounts_R=gb.tolist(),surrounding_species_amounts_S=sa.tolist(),surrounding_species_amounts_R=sb.tolist(),global_species_amount_ratios=gr.tolist(),surrounding_species_amount_ratios=sr.tolist(),initial_global_species_amount_ratios=gr[0].tolist(),final_global_species_amount_ratios=gr[-1].tolist(),maximum_global_species_log_amount_separation=np.max(abs(np.log(gr)),axis=0).tolist(),initial_paired_amount_relative_error=error)
            saved=next(x for x in summary['amount_evolution'] if all(x[k]==row[k] for k in ('seed','cell','matching_reference','physics')))
            for name,value in row.items():
                if name in ('seed','cell','matching_reference','physics'):continue
                np.testing.assert_allclose(value,saved[name],rtol=1e-12,atol=2e-13)
            rows.append(row)
    return rows



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
    assert len(groups)==9 and max(map(len,groups.values()))==16
    for level in p['levels']:
        label=level['level'];folder=root/'production'/label;r=read_json(folder/'result.json');z=read_npz(folder/'paths.npz');cp=read_npz(folder/'checkpoint.npz');commit=read_json(folder/'checkpoint-commit.json')
        assert r['protocol_sha256']==ph and r['level']==level and r['paths_sha256']==digest(folder/'paths.npz') and r['checkpoint_sha256']==digest(folder/'checkpoint.npz')
        assert commit['protocol_sha256']==ph and commit['checkpoint_sha256']==r['checkpoint_sha256']
        assert list(z['trial_keys'])==[t['key'] for t in trials] and np.array_equal(z['times'],times)
        valid(z['trajectory'],(401,96,2,16));valid(cp['state'],(2,144,16));assert np.array_equal(z['trajectory'],cp['history'])
        meta=json.loads(str(cp['metadata']));assert meta['step']==round(60/level['dt']) and commit['step']==meta['step'] and meta['protocol_sha256']==ph and meta['level']==level
        for i,t in enumerate(trials):
            key=t['key'];gi,col=mapping[key];path=z['trajectory'][:,i]
            assert np.array_equal(path[0],starts[key]) and np.array_equal(path[-1],cp['state'][:,16*gi:16*(gi+1),col]);paths[label,key]=path
        padding_error=max([0.]+[log_error(cp['state'][:,16*gi:16*(gi+1),0],cp['state'][:,16*gi:16*(gi+1),col]) for gi,group in enumerate(groups.values()) for col in range(len(group),16)])
        assert padding_error<=1e-8
        assert r['actual_paths']==96 and r['padding_paths']==48
        stages.append(dict(level=label,logical_paths=96,padding=48,padding_error=padding_error,positive_saved_states=True,checkpoint_consistent=True,wall_seconds=r['wall_seconds']))
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
            seed,cell,L=f['seed'],f['cell'],f['matching_reference'];four={}
            for name,context in f['contexts'].items():
                t=next(t for t in trials if t['context']==context and t['target']==cell);four[name]=waves[label,t['key']]
            cc={P:four['R_'+P]-four['S_'+P] for P in ('S','R')};pp={C:four[C+'_R']-four[C+'_S'] for C in ('S','R')};interaction=four['R_R']-four['S_R']-four['R_S']+four['S_S'];diagonal=four['R_R']-four['S_S']
            kinds=[('chemical','sham_physics',cc['S']),('chemical','redistributed_physics',cc['R']),('physical','sham_chemistry',pp['S']),('physical','redistributed_chemistry',pp['R']),('interaction','difference_of_differences',interaction)]
            for kind,condition,wave in kinds:
                key=f'{seed}_{cell}_{L}_fractional_{kind}_{condition}';contrasts[label,key]=wave;norms[label,key]=rms(wave,times)
            diagonals[label,seed,cell,L]=rms(diagonal,times)
            identities.extend(float(abs(x).max()) for x in [interaction-(cc['R']-cc['S']),interaction-(pp['R']-pp['S']),diagonal-cc['S']-pp['R'],diagonal-cc['R']-pp['S']])
        recorded=summary['response_summaries'][label]
        for t in trials:
            if t['target'] is None:continue
            for name,value in metrics[label,t['key']].items():compare_record(value,recorded['metrics'][t['key']][name],label+'/'+t['key']+'/'+name)
        assert len(recorded['effects'])==60
        for row in recorded['effects']:
            value=norms[label,row['key']];compare_record(value,row['normalized_response_shift_rms'],label+'/'+row['key'])
            if row['kind']=='chemical':compare_record(value>.01,row['response_shift_detected'],label+'/'+row['key']+'/decision')
            else:assert 'response_shift_detected' not in row
        for row in recorded['native_diagonal']:compare_record(diagonals[label,row['seed'],row['cell'],row['matching_reference']],row['normalized_response_shift_rms'],label+'/diagonal')
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
    parent=Path(p['parent']);old=read_json(parent/'protocol.json');native=[]
    for t in trials:
        if 'native_parent_trial' not in t:continue
        previous=t['native_parent_trial']
        for label in labels:
            if label in ('dop853','radau'):
                folder=parent/'references'/previous;rr=read_json(folder/'result.json');z=read_npz(folder/'paths.npz');other=z[label]
            else:
                folder=parent/'production'/label;rr=read_json(folder/'result.json');z=read_npz(folder/'paths.npz');other=z['trajectory'][:,list(z['trial_keys']).index(previous)]
            assert rr['paths_sha256']==digest(folder/'paths.npz') and np.array_equal(z['times'],times)
            error=log_error(paths[label,t['key']],other);native.append(dict(trial=t['key'],label=label,error=error,passed=bool(error<=1e-8)))
            if label in ('coarse','fine'):
                saved=next(r for r in summary['native_parent_comparisons'] if (r['trial'],r['level'])==(t['key'],label));compare_record(error,saved['raw_log_error'],t['key']+'/native/'+label)
    assert len(native)==0
    archive={'times':times};scientific=[];late=times>=36;baseline=[]
    parent_compact=read_json('docs/phase_carry_recipient_state_response_assessment.json')
    for f in p['factorials']:
        seed,cell,L=f['seed'],f['cell'],f['matching_reference']
        for P,condition in [('S','sham_physics'),('R','redistributed_physics')]:
            key=f'{seed}_{cell}_{L}_fractional_chemical_{condition}'
            row=dict(seed=seed,cell=cell,matching_reference=L,physics=P,full_effects={label:norms[label,key] for label in labels},above_screen=bool(norms['fine',key]>.01))
            unmatched=next(r for r in parent_compact['scientific_comparisons'] if (r['seed'],r['cell'],r['matching_reference'],r['physics'])==(seed,cell,L,P))
            row['parent_unmatched_fractional_rms']=unmatched['full_effects']['fine'];row['fine_pulse_metrics']={}
            cs=[f['contexts'][C+'_'+P] for C in ('S','R')];v=contexts[cs[0]]['volumes'];ids=contexts[cs[0]]['ids'];i=ids.tolist().index(cell)
            a,b=[paths['fine',controls[c]] for c in cs];assert np.array_equal(a[0,:,i],b[0,:,i])
            d=np.log(b/a);target=d[:,:,i];network=np.sqrt(np.sum(d**2*v,axis=(1,2))/(2*v.sum()))
            desc=dict(seed=seed,cell=cell,matching_reference=L,physics=P,initial_species_ratios=(b[0,:,i]/a[0,:,i]).tolist(),final_species_ratios=(b[-1,:,i]/a[-1,:,i]).tolist(),late_target_log_rms=rms(target[late],times[late]),late_network_log_rms=float(np.sqrt(trapezoid(network[late]**2,times[late])/24.)))
            baseline.append(desc);saved=next(r for r in summary['baseline_descriptions']['paired_chemical_states'] if (r['seed'],r['cell'],r['matching_reference'],r['physics'])==(seed,cell,L,'sham' if P=='S' else 'redistributed'))
            for akey,bkey in [('late_target_log_rms','late_two_species_state_log_rms'),('late_network_log_rms','late_network_two_species_state_log_rms')]:compare_record(desc[akey],saved[bkey],key+'/baseline')
            assert np.array_equal(desc['initial_species_ratios'],saved['initial_species_ratios']) and np.array_equal(desc['final_species_ratios'],saved['final_species_ratios'])
            for C,c in zip(('S','R'),cs):
                t=next(t for t in trials if t['context']==c and t['target']==cell);row['fine_pulse_metrics'][C]=metrics['fine',t['key']];archive[f'{seed}_{cell}_{L}_{P}_response_{C}']=waves['fine',t['key']]
            archive[key]=contrasts['fine',key];scientific.append(row)
        for kind,condition in [('physical','sham_chemistry'),('physical','redistributed_chemistry'),('interaction','difference_of_differences')]:
            key=f'{seed}_{cell}_{L}_fractional_{kind}_{condition}';archive[key]=contrasts['fine',key]
    outcomes=[]
    for seed in p['histories']:
        for L in ('S','R'):
            for P in ('S','R'):
                rows=[r for r in scientific if (r['seed'],r['matching_reference'],r['physics'])==(seed,L,P)]
                count=sum(r['above_screen'] for r in rows);outcomes.append(dict(seed=seed,matching_reference=L,physics=P,above_screen=count,recipients=2))
                saved=next(r for r in summary['history_outcomes'] if (r['seed'],r['matching_reference'],r['physics'])==(seed,L,'sham_physics' if P=='S' else 'redistributed_physics'))
                compare_record(count,saved['above_screen'],'history_outcome')
    for row in summary['reciprocal_history_outcomes']:
        count=sum(all(r['above_screen'] for r in scientific if (r['seed'],r['cell'])==(row['seed'],cell)) for cell in p['selections'][str(row['seed'])]['selected_ids'])
        compare_record(count,row['all_four_above_screen'],'allfour_outcome')
    for f,h in summary['evidence_sha256'].items():pin(f,h)
    q=read_json(root/'reference-qualification.json');assert q['protocol_sha256']==ph
    for f,h in q['evidence_sha256'].items():pin(f,h)
    overall=all(e<=1e-8 for e in noise.values()) and all(r['passed'] for r in accuracy+rawpairs+orders+comparisons+native)
    compare_record(overall,summary['passed'],'overall');compare_record(overall,status['passed'],'overallstatus')
    for f,h in list(WATCHED.items()):assert digest(f)==h,f
    resolved=[x for r in orders if r['orders'] is not None for x in r['orders']]
    descriptive=[r for r in summary['response_summaries']['fine']['effects'] if r['kind']!='chemical']
    result=dict(passed=bool(overall),protocol_sha256=ph,counts=dict(backgrounds=48,logical_starts=96,unique_initial_value_problems=p['unique_initial_value_problems'],unique_physical_initial_value_problems=p['unique_physical_initial_value_problems'],production_paths=288,cpu_paths=288,reference_trajectories=192,independent_histories=3,new_histories=0),input_audit=prepared,stage_checks=stages,reference_errors=noise,production_accuracy=accuracy,orders=orders,raw_timestep_comparisons=rawpairs,response_comparisons=comparisons,native_continuity=native,scientific_comparisons=scientific,descriptive_contrasts=descriptive,baseline_descriptions=baseline,history_outcomes=outcomes,amount_evolution=amount_curves(p,contexts,paths,controls,times,summary),inputs_unchanged=True,inputs_sha256=WATCHED,recorded_metrics_compared=len(RECORDED_DIFFERENCES),maximum_recorded_metric_difference=max(RECORDED_DIFFERENCES)[0],maxima=dict(raw_reference=max(noise.values()),cpu_torch=max(r['cpu_torch_log_error'] for r in accuracy),production_reference=max(r['raw_reference_log_error'] for r in accuracy),reference_normalized=max(r['normalized_response_error'] for r in comparisons[0]['pulse_rows']),normalized_response=max(r['normalized_response_error'] for c in comparisons[1:] for r in c['pulse_rows']),reference_signed_contrast=max(r['direct_waveform_error'] for r in comparisons[0]['contrast_rows']),signed_contrast=max(r['direct_waveform_error'] for c in comparisons[1:] for r in c['contrast_rows']),contrast_rms=max(r['rms_error'] for c in comparisons[1:] for r in c['contrast_rows']),relative_auc=max(r['relative_auc_error'] for c in comparisons for r in c['pulse_rows']),order_min=min(resolved),order_max=max(resolved),vector_identity=max(identities),native_log=None),gate_counts=dict(reference_raw_pass=sum(e<=1e-8 for e in noise.values()),accuracy_pass=sum(r['passed'] for r in accuracy),order_pass=sum(r['passed'] for r in orders),resolved_orders=sum(r['order_resolved'] for r in orders),raw_pairs_pass=sum(r['passed'] for r in rawpairs),pulse_rows_pass=sum(r['passed'] for c in comparisons for r in c['pulse_rows']),contrast_rows_pass=sum(r['passed'] for c in comparisons for r in c['contrast_rows'])),script_sha256=digest(__file__),scope='Separate savedarray reconstruction, no project helper imports or integration. Identical initial recipient pairs, equal initial surrounding and global species budgets, distinct preserved species allocation fractions, immediate pulse/free evolution; later totals endogenous. Unsaved Euler stages remain runtime/source evidence, not independently reconstructed observations.')
    np.savez_compressed(output.with_suffix('.npz'),**archive);result['waveform_sha256']=digest(output.with_suffix('.npz'))
    output.with_suffix('.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ('passed','gate_counts','maxima','history_outcomes')},indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--root',type=Path,default=Path('output/studies/phase-carry-recipient-amount-response'));parser.add_argument('--output',type=Path,default=Path('/tmp/recipient_amount_independent_review'));args=parser.parse_args();main(args.root,args.output)
