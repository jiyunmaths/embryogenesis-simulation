"""Challenge-time amount-matched saved-array review; no project helpers or new integration.

Literal log-response formulas adapted from the archived delayed-response review.
Preparation is checked by the separately archived literal endpoint auditor.
"""
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


def correction_ok(state, correction):
    assert state.dtype == correction.dtype == np.float64
    assert correction.shape == state.shape and np.isfinite(correction).all()
    assert np.all(abs(correction) <= 2*np.spacing(state))


def preparation():
    import importlib.util
    file=Path('archive/study-launches/2026-10-09-pulse-time-amount-response/audit.py')
    spec=importlib.util.spec_from_file_location('literal_endpoint_audit',file)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    p,detail=module.inputs()
    WATCHED.update(module.PINS);pin(file)
    root=module.ROOT;prefix=read_json(root/'prefix/independent-audit.json')
    assert prefix['passed'] and prefix['actual_context_timestep_cases']==1008
    for f,h in prefix['evidence_sha256'].items():pin(f,h)
    assert all(r['passed'] and r['restart_bitwise_identical'] for r in prefix['rows'])
    detail.update(source_pins=len(p['source_sha256']),input_pins=len(p['input_sha256']),
                  actual_start_prefix_cases=1008,prefix_audit_sha256=digest(root/'prefix/independent-audit.json'))
    return detail


def main(root, output):
    root=Path(root).resolve();output=Path(output)
    prepared=preparation()
    p=read_json(root/'protocol.json');ph=digest(root/'protocol.json')
    status=read_json(root/'status.json');summary=read_json(root/'summary.json')
    assert status['state']=='completed' and summary['protocol_sha256']==ph
    limits=dict(raw_log_max=.01,normalized_response_max=.01,relative_auc_max=.02,
        recovery_time_max=.3,effect_error_max=.01,reference_log_max=1e-8,
        reference_normalized_max=1e-8,cpu_torch_log_max=1e-8,order_min=1.7,
        order_max=2.3,reference_noise_factor=100.,pulse_amount_relative_max=2e-13)
    assert p['criteria']==limits and p['effect_min']==.01
    assert (p['beta'],p['da'],p['db'])==(2.,.02,.55)
    times=np.arange(401)*.15;trials=p['trials'];labels=('dop853','radau','coarse','fine','diagnostic')
    method='float64 SSP-RK2 increment with Kahan accumulation v1'
    assert summary['method']==method
    contexts={c['key']:read_npz(c['source']) for c in p['contexts']}
    for z in contexts.values():
        assert z['ids'].shape==(16,) and len(set(z['ids'].tolist()))==16
        assert z['volumes'].shape==(16,) and np.all(z['volumes']>0)
        assert z['delta'].shape==(16,16) and np.isfinite(z['delta']).all()
        off=z['delta'].copy();np.fill_diagonal(off,0.);assert np.all(off>=0)
    starts={t['key']:read_npz(t['source'])['initial'] for t in trials}
    controls={t['context']:t['key'] for t in trials if t['target'] is None}
    pulse_keys={(t['context'],t['target']):t['key'] for t in trials if t['target'] is not None}
    paths={};noise={};stages=[]
    for t in trials:
        folder=root/'references'/t['key'];r=read_json(folder/'result.json');z=read_npz(folder/'paths.npz')
        assert r['protocol_sha256']==ph and r['trial']==t and r['paths_sha256']==digest(folder/'paths.npz')
        np.testing.assert_array_equal(z['times'],times)
        for label in ('dop853','radau'):
            valid(z[label],(401,2,16));np.testing.assert_array_equal(z[label][0],starts[t['key']]);paths[label,t['key']]=z[label]
        noise[t['key']]=log_error(z['dop853'],z['radau'])
        compare_record(noise[t['key']],r['raw_solver_log_error'],t['key']+'/reference')
    groups={};mapping={};distinct={}
    for t in trials:
        z=contexts[t['context']];sig=tuple((z[k].dtype.str,z[k].shape,z[k].tobytes()) for k in ('delta','volumes','ids'))
        if sig not in groups:groups[sig]=[]
        mapping[t['key']]=(list(groups).index(sig),len(groups[sig]));groups[sig].append(t['key'])
        distinct.setdefault((sig,starts[t['key']].tobytes()),[]).append(t['key'])
    assert len(groups)==12 and {len(g) for g in groups.values()}=={28}
    assert len(distinct)==p['unique_physical_initial_value_problems']==192
    for level in p['levels']:
        label=level['level'];folder=root/'production'/label
        r=read_json(folder/'result.json');z=read_npz(folder/'paths.npz');cp=read_npz(folder/'checkpoint.npz');commit=read_json(folder/'checkpoint-commit.json')
        assert r['protocol_sha256']==ph and r['level']==level and r['method']==method
        assert r['paths_sha256']==digest(folder/'paths.npz') and r['checkpoint_sha256']==digest(folder/'checkpoint.npz')
        assert commit['protocol_sha256']==ph and commit['checkpoint_sha256']==r['checkpoint_sha256']
        np.testing.assert_array_equal(z['trial_keys'],[t['key'] for t in trials]);np.testing.assert_array_equal(z['times'],times)
        valid(z['trajectory'],(401,336,2,16));valid(cp['state'],(2,192,28));correction_ok(cp['state'],cp['correction'])
        np.testing.assert_array_equal(z['trajectory'],cp['history'])
        meta=json.loads(str(cp['metadata']));assert meta['step']==round(60/level['dt']) and commit['step']==meta['step']
        assert meta['protocol_sha256']==ph and meta['level']==level and meta['method']==method
        assert r['actual_paths']==336 and r['padding_paths']==0 and r['correction_checkpointed']
        for i,t in enumerate(trials):
            key=t['key'];gi,col=mapping[key];path=z['trajectory'][:,i]
            np.testing.assert_array_equal(path[0],starts[key]);np.testing.assert_array_equal(path[-1],cp['state'][:,16*gi:16*(gi+1),col]);paths[label,key]=path
        stages.append(dict(level=label,logical_paths=336,padding=0,positive_saved_states=True,
            checkpoint_consistent=True,correction_valid=True,wall_seconds=r['wall_seconds']))
        for t in trials:
            folder=root/'cpu-verification'/label/t['key'];r=read_json(folder/'result.json');z=read_npz(folder/'paths.npz')
            assert r['protocol_sha256']==ph and r['trial']==t and r['level']==level and r['method']==method
            assert r['paths_sha256']==digest(folder/'paths.npz');np.testing.assert_array_equal(z['times'],times)
            valid(z['trajectory'],(401,2,16));np.testing.assert_array_equal(z['trajectory'][0],starts[t['key']])
            correction_ok(z['trajectory'][-1],z['final_correction']);paths['cpu-'+label,t['key']]=z['trajectory']
    waves={};changes={};metrics={};contrasts={};norms={};diagonals={};identities=[]
    for label in labels:
        for t in trials:
            if t['target'] is None:continue
            z=contexts[t['context']];i=z['ids'].tolist().index(t['target']);key=t['key']
            waves[label,key],changes[label,key],metrics[label,key]=response(paths[label,key],paths[label,controls[t['context']]],times,z['volumes'],i,t['factor'],24.)
        for f in p['factorials']:
            seed,cell,M,B=f['seed'],f['cell'],f['matching_reference'],f['budget_reference']
            four={n:waves[label,pulse_keys[c,cell]] for n,c in f['contexts'].items()}
            cc={P:four['R_'+P]-four['S_'+P] for P in ('S','R')};pp={C:four[C+'_R']-four[C+'_S'] for C in ('S','R')}
            interaction=four['R_R']-four['S_R']-four['R_S']+four['S_S'];diagonal=four['R_R']-four['S_S']
            kinds=[('chemical','sham_physics',cc['S']),('chemical','redistributed_physics',cc['R']),
                ('physical','sham_chemistry',pp['S']),('physical','redistributed_chemistry',pp['R']),
                ('interaction','difference_of_differences',interaction)]
            for kind,condition,wave in kinds:
                key=f'{seed}_{cell}_{M}_budget-{B}_fractional_{kind}_{condition}';contrasts[label,key]=wave;norms[label,key]=rms(wave,times)
            diagonals[label,seed,cell,M,B]=rms(diagonal,times)
            identities.extend(float(abs(x).max()) for x in [interaction-(cc['R']-cc['S']),interaction-(pp['R']-pp['S']),diagonal-cc['S']-pp['R'],diagonal-cc['R']-pp['S']])
        for f in p['factorials']:
            if f['matching_reference']=='native' or f['budget_reference']=='native':continue
            for condition in ('sham_physics','redistributed_physics'):
                prefix=f"{f['seed']}_{f['cell']}_{f['matching_reference']}_budget-"
                original=prefix+'native_fractional_chemical_'+condition
                matched=prefix+f['budget_reference']+'_fractional_chemical_'+condition
                key=prefix+f['budget_reference']+'_fractional_budget_change_'+condition
                contrasts[label,key]=contrasts[label,matched]-contrasts[label,original]
                norms[label,key]=rms(contrasts[label,key],times)
        recorded=summary['response_summaries'][label]
        for t in trials:
            if t['target'] is None:continue
            for name,value in metrics[label,t['key']].items():compare_record(value,recorded['metrics'][t['key']][name],label+'/'+t['key']+'/'+name)
        assert len(recorded['effects'])==258
        for row in recorded['effects']:
            value=norms[label,row['key']];compare_record(value,row['normalized_response_shift_rms'],label+'/'+row['key'])
            if row['kind']=='chemical':compare_record(value>.01,row['response_shift_detected'],label+'/'+row['key']+'/decision')
            else:assert 'response_shift_detected' not in row
        for row in recorded['native_diagonal']:compare_record(diagonals[label,row['seed'],row['cell'],row['matching_reference'],row['budget_reference']],row['normalized_response_shift_rms'],label+'/diagonal')
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
            compare_record(limit,previous['direct_waveform_limit'],'waveform_limit');compare_record(limit,previous['rms_limit'],'rms_limit')
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
        np.testing.assert_array_equal(errors,previous['errors'])
        if order is not None:
            for a,b in zip(order,previous['orders']):compare_record(a,b,key+'/order')
        else:assert previous['orders'] is None
    duplicates=[]
    for group in distinct.values():
        for key in group[1:]:
            for label in labels:
                error=log_error(paths[label,group[0]],paths[label,key]);assert error<=1e-8
                np.testing.assert_array_equal(paths[label,group[0]],paths[label,key])
                duplicates.append(dict(first=group[0],second=key,label=label,raw_log_error=error,passed=True))
    assert len(duplicates)==720
    archive={'times':times};scientific=[];late=times>=36;baseline=[];amounts=[]
    for f in p['factorials']:
        seed,cell,M,B=f['seed'],f['cell'],f['matching_reference'],f['budget_reference']
        for P,condition in [('S','sham_physics'),('R','redistributed_physics')]:
            cs=[f['contexts'][C+'_'+P] for C in ('S','R')];v=contexts[cs[0]]['volumes'];ids=contexts[cs[0]]['ids'];i=ids.tolist().index(cell)
            a,b=[paths['fine',controls[c]] for c in cs]
            d=np.log(b/a);target=d[:,:,i];network=np.sqrt(np.sum(d**2*v,axis=(1,2))/(2*v.sum()))
            desc=dict(seed=seed,cell=cell,matching_reference=M,budget_reference=B,physics=P,initial_species_ratios=(b[0,:,i]/a[0,:,i]).tolist(),final_species_ratios=(b[-1,:,i]/a[-1,:,i]).tolist(),late_target_log_rms=rms(target[late],times[late]),late_network_log_rms=float(np.sqrt(trapezoid(network[late]**2,times[late])/24.)))
            baseline.append(desc);saved=next(r for r in summary['baseline_descriptions']['paired_chemical_states'] if (r['seed'],r['cell'],r['matching_reference'],r['budget_reference'],r['physics'])==(seed,cell,M,B,'sham' if P=='S' else 'redistributed'))
            for akey,bkey in [('late_target_log_rms','late_two_species_state_log_rms'),('late_network_log_rms','late_network_two_species_state_log_rms')]:compare_record(desc[akey],saved[bkey],str((seed,cell,M,P))+'/baseline')
            np.testing.assert_array_equal(desc['initial_species_ratios'],saved['initial_species_ratios']);np.testing.assert_array_equal(desc['final_species_ratios'],saved['final_species_ratios'])
            if M!='native':np.testing.assert_array_equal(desc['initial_species_ratios'],[1.,1.])
            key=f'{seed}_{cell}_{M}_budget-{B}_fractional_chemical_{condition}'
            row=dict(seed=seed,cell=cell,matching_reference=M,budget_reference=B,physics=P,dose='fractional',full_effects={label:norms[label,key] for label in labels},above_screen=bool(norms['fine',key]>.01),fine_pulse_metrics={})
            for C,c in zip(('S','R'),cs):
                tkey=pulse_keys[c,cell];row['fine_pulse_metrics'][C]=metrics['fine',tkey];archive[f'{seed}_{cell}_{M}_budget-{B}_{P}_response_{C}']=waves['fine',tkey]
                archive[f'{seed}_{cell}_{M}_budget-{B}_{P}_control_{C}']=paths['fine',controls[c]]
            archive[key]=contrasts['fine',key];scientific.append(row)
            others=ids!=cell;ga,gb=a@v,b@v;sa,sb=a[:,:,others]@v[others],b[:,:,others]@v[others]
            prefix=f'{seed}_{cell}_{M}_budget-{B}_{P}';archive[prefix+'_global_amounts_S']=ga;archive[prefix+'_global_amounts_R']=gb
            archive[prefix+'_surrounding_amounts_S']=sa;archive[prefix+'_surrounding_amounts_R']=sb
            amounts.append(dict(seed=seed,cell=cell,matching_reference=M,budget_reference=B,physics=P,
                initial_global_species_amount_ratios=(gb[0]/ga[0]).tolist(),final_global_species_amount_ratios=(gb[-1]/ga[-1]).tolist(),
                initial_surrounding_species_amount_ratios=(sb[0]/sa[0]).tolist(),final_surrounding_species_amount_ratios=(sb[-1]/sa[-1]).tolist(),
                maximum_global_species_log_amount_separation=np.max(abs(np.log(gb/ga)),axis=0).tolist()))
            saved=next(x for x in summary['amount_evolution'] if (x['seed'],x['cell'],x['matching_reference'],x['budget_reference'],x['physics'])==(seed,cell,M,B,P))
            for k,x in dict(times=times,global_species_amounts_S=ga,global_species_amounts_R=gb,
                surrounding_species_amounts_S=sa,surrounding_species_amounts_R=sb).items():np.testing.assert_allclose(x,saved[k],rtol=1e-12,atol=2e-13)
            for k,x in amounts[-1].items():
                if k in saved and k not in ('seed','cell','matching_reference','budget_reference','physics'):
                    np.testing.assert_allclose(x,saved[k],rtol=1e-12,atol=2e-13)
            if B!='native':
                assert M!='native'
                assert max(float(abs(gb[0]/ga[0]-1).max()),float(abs(sb[0]/sa[0]-1).max()))<=2e-13
    outcomes=[]
    for seed in p['histories']:
        for M,B in [('native','native'),*[(M,B) for M in ('S','R') for B in ('native','S','R')]]:
            for P in ('S','R'):
                rows=[r for r in scientific if (r['seed'],r['matching_reference'],r['budget_reference'],r['physics'])==(seed,M,B,P)]
                assert len(rows)==2
                count=sum(r['above_screen'] for r in rows)
                outcomes.append(dict(history=seed,matching_reference=M,budget_reference=B,physics=P,above_screen=count,recipients=2))
                saved=next(r for r in summary['history_outcomes'] if (r['history'],r['matching_reference'],r['budget_reference'],r['physical_context'])==(seed,M,B,'sham_physics' if P=='S' else 'redistributed_physics'))
                compare_record(count,saved['above_screen'],'history_outcome')
                assert [(r['cell'],r['above_screen']) for r in rows]==[(r['cell'],r['above_screen']) for r in saved['cells']]
    control_drift=[]
    for context,key in controls.items():
        path=paths['fine',key];v=contexts[context]['volumes'];z=contexts[context]
        d=np.log(path/path[0]);curve=np.sqrt(np.sum(np.mean(d**2,axis=1)*v,axis=1)/v.sum());sd=np.std(np.log(path[:,0]),axis=1)
        saved=next(x for x in summary['baseline_descriptions']['controls'] if x['context']==context)
        np.testing.assert_array_equal(saved['final_chemistry'],path[-1]);np.testing.assert_array_equal(saved['ids'],z['ids'])
        np.testing.assert_array_equal(saved['final_global_species_amounts'],path[-1]@v)
        np.testing.assert_array_equal(saved['final_per_cell_species_ratios_to_reassigned_start'],path[-1]/path[0])
        values={'final_network_two_species_log_drift_rms':curve[-1],'maximum_network_two_species_log_drift_rms':curve.max(),
            'late_network_two_species_log_drift_rms':np.sqrt(trapezoid(curve[late]**2,times[late])/24.),
            'final_log_activator_sd':sd[-1],'whole_min_log_activator_sd':sd.min(),'late_min_log_activator_sd':sd[late].min()}
        for k,x in values.items():compare_record(x,saved[k],context+'/'+k)
        control_drift.append(dict(context=context,maximum_raw_log_drift=float(np.max(abs(d))),**{k:float(x) for k,x in values.items()}))
    reference_contrasts=[];endpoint_reference_differences=[]
    for f in [f for f in p['factorials'] if f['matching_reference']=='native']:
        seed,cell=f['seed'],f['cell']
        fs={(x['matching_reference'],x['budget_reference']):x for x in p['factorials'] if (x['seed'],x['cell'])==(seed,cell)}
        for name in ('S_S','R_S','S_R','R_R'):
            for B in ('native','S','R'):
                a,b=[waves['fine',pulse_keys[fs[M,B]['contexts'][name],cell]] for M in ('S','R')]
                reference_contrasts.append(dict(seed=seed,cell=cell,budget_reference=B,chemical_physics=name,rms=rms(a-b,times),maximum_normalized_difference=float(np.max(abs(a-b)))))
            ends=[paths['fine',controls[x['contexts'][name]]][-1] for x in fs.values()]
            endpoint_reference_differences.append(dict(seed=seed,cell=cell,chemical_physics=name,maximum_final_log_separation=max(log_error(a,b) for a in ends for b in ends)))
    descriptive=[r.copy() for r in summary['response_summaries']['fine']['effects'] if r['kind']!='chemical']
    for r in descriptive:archive[r['key']]=contrasts['fine',r['key']]
    # All unmatched-budget comparators are newly computed. Check their exact
    # accepted-parent starts and complete curves, rather than assuming reuse.
    parent=Path(p['parent']);old=read_json(parent/'protocol.json')
    old_review=read_json('docs/pulse_time_recipient_response_assessment.json')
    pin(old_review['independent_review'],old_review['independent_review_sha256']);old_audit=read_json(old_review['independent_review'])
    parent_paths={};parent_trials={};parent_contexts={c['key']:c for c in old['contexts']}
    for t in old['trials']:parent_trials[t['context'],t['target']]=t
    for label in labels:
        if label in ('dop853','radau'):
            for t in old['trials']:
                file=parent/'references'/t['key']/'paths.npz';pin(file,old_audit['evidence_sha256'][str(file.resolve())])
                z=read_npz(file,('times',label));np.testing.assert_array_equal(z['times'],times);parent_paths[label,t['key']]=z[label]
        else:
            file=parent/'production'/label/'paths.npz';pin(file,old_audit['evidence_sha256'][str(file.resolve())])
            z=read_npz(file);np.testing.assert_array_equal(z['times'],times)
            for j,key in enumerate(z['trial_keys']):parent_paths[label,str(key)]=z['trajectory'][:,j]
    parent_comparisons=[]
    context_records={c['key']:c for c in p['contexts']}
    for t in trials:
        if t['budget_reference']!='native':continue
        key=context_records[t['context']]['parent_context'];previous=parent_trials[key,t['target']]
        np.testing.assert_array_equal(starts[t['key']],read_npz(previous['source'])['initial'])
        for label in labels:
            error=log_error(paths[label,t['key']],parent_paths[label,previous['key']]);assert error<=1e-8
            parent_comparisons.append(dict(trial=t['key'],parent_trial=previous['key'],label=label,raw_log_error=error,passed=True))
    assert len(parent_comparisons)==720
    q=read_json(root/'reference-qualification.json');assert q['protocol_sha256']==ph and q['passed']
    for record in (summary,q):
        for f,h in record['evidence_sha256'].items():pin(f,h)
    overall=all(e<=1e-8 for e in noise.values()) and all(r['passed'] for r in accuracy+rawpairs+orders+comparisons+duplicates)
    compare_record(overall,summary['passed'],'overall');compare_record(overall,status['passed'],'overallstatus')
    resolved=[x for r in orders if r['orders'] is not None for x in r['orders']]
    maxima=dict(raw_reference=max(noise.values()),cpu_torch=max(r['cpu_torch_log_error'] for r in accuracy),
        production_reference=max(r['raw_reference_log_error'] for r in accuracy),
        reference_normalized=max(r['normalized_response_error'] for r in comparisons[0]['pulse_rows']),
        normalized_response=max(r['normalized_response_error'] for c in comparisons[1:] for r in c['pulse_rows']),
        reference_signed_contrast=max(r['direct_waveform_error'] for r in comparisons[0]['contrast_rows']),
        signed_contrast=max(r['direct_waveform_error'] for c in comparisons[1:] for r in c['contrast_rows']),
        contrast_rms=max(r['rms_error'] for c in comparisons[1:] for r in c['contrast_rows']),
        relative_auc=max(r['relative_auc_error'] for c in comparisons for r in c['pulse_rows']),
        recovery_time=max(x for c in comparisons for r in c['pulse_rows'] for x in r['recovery_time_errors'].values()),
        order_min=min(resolved),order_max=max(resolved),vector_identity=max(identities),duplicate_log=max(r['raw_log_error'] for r in duplicates))
    gate_counts=dict(reference_raw_pass=sum(e<=1e-8 for e in noise.values()),accuracy_pass=sum(r['passed'] for r in accuracy),
        order_pass=sum(r['passed'] for r in orders),resolved_orders=sum(r['order_resolved'] for r in orders),
        raw_pairs_pass=sum(r['passed'] for r in rawpairs),pulse_rows_pass=sum(r['passed'] for c in comparisons for r in c['pulse_rows']),
        contrast_rows_pass=sum(r['passed'] for c in comparisons for r in c['contrast_rows']),duplicate_checks_pass=len(duplicates))
    assert gate_counts['resolved_orders']==336
    timing=read_json(root/'timing.json')
    preserved_scope={str(f):digest(f) for f in [Path('docs/results_ledger.md'),Path('docs/results_ledger.json'),
        *sorted(Path('manuscript').rglob('*.md')),*sorted(Path('manuscript').rglob('*.tex')),
        *sorted(Path('manuscript').rglob('*.bib')),Path('manuscript/main.pdf')] if f.is_file()}
    for f,h in preserved_scope.items():pin(f,h)
    pdf_note=p['parent_historical_export_note']['differences'][0]
    assert digest(pdf_note['path'])==pdf_note['current']
    assert digest(pdf_note['original_bytes_preserved_in'])==pdf_note['expected']
    for f,h in list(WATCHED.items()):assert digest(f)==h,f
    result=dict(passed=bool(overall),accepted_at_tested_scope=bool(overall),protocol_sha256=ph,
        counts=dict(backgrounds=168,logical_starts=336,unique_physical_initial_value_problems=len(distinct),distinct_pulse_problems=sum(any(pulse_keys.get((t['context'],t['target']))==t['key'] for t in trials if t['key'] in g) for g in distinct.values()),
            production_paths=1008,cpu_paths=1008,reference_trajectories=672,independent_histories=3,new_histories=0),
        input_audit=prepared,stage_checks=stages,reference_errors=noise,production_accuracy=accuracy,orders=orders,
        raw_timestep_comparisons=rawpairs,response_comparisons=comparisons,duplicate_continuity=duplicates,
        scientific_comparisons=scientific,descriptive_contrasts=descriptive,baseline_descriptions=baseline,
        history_outcomes=outcomes,matching_reference_response_contrasts=reference_contrasts,
        accepted_parent_comparisons=parent_comparisons,
        endpoint_reference_differences=endpoint_reference_differences,control_drift=control_drift,
        amount_evolution=amounts,historical_noninput_export_exception=p['parent_historical_export_note'],
        inputs_unchanged=True,preserved_ledger_manuscript_sha256=preserved_scope,evidence_sha256=WATCHED,recorded_metrics_compared=len(RECORDED_DIFFERENCES),
        maximum_recorded_metric_difference=max(RECORDED_DIFFERENCES)[0],maxima=maxima,gate_counts=gate_counts,
        coordinator_timing=timing,script_sha256=digest(__file__),
        scope='Independent saved-array reconstruction without project assessment or stepping helpers. Shared actual t570 endpoint graphs frozen for immediate recipient reset/pulse; initial local chemistry and surrounding/global species budgets matched in crossed measured-reference arms; later budgets evolve freely. Three existing histories; no moving reset, settled challenge, fresh development, spatial convergence, autonomy or inheritance. Unsaved SSPRK stages remain runtime/source evidence.')
    np.savez_compressed(output.with_suffix('.npz'),**archive);result['waveform_sha256']=digest(output.with_suffix('.npz'))
    output.with_suffix('.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ('passed','counts','gate_counts','maxima','history_outcomes')},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--root',type=Path,default=Path('output/studies/pulse-time-amount-response'))
    parser.add_argument('--output',type=Path,default=Path(__file__).with_suffix(''));args=parser.parse_args();main(args.root,args.output)
