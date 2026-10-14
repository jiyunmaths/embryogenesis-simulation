"""Independent reconstruction from saved arrays; imports no project modules.

No integration, mutable study output, threshold changes, or source modifications.
Temporary review artifacts only. Call after the coordinator has completed.
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


def main(root,output):
    root=Path(root).resolve();output=Path(output)
    p=read_json(root/'protocol.json');ph=WATCHED[str((root/'protocol.json').resolve())]
    for field in ('source_sha256','input_sha256'):
        for file,expected in p[field].items():pin(file,expected)
    status=read_json(root/'status.json');summary=read_json(root/'summary.json')
    if status['state'] not in ('completed','completed_with_unresolved_checks'):raise AssertionError('Coordinator not finished')
    assert summary['protocol_sha256']==ph
    assert (p['beta'],p['da'],p['db'],p['start'],p['duration'],p['interval'],p['late_window'],p['effect_min'])==(2.,.02,.55,510.,60.,.15,24.,.01)
    assert [(x['level'],x['dt']) for x in p['levels']]==[('coarse',.00375),('fine',.001875),('diagnostic',.0009375)]
    expected_criteria=dict(raw_log_max=.01,normalized_response_max=.01,relative_auc_max=.02,recovery_time_max=.3,effect_error_max=.01,reference_log_max=1e-8,reference_normalized_max=1e-8,cpu_torch_log_max=1e-8,order_min=1.7,order_max=2.3,reference_noise_factor=100.,pulse_amount_relative_max=2e-13)
    assert p['criteria']==expected_criteria
    assert p['reference_settings']==dict(rtol=1e-13,atol=1e-15,max_step=.15)
    criteria=p['criteria'];times=np.arange(401)*.15;trials=p['trials'];contexts={c['key']:c for c in p['contexts']}
    assert len(trials)==33 and len(contexts)==9 and len({t['key'] for t in trials})==33
    states={};geometry={};initial_groups=defaultdict(list);controls={};context_checks=[]
    for key,context in contexts.items():
        z=read_npz(context['source']);delta=z['delta'];volumes=z['volumes'];ids=z['ids'];initial=z['initial'];n=len(ids)
        assert n==16 and ids.shape==(n,) and len(set(ids.tolist()))==n and np.all(volumes>0)
        assert delta.shape==(n,n) and initial.dtype==np.float64;valid(initial,(2,n))
        g=volumes[:,None]*delta;off=delta.copy();np.fill_diagonal(off,0.)
        dscale=max(1.,float(np.max(np.abs(delta))));gscale=max(1e-300,float(np.max(np.sum(abs(g),axis=1))))
        assert np.all(off>=0) and np.all(np.diag(delta)<=0)
        rowres=float(abs(delta.sum(axis=1)).max()/dscale)
        massres=float(abs(volumes@delta).max()/gscale)
        reciprocity=float(abs(g-g.T).max()/gscale)
        assert rowres<1e-12 and massres<1e-12 and reciprocity<1e-12
        oldrow=read_json(context['physical_history'])[0]
        assert oldrow['time']==510. and oldrow['elapsed']==0.
        for name,array in [('chemistry',initial),('volumes',volumes),('delta',delta),('ids',ids)]:assert np.array_equal(oldrow[name],array)
        original=read_npz(context['original_source'],['activator','inhibitor','ids','metadata'])
        assert np.array_equal(np.array([original['activator'],original['inhibitor']]),initial)
        assert np.array_equal(original['ids'],ids)
        metadata=json.loads(str(original['metadata']));assert metadata['time']==510.
        assert tuple(metadata['config'][k] for k in ('signal_beta','signal_da','signal_dh'))==(2.,.02,.55)
        q=float(-np.diag(delta).min());subcounts=[max(1,int(np.ceil(l['dt']*max(1+.02*q,2+.55*q)/.2))) for l in p['levels']]
        assert subcounts==[1,1,1]
        geometry[key]=z
        context_checks.append(dict(context=key,rowsum_scaled=rowres,amount_conservation_scaled=massres,reciprocity_scaled=reciprocity,maximum_exit=q,substeps=subcounts))
    pulse_amount_errors=[]
    for trial in trials:
        key=trial['key'];z=geometry[trial['context']];initial=read_npz(trial['source'])['initial'];valid(initial,(2,16));states[key]=initial
        original=read_npz(trial['original_source'],['activator','inhibitor','ids','metadata'])
        assert np.array_equal(original['ids'],z['ids']) and np.array_equal(np.array([original['activator'],original['inhibitor']]),initial)
        assert json.loads(str(original['metadata']))['time']==510.
        expected=z['initial'].copy()
        if trial['target'] is None:
            controls[trial['context']]=key
            assert trial['factor']==1.
        else:
            index=z['ids'].tolist().index(trial['target']);expected[0,index]*=trial['factor']
            removal=float((z['initial'][0,index]-initial[0,index])*z['volumes'][index]);before=float(z['initial'][0,index]*z['volumes'][index])
            assert trial['removed_activator_amount']==removal and trial['initial_target_activator_amount']==before
            error=abs(removal/trial['requested_removed_amount']-1.);pulse_amount_errors.append(error)
            assert error<=criteria['pulse_amount_relative_max']
            assert trial['initial_log_amplitude']==float(abs(np.log(trial['factor'])))
            if trial['dose']=='fractional':assert trial['factor']==.9 and trial['requested_removed_amount']==.1*before
        assert np.array_equal(initial,expected)
        signature=(z['delta'].tobytes(),z['volumes'].tobytes(),z['ids'].tobytes(),initial.tobytes())
        initial_groups[signature].append(key)
    assert len(initial_groups)==27 and p['unique_initial_value_problems']==27 and p['logical_reference_cases']==33
    for seed in p['histories']:
        for cell in p['selections'][str(seed)]['selected_ids']:
            pair=[t for t in trials if t['seed']==seed and t['target']==cell and t['dose']=='amount']
            assert len(pair)==2
            requested=.1*min(t['initial_target_activator_amount'] for t in pair)
            assert all(t['requested_removed_amount']==requested for t in pair)
    paths={};reference_errors={};reference_records={};cpu_records={};production_records={};stage_checks=[]
    for trial in trials:
        key=trial['key'];folder=root/'references'/key;r=read_json(folder/'result.json');z=read_npz(folder/'paths.npz')
        assert r['protocol_sha256']==ph and r['trial']==trial and r['paths_sha256']==WATCHED[str((folder/'paths.npz').resolve())]
        assert np.array_equal(z['times'],times)
        for label in ('dop853','radau'):paths[label,key]=z[label];valid(z[label],(401,2,16));assert np.array_equal(z[label][0],states[key])
        error=log_error(z['dop853'],z['radau']);reference_errors[key]=error;reference_records[key]=r
        compare_record(error,r['raw_solver_log_error'],key+'/reference_error')
        compare_record(error<=criteria['reference_log_max'],r['passed'],key+'/reference_pass')
    mapping=[]
    for ci,context in enumerate(p['contexts']):
        group=[t for t in trials if t['context']==context['key']]
        assert group[0]['target'] is None
        mapping.extend((t['key'],ci,col) for col,t in enumerate(group))
    assert [x[0] for x in mapping]==[t['key'] for t in trials]
    for level in p['levels']:
        label=level['level'];folder=root/'production'/label;r=read_json(folder/'result.json');z=read_npz(folder/'paths.npz');cp=read_npz(folder/'checkpoint.npz');commit=read_json(folder/'checkpoint-commit.json')
        assert r['protocol_sha256']==ph and r['level']==level
        assert r['paths_sha256']==WATCHED[str((folder/'paths.npz').resolve())] and r['checkpoint_sha256']==WATCHED[str((folder/'checkpoint.npz').resolve())]
        assert commit['protocol_sha256']==ph and commit['checkpoint_sha256']==r['checkpoint_sha256']
        assert np.array_equal(z['times'],times) and list(z['trial_keys'])==[m[0] for m in mapping]
        assert z['trajectory'].shape==(401,33,2,16);valid(z['trajectory']);assert np.array_equal(z['trajectory'],cp['history'])
        meta=json.loads(str(cp['metadata']));assert meta['step']==int(round(60/level['dt'])) and meta['protocol_sha256']==ph and meta['level']==level
        assert commit['step']==meta['step'] and commit['level']==level
        valid(cp['state'],(2,144,5));padding_error=0.
        for i,(key,ci,col) in enumerate(mapping):
            path=z['trajectory'][:,i];paths[label,key]=path
            assert np.array_equal(path[0],states[key]) and np.array_equal(path[-1],cp['state'][:,16*ci:16*(ci+1),col])
        for ci,context in enumerate(p['contexts']):
            group=[t for t in trials if t['context']==context['key']]
            control=cp['state'][:,16*ci:16*(ci+1),0]
            for col in range(len(group),5):padding_error=max(padding_error,log_error(control,cp['state'][:,16*ci:16*(ci+1),col]))
        assert r['actual_paths']==33 and r['padding_paths']==12
        stage_checks.append(dict(level=label,paths=33,samples=401,outer_steps=meta['step'],positive_saved_concentrations=True,final_checkpoint_matches_paths=True,padding_log_error=padding_error,wall_seconds=r['wall_seconds']))
        production_records[label]=r
        for trial in trials:
            key=trial['key'];folder=root/'cpu-verification'/label/key;r=read_json(folder/'result.json');z=read_npz(folder/'paths.npz')
            assert r['protocol_sha256']==ph and r['trial']==trial and r['level']==level and r['paths_sha256']==WATCHED[str((folder/'paths.npz').resolve())]
            assert np.array_equal(z['times'],times);valid(z['trajectory'],(401,2,16));assert np.array_equal(z['trajectory'][0],states[key])
            paths['cpu-'+label,key]=z['trajectory'];cpu_records[label,key]=r
    labels=['dop853','radau','coarse','fine','diagnostic'];waves={};changes={};metrics={};effects={}
    def calculate_label(label):
        for trial in trials:
            if trial['target'] is None:continue
            key=trial['key'];geo=geometry[trial['context']];control=controls[trial['context']];index=geo['ids'].tolist().index(trial['target'])
            wave,change,metric=response(paths[label,key],paths[label,control],times,geo['volumes'],index,trial['factor'],p['late_window'])
            waves[label,key]=wave;changes[label,key]=change;metrics[label,key]=metric
        for seed in p['histories']:
            for target in p['selections'][str(seed)]['selected_ids']:
                for dose in ('fractional','amount'):
                    sham=next(t for t in trials if (t['seed'],t['background'],t['target'],t['dose'])==(seed,'sham',target,dose))
                    red=next(t for t in trials if (t['seed'],t['background'],t['target'],t['dose'])==(seed,f'redistributed-cell-{target}',target,dose))
                    effects[label,seed,target,dose]=rms(waves[label,red['key']]-waves[label,sham['key']],times)
    for label in labels:calculate_label(label)
    for label in labels:
        recorded=summary['response_summaries'][label]
        for trial in trials:
            if trial['target'] is None:continue
            for metric,value in metrics[label,trial['key']].items():compare_record(value,recorded['metrics'][trial['key']][metric],f'{label}/{trial["key"]}/{metric}')
        for row in recorded['effects']:
            value=effects[label,row['seed'],row['cell'],row['dose']]
            compare_record(value,row['normalized_response_shift_rms'],f'{label}/effect/{row["seed"]}/{row["cell"]}/{row["dose"]}')
            compare_record(value>p['effect_min'],row['response_shift_detected'],f'{label}/effect_decision')
    comparison_specs=[('dop853','radau',criteria['reference_normalized_max']),('coarse','fine',criteria['normalized_response_max']),('fine','diagnostic',criteria['normalized_response_max'])]
    comparison_specs += [(level['level'],method,criteria['normalized_response_max']) for level in p['levels'] for method in ('dop853','radau')]
    comparisons=[]
    for a,b,limit in comparison_specs:
        rows=[];effect_rows=[]
        for trial in trials:
            if trial['target'] is None:continue
            key=trial['key'];error=float(abs(changes[a,key]-changes[b,key]).max()/abs(np.log(trial['factor'])))
            left,right=metrics[a,key],metrics[b,key];auc=abs(right['target_activator_log_auc_per_log_pulse']/left['target_activator_log_auc_per_log_pulse']-1.)
            rec={}
            for name in ('target_recovery_time','network_recovery_time'):
                x,y=left[name],right[name]
                rec[name]=0. if x is None and y is None else None if x is None or y is None else abs(x-y)
            passed=error<=limit and auc<=criteria['relative_auc_max'] and all(v is not None and v<=criteria['recovery_time_max']+1e-9 for v in rec.values())
            rows.append(dict(trial=key,normalized_response_error=error,relative_auc_error=float(auc),recovery_time_errors=rec,passed=bool(passed)))
        for seed in p['histories']:
            for cell in p['selections'][str(seed)]['selected_ids']:
                for dose in ('fractional','amount'):
                    x,y=effects[a,seed,cell,dose],effects[b,seed,cell,dose];same=(x>p['effect_min'])==(y>p['effect_min']);error=abs(x-y)
                    effect_rows.append(dict(seed=seed,cell=cell,dose=dose,rms_error=error,same_decision=bool(same),passed=bool(same and error<=criteria['effect_error_max'])))
        result=dict(first=a,second=b,pulse_rows=rows,effects=effect_rows,passed=all(x['passed'] for x in rows+effect_rows));comparisons.append(result)
        recorded=next(r for r in summary['response_comparisons'] if r['first']==a and r['second']==b)
        compare_record(result['passed'],recorded['passed'],a+'/'+b+'/pass')
        for row in rows:
            saved=next(r for r in recorded['pulse_rows'] if r['trial']==row['trial'])
            for k in ('normalized_response_error','relative_auc_error','passed'):compare_record(row[k],saved[k],a+'/'+b+'/'+row['trial']+'/'+k)
            for k,v in row['recovery_time_errors'].items():compare_record(v,saved['recovery_time_errors'][k],a+'/'+b+'/'+row['trial']+'/'+k)
        for row in effect_rows:
            saved=next(r for r in recorded['effects'] if (r['seed'],r['cell'],r['dose'])==(row['seed'],row['cell'],row['dose']))
            for k in ('rms_error','same_decision','passed'):compare_record(row[k],saved[k],a+'/'+b+'/effect/'+k)
    accuracy=[];orders=[];rawpairs=[]
    for trial in trials:
        key=trial['key'];values=[]
        for level in p['levels']:
            label=level['level'];error=max(log_error(paths[label,key],paths[method,key]) for method in ('dop853','radau'));cpu=log_error(paths[label,key],paths['cpu-'+label,key]);values.append(error)
            row=dict(trial=key,level=label,raw_reference_log_error=error,cpu_torch_log_error=cpu,passed=error<=criteria['raw_log_max'] and cpu<=criteria['cpu_torch_log_max']);accuracy.append(row)
            saved=next(r for r in summary['production_accuracy'] if r['trial']==key and r['level']==label)
            for k in ('raw_reference_log_error','cpu_torch_log_error','passed'):compare_record(row[k],saved[k],key+'/'+label+'/'+k)
        noise=reference_errors[key];resolved=values[-1]>criteria['reference_noise_factor']*noise
        oo=np.log2(np.array(values[:-1])/values[1:]).tolist() if resolved and all(v>0 for v in values) else None
        passed=not resolved or (oo is not None and all(criteria['order_min']<=v<=criteria['order_max'] for v in oo))
        row=dict(trial=key,errors=values,reference_noise=noise,order_resolved=bool(resolved),orders=oo,passed=bool(passed));orders.append(row)
        saved=next(r for r in summary['orders'] if r['trial']==key)
        for k in ('reference_noise','order_resolved','passed'):compare_record(row[k],saved[k],key+'/order/'+k)
        if oo is None:assert saved['orders'] is None
        else:
            for x,y in zip(oo,saved['orders']):compare_record(x,y,key+'/order')
        for a,b in [('coarse','fine'),('fine','diagnostic')]:
            error=log_error(paths[a,key],paths[b,key]);row=dict(trial=key,first=a,second=b,raw_log_error=error,passed=error<=criteria['raw_log_max']);rawpairs.append(row)
            saved=next(r for r in summary['raw_timestep_comparisons'] if r['trial']==key and r['first']==a and r['second']==b)
            for k in ('raw_log_error','passed'):compare_record(row[k],saved[k],key+'/'+a+'/'+b+'/'+k)
    refqual=read_json(root/'reference-qualification.json')
    assert refqual['protocol_sha256']==ph
    for file,expected in refqual['evidence_sha256'].items():pin(file,expected)
    references_pass=all(x<=criteria['reference_log_max'] for x in reference_errors.values()) and comparisons[0]['passed']
    compare_record(references_pass,refqual['passed'],'reference_qualification')
    for file,expected in summary['evidence_sha256'].items():pin(file,expected)
    # Matched moving data are reconstructed separately from their raw archived histories.
    parent=Path(p['parent']);old=read_json(parent/'protocol.json');review=read_json(p['review_file'])
    assert digest(parent/'protocol.json')==p['parent_protocol_sha256']==review['protocol_sha256']
    archived=read_json(Path(review['review_archive'])/'provenance_verification.json')['evidence_sha256']
    moving_clock_max=0.;moving_start_checks=0
    for level in ('coarse','fine'):
        for trial in trials:
            job=next(j for j in old['jobs'] if (j['seed'],j['level'],j['background'],j['target'],j['dose'])==(trial['seed'],level,trial['background'],trial['target'],trial['dose']))
            item=old['reused'].get(job['key']);folder=Path(old['reuse_parent'])/item['job']['key'] if item else parent/job['key'];file=folder/'history.json'
            expected=archived.get(str(file),archived.get(str(file.resolve())));assert expected is not None
            pin(file,expected);history=read_json(file)
            assert len(history)==401
            mt=np.array([r['elapsed'] for r in history]);clock_error=float(abs(mt-times).max());moving_clock_max=max(moving_clock_max,clock_error)
            assert clock_error<1e-10 and abs(history[0]['time']-510.)<1e-10 and abs(history[-1]['time']-570.)<1e-8
            path=np.array([r['chemistry'] for r in history]);valid(path,(401,2,16));assert np.array_equal(path[0],states[trial['key']])
            geo=geometry[trial['context']]
            assert all(r['ids']==geo['ids'].tolist() for r in history)
            assert np.array_equal(history[0]['delta'],geo['delta']) and np.array_equal(history[0]['volumes'],geo['volumes']) and job['factor']==trial['factor']
            moving_start_checks+=1;paths['moving-'+level,trial['key']]=path
        calculate_label('moving-'+level)
    for recorded in review['recipients']:
        seed,cell=recorded['history'],recorded['cell']
        for level in ('coarse','fine'):
            for dose in ('fractional','amount'):
                value=effects['moving-'+level,seed,cell,dose]
                compare_record(value,recorded['by_timestep'][level][dose]['normalized_response_shift_rms'],f'moving/{level}/{seed}/{cell}/{dose}')
    scientific=[];baseline=[];archive_arrays={'times':times};early=times<=6.;late=times>=36.
    for seed in p['histories']:
        for cell in p['selections'][str(seed)]['selected_ids']:
            shcontext=next(c for c in p['contexts'] if c['seed']==seed and c['background']=='sham')['key'];rdcontext=next(c for c in p['contexts'] if c['seed']==seed and c['background']==f'redistributed-cell-{cell}')['key']
            shcontrol,rdcontrol=controls[shcontext],controls[rdcontext];index=geometry[shcontext]['ids'].tolist().index(cell)
            base=dict(seed=seed,cell=cell,description='Late baseline state log RMS is a two-species target time RMS over elapsed36..60; final ratios are concentration ratios. No new scientific gate.')
            for label in ('fine','moving-fine'):
                sh,rd=paths[label,shcontrol],paths[label,rdcontrol];contrast=np.log(rd[:,:,index]/sh[:,:,index])
                base[label]=dict(final_activator_redistribution_sham=float(rd[-1,0,index]/sh[-1,0,index]),final_inhibitor_redistribution_sham=float(rd[-1,1,index]/sh[-1,1,index]),late_target_log_rms=rms(contrast[late],times[late]),sham_contrast_min=float(np.std(np.log(sh[:,0]),axis=1).min()),redistributed_contrast_min=float(np.std(np.log(rd[:,0]),axis=1).min()),sham_late_contrast_min=float(np.std(np.log(sh[late,0]),axis=1).min()),redistributed_late_contrast_min=float(np.std(np.log(rd[late,0]),axis=1).min()),sham_final_amounts=(sh[-1]@geometry[shcontext]['volumes']).tolist() if label=='fine' else None,redistributed_final_amounts=(rd[-1]@geometry[rdcontext]['volumes']).tolist() if label=='fine' else None)
            baseline.append(base)
            for dose in ('fractional','amount'):
                sh=next(t for t in trials if t['context']==shcontext and t['target']==cell and t['dose']==dose);rd=next(t for t in trials if t['context']==rdcontext and t['target']==cell and t['dose']==dose)
                row=dict(seed=seed,cell=cell,dose=dose,full_effects={label:effects[label,seed,cell,dose] for label in ('coarse','fine','diagnostic','dop853','radau','moving-coarse','moving-fine')})
                row['frozen_detected']=row['full_effects']['fine']>p['effect_min'];row['moving_detected']=row['full_effects']['moving-fine']>p['effect_min'];row['frozen_minus_moving_effect']=row['full_effects']['fine']-row['full_effects']['moving-fine']
                differences={label:waves[label,rd['key']]-waves[label,sh['key']] for label in ('fine','moving-fine')}
                row['signed_response_contrast_waveform_change_rms']=rms(differences['fine']-differences['moving-fine'],times)
                row['sham_response_waveform_change_rms']=rms(waves['fine',sh['key']]-waves['moving-fine',sh['key']],times)
                row['redistributed_response_waveform_change_rms']=rms(waves['fine',rd['key']]-waves['moving-fine',rd['key']],times)
                row['early_descriptive']={}
                for label,dd in differences.items():
                    total=float(trapezoid(np.mean(dd**2,axis=1),times));early_area=float(trapezoid(np.mean(dd[early]**2,axis=1),times[early]))
                    row['early_descriptive'][label]=dict(first6_rms=rms(dd[early],times[early]),first6_fraction_of_integrated_squared_difference=early_area/total if total else None)
                    archive_arrays[f'seed{seed}_cell{cell}_{dose}_{label}_sham']=waves[label,sh['key']]
                    archive_arrays[f'seed{seed}_cell{cell}_{dose}_{label}_redistributed']=waves[label,rd['key']]
                row['fine_pulse_metrics']={background:metrics['fine',t['key']] for background,t in [('sham',sh),('redistributed',rd)]}
                scientific.append(row)
    duplicate_groups=[keys for keys in initial_groups.values() if len(keys)>1];duplicate_error=0.
    for keys in duplicate_groups:
        for label in labels:
            for key in keys[1:]:duplicate_error=max(duplicate_error,log_error(paths[label,keys[0]],paths[label,key]))
    history_outcomes=[]
    for seed in p['histories']:
        row=dict(seed=seed,recipients=2)
        for label in ('fine','moving-fine'):
            row[label]={dose:sum(effects[label,seed,cell,dose]>p['effect_min'] for cell in p['selections'][str(seed)]['selected_ids']) for dose in ('fractional','amount')}
            row[label]['both_doses']=sum(all(effects[label,seed,cell,dose]>p['effect_min'] for dose in ('fractional','amount')) for cell in p['selections'][str(seed)]['selected_ids'])
        saved=next(r for r in summary['history_outcomes'] if r['seed']==seed)
        compare_record(row['fine']['both_doses'],saved['both_doses_above_screen'],f'history{seed}/count')
        history_outcomes.append(row)
    all_pass=references_pass and all(r['passed'] for r in accuracy+orders+rawpairs+comparisons)
    compare_record(all_pass,summary['passed'],'overall_pass');compare_record(all_pass,status['passed'],'status_pass')
    # Compare low-dose trajectories against the separately coded, pre-freeze golden reference.
    golden_file=Path('/tmp/frozen_response_reference_probe_tight.npz');golden=None
    if golden_file.exists():
        golden_record=read_json(golden_file.with_suffix('.json'));pin(golden_file,golden_record['arrays_sha256'])
        g=read_npz(golden_file);golden={}
        assert np.array_equal(g['times'],times)
        for seed,cell in ((8,20),(9,30)):
            for dose in (None,'amount'):
                trial=next(t for t in trials if (t['seed'],t['background'],t['target'],t['dose'])==(seed,'sham',None if dose is None else cell,dose))
                short=f'seed-{seed}_'+('control' if dose is None else f'amount-cell-{cell}')
                for label,method in [('dop853','DOP853'),('radau','Radau')]:golden[short+'_'+method]=log_error(paths[label,trial['key']],g[short+'_'+method])
    for name,expected in list(WATCHED.items()):
        if digest(name)!=expected:raise AssertionError('Changed input during review: '+name)
    resolved=[v for r in orders if r['orders'] is not None for v in r['orders']]
    results=dict(scope='Independent array-level audit; no embryo module imports, integrations, project writes, or threshold changes. Actual unsaved internal Euler stages cannot be reconstructed from observations; their successful checks are runtime records/code guarantees, distinct from independently reconstructed sampled positivity and final checkpoint consistency.',protocol_sha256=ph,passed=bool(all_pass),counts=dict(contexts=9,logical_reference_cases=33,unique_initial_value_problems=27,production_paths=99,cpu_verification_paths=99,reference_trajectories=66,paired_moving_paths=66,independent_histories=3,new_histories=0),criterion_values=criteria,
        reference_pass=bool(references_pass),context_checks=context_checks,pulse_amount_relative_max=max(pulse_amount_errors),duplicate_groups=duplicate_groups,duplicate_path_log_error=duplicate_error,stage_checks=stage_checks,
        gate_counts=dict(reference_raw_pass=sum(x<=criteria['reference_log_max'] for x in reference_errors.values()),production_accuracy_pass=sum(r['passed'] for r in accuracy),order_pass=sum(r['passed'] for r in orders),resolved_orders=sum(r['order_resolved'] for r in orders),raw_timestep_pass=sum(r['passed'] for r in rawpairs),response_comparison_pass=sum(r['passed'] for r in comparisons),response_pulse_rows_pass=sum(r['passed'] for c in comparisons for r in c['pulse_rows']),effect_rows_pass=sum(r['passed'] for c in comparisons for r in c['effects'])),
        maxima=dict(reference_raw=max(reference_errors.values()),reference_normalized=max(r['normalized_response_error'] for r in comparisons[0]['pulse_rows']),production_raw_reference=max(r['raw_reference_log_error'] for r in accuracy),cpu_torch_raw=max(r['cpu_torch_log_error'] for r in accuracy),coarse_fine_raw=max(r['raw_log_error'] for r in rawpairs if r['first']=='coarse'),fine_diagnostic_raw=max(r['raw_log_error'] for r in rawpairs if r['first']=='fine'),order_min=min(resolved) if resolved else None,order_max=max(resolved) if resolved else None,normalized_response=max(r['normalized_response_error'] for c in comparisons[1:] for r in c['pulse_rows']),relative_auc=max(r['relative_auc_error'] for c in comparisons[1:] for r in c['pulse_rows']),recovery_time=max(v for c in comparisons for r in c['pulse_rows'] for v in r['recovery_time_errors'].values() if v is not None),effect_discrepancy=max(r['rms_error'] for c in comparisons for r in c['effects'])),
        moving_clock_max_absolute_error=moving_clock_max,moving_initial_sources_exact=moving_start_checks,
        scientific_comparisons=scientific,baseline_descriptions=baseline,history_outcomes=history_outcomes,production_accuracy=accuracy,orders=orders,raw_timestep_comparisons=rawpairs,response_comparisons=comparisons,
        golden_reference_raw_log_errors=golden,reference_errors=reference_errors,
        sampled_min_concentration=min(float(v.min()) for v in paths.values()),fine_all_paths_contrast_min=min(float(np.std(np.log(paths['fine',t['key']][:,0]),axis=1).min()) for t in trials),
        recorded_metrics_compared=len(RECORDED_DIFFERENCES),max_recorded_metric_difference=max(RECORDED_DIFFERENCES)[0],max_recorded_metric_difference_location=max(RECORDED_DIFFERENCES)[1],inputs_unchanged=True,inputs_sha256=WATCHED,script_sha256=digest(__file__))
    np.savez_compressed(output.with_suffix('.npz'),**archive_arrays);results['waveforms_sha256']=digest(output.with_suffix('.npz'))
    output.with_suffix('.json').write_text(json.dumps(results,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:results[k] for k in ('passed','counts','gate_counts','maxima','history_outcomes','golden_reference_raw_log_errors','inputs_unchanged')},indent=2))
    for r in scientific:print(r['seed'],r['cell'],r['dose'],'frozen',r['full_effects']['fine'],'moving',r['full_effects']['moving-fine'],'difference',r['frozen_minus_moving_effect'])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path('/home/sji/Documents/Programming/Python/embryogenesis-simulation/output/studies/phase-carry-frozen-response'))
    parser.add_argument('--output',type=Path,default=Path('/tmp/frozen_response_independent_review'))
    args=parser.parse_args();main(args.root,args.output)
