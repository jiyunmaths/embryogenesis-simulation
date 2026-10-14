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
    p=read_json(root/'protocol.json');ph=digest(root/'protocol.json')
    for field in ('source_sha256','input_sha256'):
        for file,expected in p[field].items():pin(file,expected)
    status=read_json(root/'status.json');summary=read_json(root/'summary.json')
    assert status['state'] in ('completed','completed_with_unresolved_checks')
    assert summary['protocol_sha256']==ph
    assert tuple(p[k] for k in ('beta','da','db','start','duration','interval','late_window','effect_min'))==(2.,.02,.55,510.,60.,.15,24.,.01)
    assert [(x['level'],x['dt']) for x in p['levels']]==[('coarse',.00375),('fine',.001875),('diagnostic',.0009375)]
    assert p['reference_settings']==dict(rtol=1e-13,atol=1e-15,max_step=.15)
    inherited=dict(raw_log_max=.01,normalized_response_max=.01,relative_auc_max=.02,recovery_time_max=.3,effect_error_max=.01,reference_log_max=1e-8,reference_normalized_max=1e-8,cpu_torch_log_max=1e-8,order_min=1.7,order_max=2.3,reference_noise_factor=100.,pulse_amount_relative_max=2e-13)
    for k,v in inherited.items():assert p['criteria'][k]==v,(k,p['criteria'][k])
    criteria=p['criteria'];times=np.arange(401)*.15;early=times<=6.;late=times>=36.
    parent=Path(p['parent']);old=read_json(parent/'protocol.json');oldsummary=read_json(parent/'summary.json')
    assert oldsummary['passed'] and read_json(parent/'status.json')['state']=='completed'
    for oldfile,expected in oldsummary['evidence_sha256'].items():pin(oldfile,expected)
    for field in ('source_sha256','input_sha256'):
        for file,expected in old[field].items():pin(file,expected)
    donors={}
    for c in old['contexts']:
        z=read_npz(c['source']);donors[c['key']]={**z,'context':c}
    assert len(donors)==9
    contexts={};source_checks=[];transplants=[];controls={};states={};cases={};groups=defaultdict(list)
    for c in p['contexts']:
        z=read_npz(c['source']);valid(z['initial'],(2,16))
        assert z['delta'].shape==(16,16) and z['volumes'].shape==(16,) and np.all(z['volumes']>0)
        candidates=[d for d in donors.values() if d['context']['seed']==c['seed']]
        chemistry=[d for d in candidates if np.array_equal(d['initial'],z['initial']) and np.array_equal(d['ids'],z['ids'])]
        physics=[d for d in candidates if all(np.array_equal(d[k],z[k]) for k in ('delta','volumes','ids'))]
        assert len(chemistry)==len(physics)==1,'Prepared context lacks unique exact chemical/physical donor'
        cd,pd=chemistry[0],physics[0];cb,pb=cd['context']['background'],pd['context']['background']
        target_set=set(old['selections'][str(c['seed'])]['selected_ids'])
        ct=None if cb=='sham' else int(cb.rsplit('-',1)[1]);pt=None if pb=='sham' else int(pb.rsplit('-',1)[1])
        assert (ct is None or ct in target_set) and (pt is None or pt in target_set)
        assert ct is None or pt is None or ct==pt,'Wrong-recipient factorial crossing'
        expected_targets=sorted(target_set if ct is None and pt is None else {ct if ct is not None else pt})
        assert sorted(c['targets'])==expected_targets
        assert c['chemical_context']==cd['context']['key'] and c['physical_context']==pd['context']['key']
        z.update(context=c,chemical_donor=cd['context']['key'],physical_donor=pd['context']['key'],C='S' if ct is None else 'R',P='S' if pt is None else 'R',recipient=ct if ct is not None else pt)
        contexts[c['key']]=z
        delta=z['delta'];v=z['volumes'];off=delta.copy();np.fill_diagonal(off,0.);g=v[:,None]*delta
        ds=max(1.,float(abs(delta).max()));gs=max(1e-300,float(np.sum(abs(g),axis=1).max()))
        row=float(abs(delta.sum(axis=1)).max()/ds);mass=float(abs(v@delta).max()/gs);reciprocity=float(abs(g-g.T).max()/gs)
        assert np.all(off>=0) and np.all(np.diag(delta)<=0) and row<1e-12 and mass<1e-12 and reciprocity<1e-12
        q=float(-np.diag(delta).min());substeps=[max(1,int(np.ceil(l['dt']*max(1+.02*q,2+.55*q)/.2))) for l in p['levels']];assert substeps==[1,1,1]
        source_checks.append(dict(context=c['key'],chemical_source=cd['context']['key'],physical_source=pd['context']['key'],C=z['C'],P=z['P'],recipients=expected_targets,maximum_exit=q,row_sum_scaled=row,amount_conservation_scaled=mass,reciprocity_scaled=reciprocity,substeps=substeps))
        before=cd['initial']*cd['volumes'][None,:];after=z['initial']*v[None,:]
        transplants.append(dict(context=c['key'],chemical_source=z['chemical_donor'],physical_source=z['physical_donor'],original_per_cell_amounts=before.tolist(),reassigned_per_cell_amounts=after.tolist(),per_cell_amount_change=(after-before).tolist(),original_total_amounts=before.sum(axis=1).tolist(),reassigned_total_amounts=after.sum(axis=1).tolist(),total_amount_change=(after-before).sum(axis=1).tolist(),total_amount_relative_change=(after.sum(axis=1)/before.sum(axis=1)-1.).tolist(),volume_ratio=(v/cd['volumes']).tolist()))
    assert len(contexts)==21
    assert len(p['factorials'])==6
    for factorial in p['factorials']:
        assert set(factorial['contexts'])=={'S_S','S_R','R_S','R_R'}
        for name,key in factorial['contexts'].items():
            z=contexts[key]
            assert name==z['C']+'_'+z['P'] and factorial['seed']==z['context']['seed'] and factorial['cell'] in z['context']['targets']
    trials=p['trials'];assert len(trials)==69 and len({t['key'] for t in trials})==69
    dose_errors=[];dose_details=[]
    for trial in trials:
        key=trial['key'];z=contexts[trial['context']];initial=read_npz(trial['source'])['initial'];valid(initial,(2,16));expected=z['initial'].copy()
        if trial['target'] is None:
            assert trial['context'] not in controls;controls[trial['context']]=key;assert trial['factor']==1.
        else:
            cell=trial['target'];seed=trial['seed'];assert cell in z['context']['targets']
            i=z['ids'].tolist().index(cell)
            s=next(d for d in donors.values() if d['context']['seed']==seed and d['context']['background']=='sham')
            r=next(d for d in donors.values() if d['context']['seed']==seed and d['context']['background']==f'redistributed-cell-{cell}')
            minimum=min(s['initial'][0,i],r['initial'][0,i]);current=z['initial'][0,i]
            factor=.9 if trial['dose']=='fractional' else float(1.-.1*minimum/current)
            assert trial['factor']==factor
            expected[0,i]*=factor
            requested=float(.1*current*z['volumes'][i]) if trial['dose']=='fractional' else float(.1*z['volumes'][i]*minimum)
            measured=float((z['initial'][0,i]-initial[0,i])*z['volumes'][i]);error=abs(measured/requested-1.);dose_errors.append(error);assert error<=2e-13
            assert trial['initial_log_amplitude']==float(abs(np.log(factor)))
            # Requested amounts may differ by one arithmetic-rounding unit but cannot hide a dose change.
            for field in ('requested_removed_amount','removed_activator_amount'):
                if field in trial:assert abs(trial[field]/(requested if field=='requested_removed_amount' else measured)-1.)<=2e-13
            case=(seed,cell,z['C'],z['P'],trial['dose']);assert case not in cases;cases[case]=key
            dose_details.append(dict(trial=key,seed=seed,cell=cell,C=z['C'],P=z['P'],dose=trial['dose'],factor=factor,log_amplitude=float(abs(np.log(factor))),requested_amount=requested,measured_amount=measured,relative_error=error))
        assert np.array_equal(initial,expected);states[key]=initial
        sig=(z['delta'].tobytes(),z['volumes'].tobytes(),z['ids'].tobytes(),initial.tobytes());groups[sig].append(key)
    assert len(controls)==21 and len(cases)==48 and len(groups)==57
    for seed in p['histories']:
        for cell in p['selections'][str(seed)]['selected_ids']:
            for C in ('S','R'):
                fs=[next(t['factor'] for t in trials if t['key']==cases[seed,cell,C,P,'amount']) for P in ('S','R')]
                assert fs[0]==fs[1]
            for P in ('S','R'):
                doses=[next(x['measured_amount'] for x in dose_details if (x['seed'],x['cell'],x['C'],x['P'],x['dose'])==(seed,cell,C,P,'amount')) for C in ('S','R')]
                assert abs(doses[0]/doses[1]-1.)<=2e-13
    paths={};references={};stages=[]
    for trial in trials:
        key=trial['key'];folder=root/'references'/key;r=read_json(folder/'result.json');z=read_npz(folder/'paths.npz')
        assert r['protocol_sha256']==ph and r['trial']==trial and r['paths_sha256']==digest(folder/'paths.npz') and np.array_equal(z['times'],times)
        for label in ('dop853','radau'):valid(z[label],(401,2,16));assert np.array_equal(z[label][0],states[key]);paths[label,key]=z[label]
        references[key]=log_error(z['dop853'],z['radau']);compare_record(references[key],r['raw_solver_log_error'],key+'/reference')
        compare_record(references[key]<=1e-8,r['passed'],key+'/reference-pass')
    width=max(sum(t['context']==c['key'] for t in trials) for c in p['contexts']);mapping=[]
    for ci,c in enumerate(p['contexts']):
        group=[t for t in trials if t['context']==c['key']];assert group[0]['target'] is None
        mapping.extend((t['key'],ci,col) for col,t in enumerate(group))
    assert [m[0] for m in mapping]==[t['key'] for t in trials]
    for level in p['levels']:
        label=level['level'];folder=root/'production'/label;r=read_json(folder/'result.json');z=read_npz(folder/'paths.npz');cp=read_npz(folder/'checkpoint.npz');commit=read_json(folder/'checkpoint-commit.json')
        assert r['protocol_sha256']==ph and r['level']==level and r['paths_sha256']==digest(folder/'paths.npz') and r['checkpoint_sha256']==digest(folder/'checkpoint.npz')
        assert commit['protocol_sha256']==ph and commit['checkpoint_sha256']==r['checkpoint_sha256']
        assert list(z['trial_keys'])==[m[0] for m in mapping] and np.array_equal(z['times'],times)
        valid(z['trajectory'],(401,69,2,16));assert np.array_equal(z['trajectory'],cp['history']);valid(cp['state'],(2,16*len(contexts),width))
        meta=json.loads(str(cp['metadata']));assert meta['step']==int(round(60/level['dt'])) and meta['level']==level and meta['protocol_sha256']==ph and commit['step']==meta['step']
        for index,(key,ci,col) in enumerate(mapping):
            path=z['trajectory'][:,index];assert np.array_equal(path[0],states[key]) and np.array_equal(path[-1],cp['state'][:,16*ci:16*(ci+1),col]);paths[label,key]=path
        padding_error=0.
        for ci,c in enumerate(p['contexts']):
            count=sum(t['context']==c['key'] for t in trials)
            for col in range(count,width):padding_error=max(padding_error,log_error(cp['state'][:,16*ci:16*(ci+1),0],cp['state'][:,16*ci:16*(ci+1),col]))
        stages.append(dict(level=label,paths=69,samples=401,outer_steps=meta['step'],padding_paths=len(contexts)*width-69,padding_error=padding_error,wall_seconds=r['wall_seconds'],positive_saved_states=True,final_checkpoint_consistent=True))
        for trial in trials:
            key=trial['key'];folder=root/'cpu-verification'/label/key;r=read_json(folder/'result.json');z=read_npz(folder/'paths.npz')
            assert r['protocol_sha256']==ph and r['trial']==trial and r['level']==level and r['paths_sha256']==digest(folder/'paths.npz') and np.array_equal(z['times'],times)
            valid(z['trajectory'],(401,2,16));assert np.array_equal(z['trajectory'][0],states[key]);paths['cpu-'+label,key]=z['trajectory']
    labels=['dop853','radau','coarse','fine','diagnostic'];waves={};changes={};metrics={};contrasts={};norms={};native={};identity_errors=[]
    for label in labels:
        for trial in trials:
            if trial['target'] is None:continue
            key=trial['key'];z=contexts[trial['context']];index=z['ids'].tolist().index(trial['target']);control=controls[trial['context']]
            waves[label,key],changes[label,key],metrics[label,key]=response(paths[label,key],paths[label,control],times,z['volumes'],index,trial['factor'],24.)
        for seed in p['histories']:
            for cell in p['selections'][str(seed)]['selected_ids']:
                for dose in ('fractional','amount'):
                    wave={(C,P):waves[label,cases[seed,cell,C,P,dose]] for C in ('S','R') for P in ('S','R')}
                    cc={P:wave['R',P]-wave['S',P] for P in ('S','R')};pp={C:wave[C,'R']-wave[C,'S'] for C in ('S','R')};interaction=wave['R','R']-wave['S','R']-wave['R','S']+wave['S','S'];diagonal=wave['R','R']-wave['S','S']
                    for kind,held,value in [('chemical',P,cc[P]) for P in ('S','R')]+[('physical',C,pp[C]) for C in ('S','R')]+[('interaction',None,interaction)]:
                        key=(label,seed,cell,dose,kind,held);contrasts[key]=value;norms[key]=rms(value,times)
                    native[label,seed,cell,dose]=dict(wave=diagonal,norm=rms(diagonal,times))
                    identity_errors.extend(float(abs(value).max()) for value in (diagonal-cc['S']-pp['R'],diagonal-cc['R']-pp['S'],interaction-(cc['R']-cc['S']),interaction-(pp['R']-pp['S'])))
    comparison_specs=[('dop853','radau'),('coarse','fine'),('fine','diagnostic')]+[(l['level'],r) for l in p['levels'] for r in ('dop853','radau')]
    comparisons=[]
    for a,b in comparison_specs:
        limit=1e-8 if (a,b)==('dop853','radau') else .01;pulse_rows=[];contrast_rows=[]
        for trial in trials:
            if trial['target'] is None:continue
            key=trial['key'];error=float(abs(changes[a,key]-changes[b,key]).max()/abs(np.log(trial['factor'])));ma,mb=metrics[a,key],metrics[b,key];auc=abs(mb['target_activator_log_auc_per_log_pulse']/ma['target_activator_log_auc_per_log_pulse']-1.)
            rr={}
            for name in ('target_recovery_time','network_recovery_time'):
                x,y=ma[name],mb[name];rr[name]=0. if x is None and y is None else None if x is None or y is None else abs(x-y)
            passed=error<=limit and auc<=.02 and all(v is not None and v<=.3+1e-9 for v in rr.values())
            pulse_rows.append(dict(trial=key,normalized_response_error=error,relative_auc_error=float(auc),recovery_time_errors=rr,passed=bool(passed)))
        for key,wave in contrasts.items():
            if key[0]!=a:continue
            other=(b,*key[1:]);error=float(abs(wave-contrasts[other]).max());nd=abs(norms[key]-norms[other]);same=(norms[key]>.01)==(norms[other]>.01) if key[-2]=='chemical' else None
            contrast_rows.append(dict(seed=key[1],cell=key[2],dose=key[3],kind=key[4],held=key[5],maximum_waveform_error=error,rms_error=nd,same_primary_decision=same,passed=bool(error<=limit and nd<=limit and same is not False)))
        comparisons.append(dict(first=a,second=b,pulse_rows=pulse_rows,contrast_rows=contrast_rows,passed=all(x['passed'] for x in pulse_rows+contrast_rows)))
    accuracy=[];orders=[];rawpairs=[]
    for trial in trials:
        key=trial['key'];values=[]
        for level in p['levels']:
            label=level['level'];error=max(log_error(paths[label,key],paths[method,key]) for method in ('dop853','radau'));cpu=log_error(paths[label,key],paths['cpu-'+label,key]);values.append(error)
            accuracy.append(dict(trial=key,level=label,raw_reference_error=error,cpu_torch_error=cpu,passed=error<=.01 and cpu<=1e-8))
        noise=references[key];resolved=values[-1]>100*noise;oo=np.log2(np.array(values[:-1])/values[1:]).tolist() if resolved and all(v>0 for v in values) else None
        orders.append(dict(trial=key,errors=values,reference_disagreement=noise,order_resolved=bool(resolved),orders=oo,passed=bool(not resolved or oo is not None and all(1.7<=v<=2.3 for v in oo))))
        for a,b in [('coarse','fine'),('fine','diagnostic')]:
            error=log_error(paths[a,key],paths[b,key]);rawpairs.append(dict(trial=key,first=a,second=b,error=error,passed=error<=.01))
    # Native diagonal continuity: only controls and original fractional pulses are compared.
    native_comparisons=[]
    for trial in trials:
        z=contexts[trial['context']]
        if z['chemical_donor']!=z['physical_donor']:continue
        previous=next(t for t in old['trials'] if t['context']==z['chemical_donor'] and t['target']==trial['target'] and t['dose']==trial['dose'])
        init=read_npz(previous['source'])['initial'];same=np.array_equal(init,states[trial['key']])
        assert trial['native_parent_preparation_exact']==same
        if not same:
            assert 'native_parent_trial' not in trial and trial['dose']=='amount'
            continue
        assert trial['native_parent_trial']==previous['key']
        for label in labels:
            if label in ('dop853','radau'):
                folder=parent/'references'/previous['key'];record=read_json(folder/'result.json');zz=read_npz(folder/'paths.npz');oldpath=zz[label]
            else:
                folder=parent/'production'/label;record=read_json(folder/'result.json');zz=read_npz(folder/'paths.npz');index=list(zz['trial_keys']).index(previous['key']);oldpath=zz['trajectory'][:,index]
            assert record['paths_sha256']==digest(folder/'paths.npz') and np.array_equal(zz['times'],times)
            error=log_error(paths[label,trial['key']],oldpath)
            native_comparisons.append(dict(trial=trial['key'],previous_trial=previous['key'],label=label,raw_log_error=error,within_backend_tolerance=bool(error<=1e-8)))
    assert len(native_comparisons)==135
    native_cases={r['trial'] for r in native_comparisons};assert len(native_cases)==27
    scientific=[];baseline=[];archive={'times':times};outcomes=[]
    for seed in p['histories']:
        for cell in p['selections'][str(seed)]['selected_ids']:
            for P in ('S','R'):
                tk={C:cases[seed,cell,C,P,'fractional'] for C in ('S','R')};tr={C:next(t for t in trials if t['key']==tk[C]) for C in ('S','R')};ctrl={C:controls[tr[C]['context']] for C in ('S','R')};z=contexts[tr['S']['context']];i=z['ids'].tolist().index(cell)
                row=dict(seed=seed,cell=cell,physical=P)
                for label in labels:
                    s,r=paths[label,ctrl['S']],paths[label,ctrl['R']];target=np.log(r[:,:,i]/s[:,:,i]);full=np.log(r/s);v=z['volumes'];network=np.sqrt(np.sum(full**2*v,axis=(1,2))/(2*v.sum()))
                    row[label]=dict(final_activator_ratio=float(r[-1,0,i]/s[-1,0,i]),final_inhibitor_ratio=float(r[-1,1,i]/s[-1,1,i]),late_target_log_rms=rms(target[late],times[late]),late_network_log_rms=float(np.sqrt(trapezoid(network[late]**2,times[late])/24.)),late_network_log_rms_max=float(network[late].max()),sham_whole_contrast_min=float(np.std(np.log(s[:,0]),axis=1).min()),redistributed_whole_contrast_min=float(np.std(np.log(r[:,0]),axis=1).min()),sham_final_log_drift_from_start=float(np.sqrt(np.sum(np.log(s[-1]/s[0])**2*v)/(2*v.sum()))),redistributed_final_log_drift_from_start=float(np.sqrt(np.sum(np.log(r[-1]/r[0])**2*v)/(2*v.sum()))))
                baseline.append(row)
            for dose in ('fractional','amount'):
                row=dict(seed=seed,cell=cell,dose=dose,chemical={},physical={},interaction={})
                for label in labels:
                    row['chemical'][label]={P:dict(rms=norms[label,seed,cell,dose,'chemical',P],detected=norms[label,seed,cell,dose,'chemical',P]>.01) for P in ('S','R')}
                    row['physical'][label]={C:norms[label,seed,cell,dose,'physical',C] for C in ('S','R')}
                    row['interaction'][label]=norms[label,seed,cell,dose,'interaction',None]
                row['native_current_dose_fine_rms']=native['fine',seed,cell,dose]['norm']
                row['fine_pulse_metrics']={C+P:metrics['fine',cases[seed,cell,C,P,dose]] for C in ('S','R') for P in ('S','R')}
                row['early_descriptive']={}
                for kind,held in [('chemical','S'),('chemical','R'),('physical','S'),('physical','R'),('interaction',None)]:
                    wave=contrasts['fine',seed,cell,dose,kind,held];energy=float(trapezoid(np.mean(wave**2,axis=1),times));early_energy=float(trapezoid(np.mean(wave[early]**2,axis=1),times[early]))
                    row['early_descriptive'][kind+'_'+str(held)]=dict(first6_rms=rms(wave[early],times[early]),first6_energy_fraction=early_energy/energy if energy else None)
                    archive[f'seed{seed}_cell{cell}_{dose}_{kind}_{held}']=wave
                for C in ('S','R'):
                    for P in ('S','R'):archive[f'seed{seed}_cell{cell}_{dose}_response_{C}{P}']=waves['fine',cases[seed,cell,C,P,dose]]
                scientific.append(row)
        row=dict(seed=seed,recipients=2,by_physics={})
        for P in ('S','R'):
            row['by_physics'][P]={dose:sum(norms['fine',seed,cell,dose,'chemical',P]>.01 for cell in p['selections'][str(seed)]['selected_ids']) for dose in ('fractional','amount')}
            row['by_physics'][P]['both_doses']=sum(all(norms['fine',seed,cell,dose,'chemical',P]>.01 for dose in ('fractional','amount')) for cell in p['selections'][str(seed)]['selected_ids'])
        row['both_doses_both_physics']=sum(all(norms['fine',seed,cell,dose,'chemical',P]>.01 for dose in ('fractional','amount') for P in ('S','R')) for cell in p['selections'][str(seed)]['selected_ids']);outcomes.append(row)
    duplicates=[v for v in groups.values() if len(v)>1];assert len(duplicates)==12
    duplicate_error=max(log_error(paths[label,keys[0]],paths[label,key]) for keys in duplicates for key in keys[1:] for label in labels)
    all_pass=all(x<=1e-8 for x in references.values()) and all(r['passed'] for r in accuracy+orders+rawpairs+comparisons) and all(r['within_backend_tolerance'] for r in native_comparisons if r['label'] in ('coarse','fine'))
    compare_record(all_pass,summary['passed'],'overall_pass');compare_record(all_pass,status['passed'],'status_pass')
    # Match every saved numerical/scientific report to the independently derived arrays.
    condition_to_held={'sham_physics':'S','redistributed_physics':'R','sham_chemistry':'S','redistributed_chemistry':'R','difference_of_differences':None}
    def report_key(row):return (row['seed'],row['cell'],row['dose'],row['kind'],condition_to_held[row['condition']])
    for label in labels:
        recorded=summary['response_summaries'][label]
        for trial in trials:
            if trial['target'] is None:continue
            for k,value in metrics[label,trial['key']].items():compare_record(value,recorded['metrics'][trial['key']][k],label+'/'+trial['key']+'/'+k)
        assert len(recorded['effects'])==60
        for row in recorded['effects']:
            key=(label,*report_key(row));value=norms[key]
            compare_record(value,row['normalized_response_shift_rms'],label+'/'+row['key'])
            if row['kind']=='chemical':compare_record(value>.01,row['response_shift_detected'],label+'/'+row['key']+'/decision')
            else:assert 'response_shift_detected' not in row
        for row in recorded['native_diagonal']:
            compare_record(native[label,row['seed'],row['cell'],row['dose']]['norm'],row['normalized_response_shift_rms'],label+'/native_diagonal')
    for comparison in comparisons:
        a,b=comparison['first'],comparison['second'];saved=next(r for r in summary['response_comparisons'] if r['first']==a and r['second']==b)
        compare_record(comparison['passed'],saved['passed'],a+'/'+b+'/pass')
        for row in comparison['pulse_rows']:
            previous=next(r for r in saved['pulse_rows'] if r['trial']==row['trial'])
            for k in ('normalized_response_error','relative_auc_error','passed'):compare_record(row[k],previous[k],a+'/'+b+'/'+row['trial']+'/'+k)
            for k,value in row['recovery_time_errors'].items():compare_record(value,previous['recovery_time_errors'][k],a+'/'+b+'/'+row['trial']+'/'+k)
        for row in comparison['contrast_rows']:
            held=row['held'];condition=('sham_physics' if held=='S' else 'redistributed_physics') if row['kind']=='chemical' else ('sham_chemistry' if held=='S' else 'redistributed_chemistry') if row['kind']=='physical' else 'difference_of_differences'
            key=f'{row["seed"]}_{row["cell"]}_{row["dose"]}_{row["kind"]}_{condition}'
            previous=next(r for r in saved['factorial_contrasts'] if r['key']==key)
            for mine,theirs in [('maximum_waveform_error','direct_waveform_error'),('rms_error','rms_error'),('same_primary_decision','same_scientific_decision'),('passed','passed')]:compare_record(row[mine],previous[theirs],a+'/'+b+'/'+key+'/'+mine)
    for row in accuracy:
        saved=next(r for r in summary['production_accuracy'] if r['trial']==row['trial'] and r['level']==row['level'])
        for mine,theirs in [('raw_reference_error','raw_reference_log_error'),('cpu_torch_error','cpu_torch_log_error'),('passed','passed')]:compare_record(row[mine],saved[theirs],row['trial']+'/'+row['level']+'/'+mine)
    for row in rawpairs:
        saved=next(r for r in summary['raw_timestep_comparisons'] if r['trial']==row['trial'] and r['first']==row['first'] and r['second']==row['second'])
        compare_record(row['error'],saved['raw_log_error'],row['trial']+'/rawpair');compare_record(row['passed'],saved['passed'],row['trial']+'/rawpairpass')
    for row in orders:
        saved=next(r for r in summary['orders'] if r['trial']==row['trial'])
        for mine,theirs in [('reference_disagreement','reference_noise'),('order_resolved','order_resolved'),('passed','passed')]:compare_record(row[mine],saved[theirs],row['trial']+'/order/'+mine)
        if row['orders'] is None:assert saved['orders'] is None
        else:
            for x,y in zip(row['orders'],saved['orders']):compare_record(x,y,row['trial']+'/orders')
    for row in native_comparisons:
        if row['label'] not in ('coarse','fine'):continue
        saved=next(r for r in summary['native_parent_comparisons'] if r['trial']==row['trial'] and r['level']==row['label'])
        compare_record(row['raw_log_error'],saved['raw_log_error'],row['trial']+'/native_parent');compare_record(row['within_backend_tolerance'],saved['passed'],row['trial']+'/native_parent_pass')
    for row in transplants:
        saved=p['assignments'][row['context']]
        mapping={'original_per_cell_amounts':'original_per_cell_species_amounts','reassigned_per_cell_amounts':'assigned_per_cell_species_amounts','per_cell_amount_change':'added_per_cell_species_amounts','original_total_amounts':'original_global_species_amounts','reassigned_total_amounts':'assigned_global_species_amounts','total_amount_change':'added_global_species_amounts','volume_ratio':'volume_ratios'}
        for mine,theirs in mapping.items():
            aa,bb=np.asarray(row[mine]),np.asarray(saved[theirs]);assert aa.shape==bb.shape
            for x,y in zip(aa.ravel(),bb.ravel()):compare_record(x,y,row['context']+'/assignment/'+mine)
    for row in baseline:
        saved=next(r for r in summary['baseline_descriptions']['paired_chemical_states'] if r['seed']==row['seed'] and r['cell']==row['cell'] and r['physics']==('sham' if row['physical']=='S' else 'redistributed'))
        values=row['fine']
        compare_record(values['final_activator_ratio'],saved['final_species_ratios'][0],'late_state/a_ratio')
        compare_record(values['final_inhibitor_ratio'],saved['final_species_ratios'][1],'late_state/b_ratio')
        compare_record(values['late_target_log_rms'],saved['late_two_species_state_log_rms'],'late_state/target_log_rms')
        compare_record(values['late_network_log_rms'],saved['late_network_two_species_state_log_rms'],'late_state/network_log_rms')
    control_drifts=[]
    for context,control in controls.items():
        path=paths['fine',control];z=contexts[context];v=z['volumes'];diff=np.log(path/path[0]);curve=np.sqrt(np.sum(np.mean(diff**2,axis=1)*v[None,:],axis=1)/v.sum());sd=np.std(np.log(path[:,0]),axis=1)
        row=dict(context=context,final_network_two_species_log_drift_rms=float(curve[-1]),maximum_network_two_species_log_drift_rms=float(curve.max()),late_network_two_species_log_drift_rms=float(np.sqrt(trapezoid(curve[late]**2,times[late])/24.)),final_log_activator_sd=float(sd[-1]),whole_min_log_activator_sd=float(sd.min()),late_min_log_activator_sd=float(sd[late].min()),final_global_species_amounts=(path[-1]@v).tolist())
        saved=next(r for r in summary['baseline_descriptions']['controls'] if r['context']==context)
        for k,value in row.items():
            if k in ('context','final_global_species_amounts'):continue
            compare_record(value,saved[k],context+'/control_drift/'+k)
        assert np.array_equal(np.array(saved['final_chemistry']),path[-1])
        for x,y in zip(row['final_global_species_amounts'],saved['final_global_species_amounts']):compare_record(x,y,context+'/final_amount')
        control_drifts.append(row)
    for row in outcomes:
        for P in ('S','R'):
            saved=next(r for r in summary['history_outcomes'] if r['seed']==row['seed'] and r['physics']==('sham_physics' if P=='S' else 'redistributed_physics'))
            compare_record(row['by_physics'][P]['both_doses'],saved['both_doses_above_screen'],'history_outcome')
            for dose in ('fractional','amount'):
                count=next(x['above_screen'] for x in saved['dose_counts'] if x['dose']==dose)
                compare_record(row['by_physics'][P][dose],count,'history_outcome/'+dose)
        saved=next(r for r in summary['reciprocal_history_outcomes'] if r['seed']==row['seed'])
        compare_record(row['both_doses_both_physics'],saved['both_doses_both_physical_contexts_above_screen'],'reciprocal_history_outcome')

    for file,expected in summary.get('evidence_sha256',{}).items():pin(file,expected)
    refqual=read_json(root/'reference-qualification.json')
    assert refqual['protocol_sha256']==ph
    for file,expected in refqual.get('evidence_sha256',{}).items():pin(file,expected)
    for name,expected in list(WATCHED.items()):
        if digest(name)!=expected:raise AssertionError('Changed consumed input during review: '+name)
    resolved=[v for r in orders if r['orders'] is not None for v in r['orders']]
    result=dict(protocol_sha256=ph,passed=bool(all_pass),scope='Independent saved-array reconstruction; no project analysis functions imported and no new simulations or project writes. The factorial preserves concentration vectors, not amounts; recipient baseline chemistry is not matched. Unsaved internal stage positivity remains runtime evidence, while saved positivity/checkpoint consistency are independently reconstructed.',
        counts=dict(contexts=21,logical_cases=69,unique_initial_value_problems=57,production_paths=207,cpu_paths=207,reference_trajectories=138,independent_histories=3,new_histories=0),source_checks=source_checks,transplant_amounts=transplants,dose_accounting=dose_details,maximum_dose_relative_error=max(dose_errors),stage_checks=stages,reference_errors=references,production_accuracy=accuracy,orders=orders,raw_timestep_comparisons=rawpairs,response_comparisons=comparisons,
        gate_counts=dict(reference_raw_pass=sum(v<=1e-8 for v in references.values()),production_accuracy_pass=sum(r['passed'] for r in accuracy),order_pass=sum(r['passed'] for r in orders),resolved_order_cases=sum(r['order_resolved'] for r in orders),raw_timestep_pass=sum(r['passed'] for r in rawpairs),comparison_sets_pass=sum(r['passed'] for r in comparisons),pulse_comparison_rows_pass=sum(r['passed'] for c in comparisons for r in c['pulse_rows']),contrast_comparison_rows_pass=sum(r['passed'] for c in comparisons for r in c['contrast_rows'])),
        maxima=dict(reference_raw=max(references.values()),reference_normalized=max(r['normalized_response_error'] for r in comparisons[0]['pulse_rows']),reference_contrast_waveform=max(r['maximum_waveform_error'] for r in comparisons[0]['contrast_rows']),reference_contrast_rms=max(r['rms_error'] for r in comparisons[0]['contrast_rows']),production_raw_reference=max(r['raw_reference_error'] for r in accuracy),cpu_torch_raw=max(r['cpu_torch_error'] for r in accuracy),normalized_response=max(r['normalized_response_error'] for c in comparisons[1:] for r in c['pulse_rows']),contrast_waveform=max(r['maximum_waveform_error'] for c in comparisons[1:] for r in c['contrast_rows']),contrast_rms=max(r['rms_error'] for c in comparisons[1:] for r in c['contrast_rows']),relative_auc=max(r['relative_auc_error'] for c in comparisons for r in c['pulse_rows']),recovery_time=max(v for c in comparisons for r in c['pulse_rows'] for v in r['recovery_time_errors'].values() if v is not None),order_min=min(resolved) if resolved else None,order_max=max(resolved) if resolved else None,vector_identity=max(identity_errors)),
        duplicate_groups=duplicates,duplicate_path_error=duplicate_error,native_continuity=native_comparisons,native_continuity_max_log_error=max(r['raw_log_error'] for r in native_comparisons),scientific_comparisons=scientific,baseline_descriptions=baseline,control_drifts=control_drifts,history_outcomes=outcomes,recorded_metrics_compared=len(RECORDED_DIFFERENCES),max_recorded_metric_difference=max(RECORDED_DIFFERENCES)[0],max_recorded_metric_difference_location=max(RECORDED_DIFFERENCES)[1],sampled_minimum_concentration=min(float(v.min()) for v in paths.values()),fine_all_paths_contrast_min=min(float(np.std(np.log(paths['fine',t['key']][:,0]),axis=1).min()) for t in trials),inputs_unchanged=True,inputs_sha256=WATCHED,script_sha256=digest(__file__))
    np.savez_compressed(output.with_suffix('.npz'),**archive);result['waveform_sha256']=digest(output.with_suffix('.npz'))
    output.with_suffix('.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ('passed','counts','gate_counts','maxima','history_outcomes','native_continuity_max_log_error','inputs_unchanged')},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path('/home/sji/Documents/Programming/Python/embryogenesis-simulation/output/studies/phase-carry-shared-context-response'))
    parser.add_argument('--output',type=Path,default=Path('/tmp/shared_context_independent_review'))
    args=parser.parse_args();main(args.root,args.output)
