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


def precheck(root,output):
    root=Path(root).resolve();output=Path(output)
    p=read_json(root/'protocol.json');ph=digest(root/'protocol.json')
    for field in ('source_sha256','input_sha256'):
        for file,expected in p[field].items():pin(file,expected)
    status=read_json(root/'status.json')
    assert status['state'] in ('prepared','running')
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
    for row in transplants:
        saved=p['assignments'][row['context']]
        for mine,theirs in [('original_per_cell_amounts','original_per_cell_species_amounts'),('reassigned_per_cell_amounts','assigned_per_cell_species_amounts'),('per_cell_amount_change','added_per_cell_species_amounts'),('original_total_amounts','original_global_species_amounts'),('reassigned_total_amounts','assigned_global_species_amounts'),('total_amount_change','added_global_species_amounts'),('volume_ratio','volume_ratios')]:
            assert np.array_equal(np.asarray(row[mine]),np.asarray(saved[theirs])),(row['context'],mine)
    # Coordinator status may legitimately advance while preparation hashes stay fixed.
    WATCHED.pop(str((root/'status.json').resolve()))
    for name,expected in list(WATCHED.items()):assert digest(name)==expected,name
    result=dict(passed=True,protocol_sha256=ph,contexts=len(contexts),logical_cases=len(trials),unique_initial_value_problems=len(groups),controls=len(controls),pulses=len(cases),maximum_dose_relative_error=max(dose_errors),source_checks=source_checks,transplants=transplants,doses=dose_details,duplicate_groups=[x for x in groups.values() if len(x)>1],inputs_unchanged=True,inputs_sha256=WATCHED,script_sha256=digest(__file__))
    Path(output).write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ('passed','protocol_sha256','contexts','logical_cases','unique_initial_value_problems','controls','pulses','maximum_dose_relative_error','inputs_unchanged')},indent=2))
    print('Immutable files',len(WATCHED))

if __name__=='__main__':
    precheck('/home/sji/Documents/Programming/Python/embryogenesis-simulation/output/studies/phase-carry-shared-context-response','/tmp/shared_context_prelaunch_audit.json')
