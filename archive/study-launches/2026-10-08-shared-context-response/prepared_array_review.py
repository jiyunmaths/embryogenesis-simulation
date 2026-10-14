"""Independent literal validation of every saved reciprocal-context preparation."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib, json, math, sys
sys.dont_write_bytecode=True
import numpy as np
ROOT=Path('/home/sji/Documents/Programming/Python/embryogenesis-simulation')
OLD=ROOT/'output/studies/phase-carry-frozen-response'
NEW=ROOT/'output/studies/phase-carry-shared-context-response'
REPORT=Path('/tmp/shared-context-prepared-array-review.json')
read=lambda p:json.loads(Path(p).read_text())
pins={};stats={}

def metadata(p):
    s=Path(p).stat();return dict(device=s.st_dev,inode=s.st_ino,bytes=s.st_size,mtime_ns=s.st_mtime_ns,mode=s.st_mode)

def sha(p,expected=None):
    p=Path(p);before=metadata(p)
    with p.open('rb') as f:value=hashlib.file_digest(f,'sha256').hexdigest()
    assert metadata(p)==before and (expected is None or expected==value),str(p)
    assert str(p) not in pins or pins[str(p)]==value
    pins[str(p)]=value;stats[str(p)]=before;return value

def values(p):
    with np.load(p,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}

def signature(*arrays):
    h=hashlib.sha256()
    for a in arrays:
        a=np.ascontiguousarray(a);h.update(a.dtype.str.encode());h.update(str(a.shape).encode());h.update(a.tobytes())
    return h.hexdigest()

report=dict(passed=False,started_utc=datetime.now(timezone.utc).isoformat(),
    scope='Read-only all-prepared-array/pin audit using literal NumPy construction independent of preparation/analysis helpers. No integration or permanent restoration.',
    integrations_performed=0,phase_arrays_loaded_or_copied=False,files_written_under_outputs=0)
try:
    p,old=read(NEW/'protocol.json'),read(OLD/'protocol.json')
    ph=sha(NEW/'protocol.json');oh=sha(OLD/'protocol.json')
    assert p['parent_protocol_sha256']==oh
    assert p['source_sha256'][str(ROOT/'embryo/phase_carry_shared_context_response.py')]=='ea63ff96ff9d54dffa758cdedeeb0e701b1647d8a08246ce8a7df32b173c820d'
    for key,expected in p['source_sha256'].items():sha(key,expected)
    for key,expected in p['input_sha256'].items():sha(key,expected)
    for key,expected in old['source_sha256'].items():assert p['source_sha256'][key]==expected
    for key in ('levels','beta','da','db','start','duration','interval','checkpoint_interval','late_window','effect_min','reference_settings','runtime','device','criteria','selections'):
        assert p[key]==old[key],key
    assert p['histories']==[7,8,9] and p['new_histories']==0 and p['independent_histories']==3 and p['parent_paths_reused']==0
    assert (p['primary_paths'],p['diagnostic_paths'],p['total_production_paths'],p['cpu_verification_paths'],p['logical_reference_cases'],p['reference_trajectories'])==(138,69,207,207,69,138)
    saved={c['key']:values(c['source']) for c in old['contexts']}
    contexts={c['key']:c for c in p['contexts']};loaded={};assignments=[]
    assert len(contexts)==21 and len(p['trials'])==69 and len(p['factorials'])==6
    for c in p['contexts']:
        chemical=saved[c['chemical_context']];physical=saved[c['physical_context']]
        x=values(c['source']);assert set(x)=={'initial','delta','volumes','ids'}
        for k in ('initial','ids'):assert np.array_equal(x[k],chemical[k])
        for k in ('delta','volumes','ids'):assert np.array_equal(x[k],physical[k])
        assert x['initial'].dtype==x['delta'].dtype==x['volumes'].dtype==np.float64
        assert x['initial'].shape==(2,16) and x['delta'].shape==(16,16) and x['volumes'].shape==x['ids'].shape==(16,)
        assert np.all(x['initial']>0) and np.all(x['volumes']>0)
        a=p['assignments'][c['key']];original=x['initial']*chemical['volumes'][None,:];assigned=x['initial']*x['volumes'][None,:]
        expected=dict(original_per_cell_species_amounts=original,assigned_per_cell_species_amounts=assigned,
            added_per_cell_species_amounts=assigned-original,original_global_species_amounts=original.sum(axis=1),
            assigned_global_species_amounts=assigned.sum(axis=1),added_global_species_amounts=(assigned-original).sum(axis=1),
            original_volumes=chemical['volumes'],assigned_volumes=x['volumes'],volume_ratios=x['volumes']/chemical['volumes'])
        assert a['concentrations_exact']
        for k,v in expected.items():assert np.array_equal(np.asarray(a[k]),v),(c['key'],k)
        assignments.append(dict(context=c['key'],exact_concentration_operator_volume_ids_and_accounting=True))
        loaded[c['key']]=x
    factorial_checks=[]
    for f in p['factorials']:
        seed,target=f['seed'],f['cell'];assert target in p['selections'][str(seed)]['selected_ids']
        assert set(f['contexts'])=={'S_S','S_R','R_S','R_R'}
        for factors,key in f['contexts'].items():
            C,G=factors.split('_');c=contexts[key]
            cn='sham' if C=='S' else f'redistributed-cell-{target}'
            gn='sham' if G=='S' else f'redistributed-cell-{target}'
            assert c['seed']==seed and c['chemical_context']==f'seed-{seed}_{cn}' and c['physical_context']==f'seed-{seed}_{gn}'
        factorial_checks.append(dict(seed=seed,target=target,all_four_factor_mappings_exact=True))
    unique={};trial_checks=[];errors=[]
    for t in p['trials']:
        c=contexts[t['context']];x=loaded[t['context']];initial=x['initial'];expected=initial.copy()
        record=dict(factor=1.,initial_log_amplitude=0.,requested_removed_amount=0.,removed_activator_amount=0.,
            initial_target_activator_amount=None,initial_target_volume=None,pulse_amount_relative_error=0.)
        if t['target'] is not None:
            i=x['ids'].tolist().index(t['target']);before=float(initial[0,i]*x['volumes'][i])
            aa=[saved[f'seed-{t["seed"]}_{name}']['initial'][0,i] for name in ('sham',f'redistributed-cell-{t["target"]}')]
            factor=.9 if t['dose']=='fractional' else float(1-.1*float(min(aa))/float(initial[0,i]))
            expected[0,i]*=factor
            requested=.1*before if t['dose']=='fractional' else float(.1*float(min(aa))*x['volumes'][i])
            measured=float((initial[0,i]-expected[0,i])*x['volumes'][i]);error=abs(measured/requested-1.)
            assert error<=2e-13;errors.append(error)
            record.update(factor=factor,initial_log_amplitude=float(abs(np.log(factor))),requested_removed_amount=requested,
                removed_activator_amount=measured,initial_target_activator_amount=before,initial_target_volume=float(x['volumes'][i]),pulse_amount_relative_error=error)
        for k,v in record.items():assert t[k]==v,(t['key'],k,t[k],v)
        actual=values(t['source']);assert set(actual)=={'initial'} and np.array_equal(actual['initial'],expected)
        sig=signature(x['delta'],x['volumes'],x['ids'],expected);unique.setdefault(sig,[]).append(t['key'])
        trial_checks.append(dict(trial=t['key'],literal_saved_array_and_dose_metadata_exact=True))
    assert len(unique)==p['unique_initial_value_problems']==57
    assert len([g for g in unique.values() if len(g)>1])==12
    mapping=[t['key'] for c in p['contexts'] for t in p['trials'] if t['context']==c['key']]
    assert mapping==[t['key'] for t in p['trials']]
    width=max(sum(t['context']==c['key'] for t in p['trials']) for c in p['contexts'])
    assert width==5 and 21*width-len(mapping)==36
    step_counts={}
    for c in p['contexts']:
        x=loaded[c['key']];loss=max(1+p['da']*float(-np.diag(x['delta']).min()),p['beta']+p['db']*float(-np.diag(x['delta']).min()))
        step_counts[c['key']]={l['level']:max(1,math.ceil(l['dt']*loss/.2)) for l in p['levels']}
        assert set(step_counts[c['key']].values())=={1}
    before=dict(pins)
    for file,value in list(before.items()):sha(file,value)
    for file,value in stats.items():assert metadata(file)==value
    report.update(passed=True,protocol_sha256=ph,parent_protocol_sha256=oh,
        source_hashes_verified=len(p['source_sha256']),input_hashes_verified=len(p['input_sha256']),
        contexts_verified=21,trials_verified=69,factorials_verified=6,unique_initial_value_problems=57,
        maximum_amount_relative_error=max(errors),assignment_checks=assignments,trial_checks=trial_checks,
        factorial_checks=factorial_checks,context_grouped_mapping_exact=True,batch_blocks=21,batch_matrix_shape=[336,336],
        batch_state_shape=[2,336,5],logical_paths=69,padding_paths=36,positive_substep_counts=step_counts,
        immutable_before_sha256=before,immutable_after_sha256=dict(pins),immutable_stat_guards=stats,
        immutable_before_after_equal=before==pins,all_mismatches=[],script_path=str(Path(__file__)),script_sha256=sha(__file__))
except BaseException as exc:
    import traceback
    report.update(error=str(exc),traceback=traceback.format_exc(),all_mismatches=[str(exc)]);raise
finally:
    report['finished_utc']=datetime.now(timezone.utc).isoformat()
    REPORT.write_text(json.dumps(report,indent=2)+'\n')
    print('Prepared arrays review passed',report['passed'],REPORT)
