"""Read-only frozen-parent pins, shared-context dose/equality, and readiness audit."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib, json, math, os, sys, time
sys.dont_write_bytecode = True
import numpy as np
ROOT = Path('/home/sji/Documents/Programming/Python/embryogenesis-simulation')
PARENT = ROOT/'output/studies/phase-carry-frozen-response'
REPORT = Path('/tmp/shared-context-input-review.json')
read = lambda p: json.loads(Path(p).read_text())
evidence, stats = {}, {}

def stat(p):
    s = Path(p).stat()
    return dict(device=s.st_dev, inode=s.st_ino, bytes=s.st_size, mode=s.st_mode,
                mtime_ns=s.st_mtime_ns, nlink=s.st_nlink)

def sha(p, expected=None):
    p = Path(p); before = stat(p)
    with p.open('rb') as f: value = hashlib.file_digest(f, 'sha256').hexdigest()
    assert before == stat(p), str(p)
    assert expected is None or expected == value, (str(p), expected, value)
    assert str(p) not in evidence or evidence[str(p)] == value, str(p)
    evidence[str(p)], stats[str(p)] = value, before
    return value

def fingerprint(*arrays):
    h = hashlib.sha256()
    for a in arrays:
        a = np.ascontiguousarray(a)
        h.update(a.dtype.str.encode()); h.update(str(a.shape).encode()); h.update(a.tobytes())
    return h.hexdigest()

started = time.monotonic()
report = dict(passed=False, started_utc=datetime.now(timezone.utc).isoformat(),
    scope='Read-only preparation audit. New backgrounds retain full concentrations and inherit measured graph/volumes independently; no simulation or permanent restore.',
    files_written_under_outputs=0, phase_arrays_loaded_or_copied=False,
    permanent_decompression=False)
try:
    assert (ROOT/'outputs').is_symlink()
    assert os.readlink(ROOT/'outputs') == '/mnt/data/Documents/Python/embryogenesis-simulation/outputs'
    sha(ROOT/'archive/output-relocation/2026-10-07/manifest.json',
        'edcf6b433ed7d3db70adb89c8ca701ec49b0d73a5b981380c61a966bf89a3c02')
    p = read(PARENT/'protocol.json'); ph = sha(PARENT/'protocol.json')
    status, summary = read(PARENT/'status.json'), read(PARENT/'summary.json')
    compactfile = ROOT/'docs/phase_carry_frozen_response_assessment.json'
    compact = read(compactfile); sha(compactfile)
    assert status['state'] == 'completed' and status['passed'] and summary['passed']
    assert compact['full_study_accepted'] and compact['protocol_sha256'] == ph == summary['protocol_sha256']
    assert status['total_production_paths'] == summary['total_production_paths'] == 99
    assert compact['independent_histories'] == 3 and compact['new_histories'] == 0
    sha(PARENT/'status.json', compact['original_status_sha256'])
    sha(PARENT/'summary.json', compact['original_summary_sha256'])
    rr = ROOT/compact['review_archive']; rp = rr/'review-provenance.json'
    review = read(rp); sha(rp)
    assert review['passed'] and review['source_and_input_pins_unchanged'] and review['protocol_sha256'] == ph
    for k in ('original_summary_sha256', 'original_status_sha256'):
        assert review[k] == compact[k]
    exclusions = []
    for name, expected in review['files_sha256'].items():
        if name in ('README.md', 'docs/plan.md'):
            exclusions.append(name); continue
        sha(ROOT/name, expected)
    proof = read(compact['independent_review'])
    sha(compact['independent_review'], compact['independent_review_sha256'])
    assert proof['passed'] and proof['inputs_unchanged'] and proof['protocol_sha256'] == ph
    for name, expected in p['source_sha256'].items(): sha(name, expected)
    print('Verified accepted review and56 source pins', flush=True)
    for name, expected in p['input_sha256'].items(): sha(name, expected)
    print('Verified all115 parent input pins without loading phase fields', flush=True)
    old = {}; minimal = {}
    for c in p['contexts']:
        assert c['source'] in proof['inputs_sha256'] and c['source'] in p['input_sha256']
        sha(c['source'], proof['inputs_sha256'][c['source']]); minimal[c['source']] = evidence[c['source']]
        with np.load(c['source'], allow_pickle=False) as z:
            assert set(z.files) == {'delta', 'volumes', 'initial', 'ids'}
            x = {name: z[name].copy() for name in z.files}
        d, v, u, ids = [x[name] for name in ('delta','volumes','initial','ids')]
        assert d.dtype == v.dtype == u.dtype == np.float64
        assert d.shape == (16,16) and v.shape == (16,) and u.shape == (2,16) and ids.shape == (16,)
        assert all(np.isfinite(a).all() for a in (d,v,u)) and np.all(v > 0) and np.all(u > 0)
        assert len(set(ids.tolist())) == 16
        assert np.max(np.abs(d@np.ones(16))) < 1e-11
        assert np.max(np.abs(v@d)) < 1e-12
        k = -v[:,None]*d
        assert np.max(np.abs(k-k.T)) < 1e-12
        assert np.min(np.linalg.eigvalsh((k+k.T)/2)) > -1e-11
        control = next(t for t in p['trials'] if t['context'] == c['key'] and t['target'] is None)
        sha(control['source'], proof['inputs_sha256'][control['source']])
        with np.load(control['source'],allow_pickle=False) as z:
            assert np.array_equal(z['initial'],u)
        minimal[control['source']] = evidence[control['source']]
        old[c['key']] = x
    native = {}
    for t in p['trials']:
        sha(t['source'],proof['inputs_sha256'][t['source']])
        with np.load(t['source'],allow_pickle=False) as z: initial = z['initial'].copy()
        x = old[t['context']]
        native[t['key']] = dict(initial=initial, signature=fingerprint(x['delta'],x['volumes'],x['ids'],initial))
    contexts = {}; trials = []; doses = []; continuity = []; one_ulp = []
    for seed in p['histories']:
        s = old[f'seed-{seed}_sham']
        for target in p['selections'][str(seed)]['selected_ids']:
            r = old[f'seed-{seed}_redistributed-cell-{target}']
            assert np.array_equal(s['ids'],r['ids'])
            i = s['ids'].tolist().index(target)
            for G,g in [('sham',s),(f'redistributed-cell-{target}',r)]:
                for C,c in [('sham',s),(f'redistributed-cell-{target}',r)]:
                    key = f'seed-{seed}_chem-{C}_physical-{G}'
                    x = dict(delta=g['delta'], volumes=g['volumes'], ids=g['ids'], initial=c['initial'])
                    if key not in contexts:
                        contexts[key] = x
                        trials.append(dict(key=key+'_control', context=key, seed=seed,
                            target=None, dose=None, factor=1., initial=c['initial'].copy()))
                    aS,aR = float(s['initial'][0,i]),float(r['initial'][0,i])
                    a,v = float(c['initial'][0,i]),float(g['volumes'][i])
                    requested = .1*min(v*aS,v*aR)
                    fA = 1-requested/(v*a)
                    fB = 1-.1*min(aS,aR)/a
                    if fA != fB:
                        one_ulp.append(dict(seed=seed,target=target,physical=G,chemistry=C,
                            amount_first_factor=fA,concentration_first_factor=fB,
                            target_state_equal=(a*fA==a*fB)))
                    for dose, factor in [('fractional',.9),('amount',fB)]:
                        changed = c['initial'].copy(); changed[0,i] *= factor
                        removal = float((a-changed[0,i])*v)
                        desired = .1*v*a if dose == 'fractional' else requested
                        relerr = abs(removal/desired-1.)
                        assert 0 < factor < 1 and relerr < 1e-10
                        trial = dict(key=key+f'_cell-{target}_{dose}', context=key, seed=seed,
                            target=target, dose=dose, factor=factor, initial=changed)
                        trials.append(trial)
                        if dose == 'amount':
                            fractional = c['initial'].copy(); fractional[0,i] *= .9
                            doses.append(dict(seed=seed,target=target,physical=G,chemistry=C,
                                target_concentration=a,target_volume=v,common_removed_amount=requested,
                                actual_removed_amount=removal,relative_amount_error=relerr,
                                amount_factor=factor,amount_first_factor=fA,
                                lower_concentration=a==min(aS,aR),
                                byte_identical_to_fractional=np.array_equal(changed,fractional)))
    assert len(contexts)==21 and len(trials)==69
    groups = {}
    for t in trials:
        x = contexts[t['context']]
        sig = fingerprint(x['delta'],x['volumes'],x['ids'],t['initial'])
        groups.setdefault(sig,[]).append(t['key'])
        matches = [key for key,n in native.items() if n['signature'] == sig]
        if matches:
            continuity.append(dict(key=t['key'],target=t['target'],dose=t['dose'],
                parent_exact_matches=matches,identical_full_initial_graph_volume_ids=True))
    duplicate_groups = [names for names in groups.values() if len(names)>1]
    assert len(groups)==57 and len(duplicate_groups)==12
    assert len(continuity)==27
    assert sum(t['target'] is None for t in continuity)==9
    assert sum(t['dose']=='fractional' for t in continuity)==12
    assert sum(t['dose']=='amount' for t in continuity)==6
    assert sum(d['lower_concentration'] for d in doses)==12
    assert all(d['byte_identical_to_fractional'] == d['lower_concentration'] for d in doses)
    for l in p['levels']:
        folder=PARENT/'production'/l['level']
        for filename in ('paths.npz','checkpoint.npz','result.json'):
            file=folder/filename; sha(file,proof['inputs_sha256'][str(file)])
        result=read(folder/'result.json')
        assert result['protocol_sha256']==ph and result['level']==l
        with np.load(folder/'paths.npz',allow_pickle=False) as z:
            assert len(z['times'])==401 and z['times'][0]==0. and z['times'][-1]==60.
            assert z['trajectory'].shape==(401,33,2,16)
            assert list(z['trial_keys'])==[t['key'] for t in p['trials']]
            for i,t in enumerate(p['trials']): assert np.array_equal(z['trajectory'][0,i],native[t['key']]['initial'])
    assert (p['beta'],p['da'],p['db'],p['start'],p['duration'],p['interval'])==(2.,.02,.55,510.,60.,.15)
    hardware=read('/tmp/shared-context-hardware-readiness.json')
    assert hardware['cuda_available'] and hardware['logical_gpu_count']==1
    assert hardware['device_name']==p['runtime']['expected_gpu_name']=='NVIDIA GeForce GTX 1080 Ti'
    assert hardware['torch']==p['runtime']['torch'] and hardware['torch_cuda']==p['runtime']['torch_cuda']
    assert hardware['compute_apps']=='' and hardware['processes']==[]
    print('Verified21 backgrounds,69 logical starts,57 exactIVPs,27 native equalitychecks',flush=True)
    before = dict(evidence)
    for name, expected in list(before.items()): sha(name,expected)
    for name, original in stats.items(): assert stat(name)==original
    accepted_pins = {str(PARENT/f):evidence[str(PARENT/f)] for f in ('protocol.json','status.json','summary.json')}
    accepted_pins.update({str(compactfile):evidence[str(compactfile)],str(rp):evidence[str(rp)]})
    accepted_pins.update({str(ROOT/n):evidence[str(ROOT/n)] for n in review['files_sha256'] if n not in exclusions})
    context_descriptions=[]; step_counts = {}
    for key,x in contexts.items():
        context_descriptions.append(dict(key=key,arrays_sha256=fingerprint(x['delta'],x['volumes'],x['ids'],x['initial']),
            total_activator_amount=float(x['initial'][0]@x['volumes']),
            total_inhibitor_amount=float(x['initial'][1]@x['volumes']),
            small_arrays_bytes=sum(a.nbytes for a in x.values())))
        loss = max(1+p['da']*float(-np.diag(x['delta']).min()),p['beta']+p['db']*float(-np.diag(x['delta']).min()))
        step_counts[key] = {l['level']:max(1,math.ceil(l['dt']*loss/.2)) for l in p['levels']}
        assert set(step_counts[key].values())=={1}
    report.update(passed=True,protocol_sha256=ph,source_hashes_verified=len(p['source_sha256']),
        parent_input_hashes_verified=len(p['input_sha256']),parent_review_artifacts_verified=len(review['files_sha256'])-2,
        minimal_accepted_review_pins=accepted_pins,minimal_context_control_pins=minimal,
        inherited_source_pins=p['source_sha256'],source_and_input_pins_unchanged=True,
        excluded_reporting_guards=exclusions,historical_document_guards_not_assay_inputs=True,
        parent_compact_schema='independent_review + independent_review_sha256; no full_assessment',
        parent_independent_evidence_schema='inputs_sha256; no evidence_sha256',
        parent_contexts_verified=9,shared_contexts=context_descriptions,shared_contexts_count=21,
        recipient_factorial_context_labels=24,logical_starts_per_level=69,
        unique_initial_value_problems_per_level=57,duplicate_groups=duplicate_groups,
        production_paths_all_new_at_three_levels=207,logical_reference_cases=69,reference_trajectories_two_methods=138,
        native_exact_continuity_cases=continuity,native_exact_matches_count=27,
        native_controls=9,native_fractional_pulses=12,native_lower_common_amount_pulses=6,
        native_trajectory_reuse_recommended=False,
        rationale='All69 starts should be rerun perlevel;27 exact native cases support numerical continuity. Do not reuse old upper common-amount curves because their dose was matched across different volumes.',
        dose_definition='fractional factor0.9; common factor1-0.1*min(a_sham,a_red)/a_chemical with separately requested0.1*min(V*a_sham,V*a_red); preserve all initial concentrations.',
        amount_design_checks=doses,operation_order_factor_differences=one_ulp,
        maximum_amount_relative_error=max(d['relative_amount_error'] for d in doses),
        total_small_context_arrays_bytes=sum(c['small_arrays_bytes'] for c in context_descriptions),
        inherited_positive_substep_counts=step_counts,observation_clock='elapsed0..60 at interval0.15, absolute t510..570',
        raw_three_level_trajectory_bytes=207*401*2*16*8,
        minimal_context_operator_state_only=True,output_path_within_writable_workspace=str(ROOT/'output/studies/phase-carry-shared-context-response'),
        original_outputs_alias=os.readlink(ROOT/'outputs'),hardware=hardware,
        hardware_command='env CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1 PYTHONDONTWRITEBYTECODE=1 /home/sji/miniforge3/envs/deeplearning/bin/python - (host require_escalated read-only)',
        immutable_before_sha256=before,immutable_after_sha256=dict(evidence),immutable_stat_guards=stats,
        immutable_before_after_equal=before==evidence,all_mismatches=[],script_path=str(Path(__file__)),script_sha256=sha(Path(__file__)))
except BaseException as exc:
    import traceback
    report.update(error=str(exc),traceback=traceback.format_exc(),all_mismatches=[str(exc)])
    raise
finally:
    report['finished_utc']=datetime.now(timezone.utc).isoformat()
    report['elapsed_seconds']=time.monotonic()-started
    REPORT.write_text(json.dumps(report,indent=2)+'\n')
    print('Shared context input review passed',report['passed'],REPORT,flush=True)
