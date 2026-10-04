"""Matched initiation controls across three histories using phase-update carry.

Two directional-tension contrasts and two timesteps; reuse the completed
history-9 zero-contrast pair read-only. No new histories or pattern starts.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
import fcntl
import multiprocessing
from pathlib import Path
import shutil
import time

import numpy as np
import scipy
import torch

from .attribute_development import AttributeSimulation
from .feedback_long import digest, restore_checkpoint, save_checkpoint
from .geometry_precision import checked_history, implementation_gate
from .gpu_backend import library as accepted_library
from .gpu_precision_control import PrecisionSimulation, library
from .gpu_response_runner import require_validation
from .mechanics_float64_reference import run as reference_run, CRITERIA as REFERENCE_CRITERIA
from .neighbor_context import read
from .parameter_robustness import verify
from .parameter_robustness_moving import initialize, prefix, retain, endpoint_assay
from .phase_carry_convergence import completed_evidence, reference_evidence
from .polarity_robustness import assert_same_physical_start, first_crossing, refinement_comparison
from .polarity_robustness_refinement import check_halving
from .resolution import _steps, write_json


def job_design(sources):
    jobs = []
    # The previously initiating history is checked first, before unselected repeats.
    for seed in (9, 7, 8):
        for chi in (0., .35):
            for level, dt in (('coarse', .00375), ('fine', .001875)):
                jobs.append(dict(key=f'seed-{seed}_{level}_chi-{chi:g}', seed=seed, level=level,
                    arm='phase_carry', dt=dt, family='uniform', point=dict(key=f'polarity-{chi:g}', ratio=27.5, chi=chi),
                    source=sources[seed][level], chemical_file=sources[seed]['chemical']))
    return jobs


def parent_prefix(parent, p, job):
    folder = Path(parent)/job['key']; gate = read(folder/'prefix/comparison.json')
    if gate['protocol_sha256'] != digest(Path(parent)/'protocol.json') or gate['job'] != job or not gate['passed']:
        raise ValueError('Unaccepted historical context gate')
    for file, h in gate['evidence_sha256'].items():
        if digest(file) != h: raise ValueError('Changed historical context evidence')
    return [folder/'prefix/comparison.json', *[Path(f) for f in gate['evidence_sha256']]]


def prepare(root, qualification=Path('outputs/phase-carry-convergence'), parent=Path('outputs/polarity-robustness')):
    root, qualification, parent = [Path(x).resolve() for x in (root, qualification, parent)]
    if root.exists(): raise FileExistsError(root)
    q, old = read(qualification/'protocol.json'), read(parent/'protocol.json')
    verify(q); verify(old)
    qs = read(qualification/'summary.json')
    if (read(qualification/'status.json')['state'] != 'completed' or not qs['passed'] or
            qs['protocol_sha256'] != digest(qualification/'protocol.json')):
        raise ValueError('Complete the three-timestep qualification first')
    reference_evidence(qualification, q)
    files = [qualification/f for f in ('protocol.json', 'summary.json', 'status.json', 'reference-verification.json', 'float64-reference/result.json')]
    files += [Path(f) for f in qs['reference']['field_sha256']]
    files += [parent/f for f in ('protocol.json', 'summary.json', 'status.json')]
    for level, validation in old['validation_roots'].items():
        config, evidence = require_validation(validation)
        if asdict(config) != old['accepted_configs'][level]: raise ValueError('Changed accepted context configuration')
        files += evidence
    root.mkdir(parents=True); sources = {}
    for seed in (7, 8, 9):
        folder = root/f'seed-{seed}'; folder.mkdir()
        sources[seed] = {}
        for label, filename in (('coarse', 'coarse-source.npz'), ('fine', 'fine-source.npz'), ('chemical', 'initial_states.npz')):
            original, destination = parent/f'seed-{seed}'/filename, folder/filename
            shutil.copy2(original, destination); files += [original, destination]
            sources[seed][label] = str(destination)
        check_halving(AttributeSimulation.restore(sources[seed]['coarse']), AttributeSimulation.restore(sources[seed]['fine']))
    jobs = job_design(sources)
    # Completed zero-contrast history-9 paths are numerical reuse, not new trials.
    reused = {}
    candidate_reuse = {'coarse':dict(root=str(qualification), job=q['jobs'][0]),
                      'fine':dict(root=q['parent'], job=q['reused_references']['fine']['job'])}
    for level, item in candidate_reuse.items():
        op = read(Path(item['root'])/'protocol.json'); verify(op)
        completed_evidence(item['root'], op, item['job'])
        job = next(j for j in jobs if j['seed'] == 9 and j['level'] == level and j['point']['chi'] == 0.)
        a = initialize(job['source'], job['chemical_file'], 'uniform', job['point'], old['accepted_configs'][level])
        b = initialize(item['job']['source'], item['job']['chemical_file'], 'uniform', item['job']['point'], op['accepted_configs'][level])
        assert_same_physical_start(a, b)
        if asdict(a.config) != asdict(b.config): raise ValueError('Reused start has changed parameters')
        reused[job['key']] = item
        folder = Path(item['root'])/item['job']['key']
        ep = read(folder/'endpoint/protocol.json'); verify(ep); e = read(folder/'endpoint/assay/result.json')
        if (not e['numerical_pass'] or not e['all_trials_settled'] or
                e['paths_sha256'] != digest(folder/'endpoint/assay/paths.npz') or
                e['protocol_sha256'] != digest(folder/'endpoint/protocol.json') or
                ep['moving_protocol_sha256'] != digest(Path(item['root'])/'protocol.json') or
                ep['moving_history_sha256'] != digest(folder/'history.json')):
            raise ValueError('Changed reusable endpoint')
        files += [Path(item['root'])/'protocol.json']
        files += [folder/f for f in ('result.json', 'history.json', 'latest_state.npz', 'endpoint/protocol.json', 'endpoint/source.npz', 'endpoint/assay/result.json', 'endpoint/assay/paths.npz')]
    historical_gates = {}
    for job in jobs:
        if job['key'] in reused: continue
        match = next((j for j in old['jobs'] if j['seed'] == job['seed'] and j['level'] == job['level'] and
                      j['point'] == job['point'] and j['family'] == 'uniform'), None)
        if match:
            files += parent_prefix(parent, old, match)
            historical_gates[job['key']] = dict(root=str(parent), job=match)
    # Prepare a true positive-chi state for the independent component check.
    positive_job = next(j for j in jobs if j['seed'] == 9 and j['level'] == 'coarse' and j['point']['chi'] == .35)
    host = initialize(positive_job['source'], positive_job['chemical_file'], 'uniform', positive_job['point'], old['accepted_configs']['coarse'])
    reference_source = root/'positive-chi-reference-source.npz'; host.checkpoint(reference_source); files.append(reference_source)
    binary, original_binary = Path(library()._name), Path(accepted_library()._name)
    if (digest(binary) != q['device']['experimental_binary_sha256'] or
            digest(original_binary) != q['device']['accepted_binary_sha256'] or
            torch.cuda.get_device_name(0) != q['device']['name']): raise ValueError('Changed qualified GPU')
    files += [binary, original_binary]
    sources_sha = {**old['source_sha256'], **q['source_sha256'], str(Path(__file__).resolve()):digest(__file__)}
    p = dict(jobs=jobs, reused=reused, historical_gates=historical_gates,
        parent=str(parent), qualification=str(qualification), accepted_configs=old['accepted_configs'], device=q['device'],
        start=150., duration=240., pilot_duration=60., interval=.15, checkpoint_interval=3.,
        late_window=24., pilot_late_window=12., independent_histories=3, new_histories=0,
        histories=[7,8,9], contrasts=[0.,.35], dts=[.00375,.001875],
        new_moving_jobs=10, reused_trajectories=2, endpoint_workers=2,
        prefix_duration=.6, prefix_native_threads=4, prefix_criteria=old['prefix_criteria'],
        criteria=old['criteria'], refinement_criteria=old['refinement_criteria'],
        endpoint_horizons=old['endpoint_horizons'], endpoint_interval=old['endpoint_interval'], endpoint_criteria=old['endpoint_criteria'],
        reference_contexts=[dict(key='history-9-positive-chi', source=str(reference_source))],
        reference_dts=[.00375,.001875,.0009375], reference_horizon=.15, reference_criteria=REFERENCE_CRITERIA,
        reference_runtime=dict(numpy=np.__version__, scipy=scipy.__version__),
        estimated_new_moving_seconds=5*(qs['results'][-1]['result']['audit']['wall_seconds']+
            next(r['result']['audit']['wall_seconds'] for r in qs['results'] if r['level']=='fine')),
        design='Matched near-uniform chemical restarts on existing mature t=150 geometries from histories 7,8,9. At fixed beta=2, D_a=.02, D_b=.55, compare chi=0 and .35 while activity tension/adhesion remain .25/.35 and polarity dynamics continue. Each context is checked at dt=.00375 and .001875 to elapsed240. No new zygote histories, noise realizations, pattern starts or nutrient/oxygen model.',
        gates='Verify completed three-timestep carry qualification, full accepted backend evidence, and matched starts. Require a positive-chi float64 mechanics component check before moving jobs; per-context baseline native/GPU evidence, exact accepted/experimental baseline checks and residual restart checks. Each paired context must pass the original 60-unit refinement screen before its NEW long continuations. A failed pair stays at its pilot and does not prevent independent contexts from being tested.',
        interpretation='Primary scientific question: does directional polarity mechanics still reduce formation opportunity after the precision correction? A history supports suppression only if its zero-contrast control initiates sustained contrast and the .35 branch does not, with passing numerical checks. Both uniform means no successful initiating control, not evidence for polarity-specific suppression. Outcomes need not support the hypothesis to pass numerical quality. Three reused developmental histories; twelve nested numerical paths. Mature near-uniform restart is not initiation from fresh zygotes, maintenance, inheritance, autonomous identity, spatial convergence, or arbitrary GPU/backend acceptance.',
        source_sha256=sources_sha, input_sha256={**q['input_sha256'], **{str(f.resolve()):digest(f) for f in files}})
    write_json(root/'protocol.json', p)
    write_json(root/'status.json', dict(state='prepared', completed=0, total=10, reused_trajectories=2,
        independent_histories=3, new_histories=0, estimated_new_moving_seconds=p['estimated_new_moving_seconds']))
    return p


def context_gate(root, job, p):
    folder = Path(root)/job['key']; folder.mkdir(exist_ok=True); result_file = folder/'context-gate.json'
    ph = digest(Path(root)/'protocol.json')
    if result_file.exists():
        r = read(result_file)
        if not r['passed'] or r['protocol_sha256'] != ph or r['job'] != job: raise ValueError('Changed context gate')
        for f,h in r['evidence_sha256'].items():
            if digest(f) != h: raise ValueError('Changed context gate evidence')
        return r
    historical = p['historical_gates'].get(job['key']); files = []
    if historical:
        op = read(Path(historical['root'])/'protocol.json'); verify(op)
        old_job = historical['job']; files += parent_prefix(historical['root'], op, old_job)
        a = initialize(job['source'],job['chemical_file'],'uniform',job['point'],p['accepted_configs'][job['level']])
        b = initialize(old_job['source'],old_job['chemical_file'],'uniform',old_job['point'],op['accepted_configs'][old_job['level']])
        assert_same_physical_start(a,b)
        if asdict(a.config) != asdict(b.config): raise ValueError('Historical context parameters changed')
        baseline_gate = read(Path(historical['root'])/old_job['key']/'prefix/comparison.json')
    else:
        view = dict(p, accepted_config=p['accepted_configs'][job['level']], dt=job['dt'])
        baseline_gate = prefix(root,job,view)
        files += [folder/'prefix/comparison.json', *[Path(f) for f in baseline_gate['evidence_sha256']]]
    gate_folder = folder/'implementation'; gate_folder.mkdir(exist_ok=True)
    implementation = implementation_gate(gate_folder, job['source'], job['chemical_file'], job['point'], p['accepted_configs'][job['level']])
    files += [gate_folder/'gate-checkpoint.npz']
    # Check the actual contact32 carry arm, independently of the existing both-arm gate.
    sim = PrecisionSimulation(initialize(job['source'],job['chemical_file'],'uniform',job['point'],p['accepted_configs'][job['level']]), 'phase_carry')
    for _ in range(4): sim.step()
    cp = gate_folder/'phase-carry-restart.npz'; sim.checkpoint(cp); restored = PrecisionSimulation.restore(cp)
    for _ in range(2):
        sim.step(); restored.step()
        for key in ('phi','activator','inhibitor','polarity','geometry','phase_carry','rounding'):
            if not torch.equal(getattr(sim,key),getattr(restored,key)): raise ValueError('Context carry restart differs: '+key)
    del sim,restored; torch.cuda.empty_cache(); files.append(cp)
    result = dict(protocol_sha256=ph,job=job,passed=True,baseline_native_gpu=baseline_gate,
        historical_baseline_gate=historical is not None,implementation=implementation,actual_phase_carry_restart_exact_steps=2,
        evidence_sha256={str(f.resolve()):digest(f) for f in files})
    write_json(result_file,result); return result


def current_evidence(root,p,job):
    """Validate new results with a tolerance for equivalent floating clocks."""
    folder=Path(root)/job['key'];ph=digest(Path(root)/'protocol.json')
    r,h=read(folder/'result.json'),read(folder/'history.json')
    if (r['job']!=job or r['protocol_sha256']!=ph or not r['quality_pass'] or
            r['history_sha256']!=digest(folder/'history.json') or
            r['checkpoint_sha256']!=digest(folder/'latest_state.npz')):raise ValueError('Changed completed carry evidence')
    checked_history(h,p,job,p['duration'])
    host,audit,saved=restore_checkpoint(folder/'latest_state.npz',ph,job)
    if (saved!=h or audit!=r['audit'] or abs(host.time-p['start']-p['duration'])>1e-9 or
            host.step_number!=_steps(p['start']+p['duration'],job['dt']) or
            not np.array_equal([host.activator,host.inhibitor],h[-1]['chemistry'])):
        raise ValueError('Carry checkpoint/history mismatch')
    with np.load(folder/'latest_state.npz') as z:
        if (str(z['precision_arm'])!='phase_carry' or z['phase_carry'].dtype!=np.float64 or
                z['phase_carry'].shape!=host.phi.shape or not np.isfinite(z['phase_carry']).all()):
            raise ValueError('Invalid completed carry state')
    return r,h


def advance(root, job, p, horizon):
    """Stop/resume at declared horizons without publishing a pilot as a result."""
    if horizon not in (p['pilot_duration'],p['duration']): raise ValueError('Undeclared horizon')
    root=Path(root); reuse=p['reused'].get(job['key'])
    if reuse:
        op=read(Path(reuse['root'])/'protocol.json')
        _,h=completed_evidence(reuse['root'],op,reuse['job'])
        return checked_history([r for r in h if r['elapsed']<=horizon+1e-9],p,job,horizon)
    context_gate(root,job,p)
    folder=root/job['key']; cp=folder/'latest_state.npz'; ph=digest(root/'protocol.json')
    if (folder/'result.json').exists():
        _,h=current_evidence(root,p,job)
        return checked_history([r for r in h if r['elapsed']<=horizon+1e-9],p,job,horizon)
    if cp.exists():
        host,audit,history=restore_checkpoint(cp,ph,job); sim=PrecisionSimulation.restore(cp)
        if sim.precision_arm != 'phase_carry': raise ValueError('Carry arm changed')
        checked_history(history,p,job,sim.time-p['start'])
        if not np.array_equal(history[-1]['chemistry'],[host.activator,host.inhibitor]): raise ValueError('Checkpoint chemistry changed')
    else:
        sim=PrecisionSimulation(initialize(job['source'],job['chemical_file'],'uniform',job['point'],p['accepted_configs'][job['level']]),'phase_carry')
        history=[];audit=dict(max_volume_error=0.,min_radius=1e100,max_clipping=0.,boundary_max=0.,dilution_error_max=0.,wall_seconds=0.)
    start=_steps(p['start'],job['dt']); stop=_steps(p['start']+horizon,job['dt'])
    every=_steps(p['interval'],job['dt']); save_every=_steps(p['checkpoint_interval'],job['dt'])
    if not start<=sim.step_number<=_steps(p['start']+p['duration'],job['dt']) or abs(sim.time-sim.step_number*job['dt'])>1e-9:
        raise ValueError('Invalid carry checkpoint clock')
    if sim.step_number>stop:
        del sim;torch.cuda.empty_cache()
        return checked_history([r for r in history if r['elapsed']<=horizon+1e-9],p,job,horizon)
    clock=time.perf_counter();base_wall=audit['wall_seconds']
    while sim.step_number<=stop:
        quality=sim.audit()
        for key in ('max_volume_error','max_clipping'):audit[key]=max(audit[key],quality[key])
        audit['min_radius']=min(audit['min_radius'],quality['min_radius'])
        audit['dilution_error_max']=max(audit['dilution_error_max'],quality['dilution_amount_error'])
        if audit['dilution_error_max']>p['criteria']['dilution_error_max']:raise RuntimeError('Dilution failure')
        offset=sim.step_number-start
        if offset%every==0:
            elapsed=offset*job['dt']
            if not history or abs(history[-1]['elapsed']-elapsed)>1e-10:
                row=retain(sim.observe(elapsed),job);row['precision']=sim.precision_diagnostics();history.append(row)
            audit['boundary_max']=max(audit['boundary_max'],history[-1]['boundary_occupancy'])
            if audit['boundary_max']>=p['criteria']['boundary_max']:raise RuntimeError('Boundary failure')
            audit['wall_seconds']=base_wall+time.perf_counter()-clock
            write_json(folder/'history.json',history)
            write_json(folder/'status.json',dict(state='running',elapsed=elapsed,target_horizon=horizon,audit=audit))
            if offset%save_every==0 or sim.step_number==stop:
                save_checkpoint(sim,cp,audit,history,ph,job)
                print(f'{job["key"]}: elapsed={elapsed:g}, target={horizon:g}, wall={audit["wall_seconds"]:.1f}s',flush=True)
        if sim.step_number==stop:break
        sim.step()
    checked_history(history,p,job,horizon)
    if horizon==p['duration']:
        write_json(folder/'result.json',dict(job=job,protocol_sha256=ph,history_sha256=digest(folder/'history.json'),
            checkpoint_sha256=digest(cp),quality_pass=True,audit=audit,final_log_sd=history[-1]['log_activator_sd']))
        write_json(folder/'status.json',dict(state='completed',elapsed=horizon))
    else:write_json(folder/'status.json',dict(state='awaiting_pair_refinement',elapsed=horizon,audit=audit))
    del sim;torch.cuda.empty_cache();return history


def pair_report(root,p,seed,chi,horizon,histories):
    folder=Path(root)/'pairs'/f'seed-{seed}_chi-{chi:g}';folder.mkdir(parents=True,exist_ok=True)
    label='pilot' if horizon==p['pilot_duration'] else 'full'
    file=folder/f'{label}-refinement.json';ph=digest(Path(root)/'protocol.json')
    if file.exists():
        r=read(file)
        if r['protocol_sha256']!=ph:raise ValueError('Changed pair report protocol')
        for f,h in r['history_sha256'].items():
            if digest(f)!=h:raise ValueError('Changed immutable pair histories')
        for level in ('coarse','fine'):
            if read(folder/f'{label}-{level}.json')!=histories[level]:raise ValueError('Changed compared histories')
        return r
    files=[]
    for level in ('coarse','fine'):
        target=folder/f'{label}-{level}.json';write_json(target,histories[level]);files.append(target)
    r=dict(refinement_comparison(histories['coarse'],histories['fine'],p,horizon),seed=seed,chi=chi,protocol_sha256=ph,
        history_sha256={str(f.resolve()):digest(f) for f in files})
    write_json(file,r);return r


def initiation_comparison(zero,positive,numerical_pass):
    if not numerical_pass:return 'unresolved_numerical_agreement'
    if zero and not positive:return 'suppression_supported_in_this_history'
    if zero and positive:return 'both_initiate'
    if not zero and positive:return 'positive_contrast_only_initiates'
    return 'neither_initiates_no_successful_control'


def trajectory_summary(history,p):
    t=np.array([r['elapsed'] for r in history]);sd=np.array([r['log_activator_sd'] for r in history])
    growth=np.array([r['uniform_growth_max'] for r in history])
    late=float(sd[t>=p['duration']-p['late_window']-1e-9].min())
    return dict(persistent_contrast=late>p['criteria']['late_log_sd_min'],late_min_log_sd=late,final_log_sd=float(sd[-1]),
        onset_elapsed=first_crossing(t,sd,p['criteria']['late_log_sd_min'],'up'),
        uniform_growth_crossing=first_crossing(t,growth),
        positive_frozen_growth_integral=float(np.sum(.5*(np.maximum(growth[:-1],0)+np.maximum(growth[1:],0))*np.diff(t))),
        final_axis_ratio=history[-1]['axis_ratio'],initial_axis_ratio=history[0]['axis_ratio'])


def reused_endpoint(p,job):
    reuse=p['reused'][job['key']];folder=Path(reuse['root'])/reuse['job']['key']/'endpoint'
    ep=read(folder/'protocol.json');verify(ep);e=read(folder/'assay/result.json')
    if e['protocol_sha256']!=digest(folder/'protocol.json') or e['paths_sha256']!=digest(folder/'assay/paths.npz'):
        raise ValueError('Changed reused endpoint')
    return e


def assess(root):
    root=Path(root);p=read(root/'protocol.json');verify(p);histories={};rows=[];missing=[];endpoints=[]
    if not (root/'positive-reference-verification.json').exists():raise ValueError('Missing positive-chi reference gate')
    reference=positive_reference(root,p)
    for job in p['jobs']:
        if job['key'] not in p['reused'] and not (root/job['key']/'result.json').exists():missing.append(job['key']);continue
        h=advance(root,job,p,p['duration']);histories[job['seed'],job['point']['chi'],job['level']]=h
        e=reused_endpoint(p,job) if job['key'] in p['reused'] else endpoint_assay(root,job,p,h[-1])
        rows.append(dict(job=job,reused=job['key'] in p['reused'],**trajectory_summary(h,p)))
        endpoints.append(dict(job=job['key'],numerical_pass=e['numerical_pass'],all_trials_settled=e['all_trials_settled'],
            phase=e['phase'],local_bistability=e['local_bistability_supported']))
    pairs=[];by_history=[]
    for seed in p['histories']:
        numerical=True;branches={}
        for chi in p['contrasts']:
            if any((seed,chi,level) not in histories for level in ('coarse','fine')):numerical=False;continue
            h={level:histories[seed,chi,level] for level in ('coarse','fine')}
            r=pair_report(root,p,seed,chi,p['duration'],h);pairs.append(r);numerical &= r['passed']
            branches[chi]=trajectory_summary(h['fine'],p)
        by_history.append(dict(seed=seed,numerical_pass=bool(numerical),branches=branches,
            conclusion=initiation_comparison(branches.get(0.,{}).get('persistent_contrast',False),
                branches.get(.35,{}).get('persistent_contrast',False),numerical)))
    endpoint_pass=all(e['numerical_pass'] and e['all_trials_settled'] for e in endpoints)
    result=dict(protocol_sha256=digest(root/'protocol.json'),completed=sum(not r['reused'] for r in rows),total=10,
        independent_histories=3,new_histories=0,reused_trajectories=2,results=rows,missing=missing,
        full_refinement=pairs,by_history=by_history,endpoints=endpoints,endpoint_numerical_pass=bool(endpoint_pass),
        reference_pass=reference['passed'],quality_pass=not missing,
        passed=not missing and all(h['numerical_pass'] for h in by_history) and endpoint_pass,scope=p['interpretation'])
    write_json(root/'summary.json',result);return result


def positive_reference(root,p):
    root=Path(root);file=root/'positive-reference-verification.json';ph=digest(root/'protocol.json')
    target=root/'positive-float64-reference/result.json'
    if not file.exists():
        result=reference_run(target.parent,p['reference_contexts'],p['reference_criteria'],p['reference_horizon'],p['reference_dts'])
        write_json(file,dict(protocol_sha256=ph,result_sha256=digest(target),passed=result['passed']))
    gate=read(file);r=read(target)
    if (not gate['passed'] or gate['protocol_sha256']!=ph or gate['result_sha256']!=digest(target) or
            not r['passed'] or r['contexts']!=p['reference_contexts'] or r['criteria']!=p['reference_criteria'] or
            r['dts']!=p['reference_dts'] or r['horizon']!=p['reference_horizon']):raise ValueError('Positive-chi reference failed or changed')
    for f,h in r['field_sha256'].items():
        if digest(f)!=h:raise ValueError('Changed positive-chi reference fields')
    return r


def endpoint_task(root,job,p,final):
    return endpoint_assay(root,job,p,final)


def run_pairs(root,p,executor):
    """A failed pilot blocks its own long pair, without censoring other histories."""
    root=Path(root);done=0;pending=[];failed_pilots=[]
    for seed in (9,7,8):
        for chi in p['contrasts']:
            jobs=[j for j in p['jobs'] if j['seed']==seed and j['point']['chi']==chi]
            histories={}
            for job in jobs:
                write_json(root/'status.json',dict(state='running',stage='matched_60_unit_pilot',current_job=job['key'],completed=done,total=10))
                histories[job['level']]=advance(root,job,p,p['pilot_duration'])
            pilot=pair_report(root,p,seed,chi,p['pilot_duration'],histories)
            if not pilot['passed']:
                failed_pilots.append(dict(seed=seed,chi=chi));continue
            for job in jobs:
                write_json(root/'status.json',dict(state='running',stage='long_moving_initiation',current_job=job['key'],completed=done,total=10,
                    failed_pilot_contexts=failed_pilots))
                h=advance(root,job,p,p['duration'])
                if job['key'] not in p['reused']:
                    done+=1;pending.append(executor.submit(endpoint_task,root,job,p,h[-1]))
            pair_report(root,p,seed,chi,p['duration'],{j['level']:advance(root,j,p,p['duration']) for j in jobs})
    return done,pending,failed_pilots


def run(root):
    root=Path(root).resolve()
    with (root/'coordinator.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        p=read(root/'protocol.json');verify(p);torch.set_num_threads(1);device=p['device']
        if (device['name']!=torch.cuda.get_device_name(0) or device['torch']!=torch.__version__ or
                device['torch_cuda']!=torch.version.cuda or digest(library()._name)!=device['experimental_binary_sha256'] or
                digest(accepted_library()._name)!=device['accepted_binary_sha256'] or
                p['reference_runtime']!=dict(numpy=np.__version__,scipy=scipy.__version__)):
            raise ValueError('Changed qualified device/software')
        try:
            write_json(root/'status.json',dict(state='running',stage='positive_chi_float64_reference',completed=0,total=10))
            positive_reference(root,p)
            with ProcessPoolExecutor(max_workers=p['endpoint_workers'],mp_context=multiprocessing.get_context('spawn')) as executor:
                done,pending,failed_pilots=run_pairs(root,p,executor)
                write_json(root/'status.json',dict(state='running',stage='endpoint_assays',completed=done,total=10,failed_pilot_contexts=failed_pilots))
                for future in pending:future.result()
            summary=assess(root)
            write_json(root/'status.json',dict(state='completed' if summary['passed'] else 'completed_with_unresolved_checks',
                completed=summary['completed'],total=10,passed=summary['passed'],failed_pilot_contexts=failed_pilots))
            return summary
        except Exception as error:
            write_json(root/'status.json',dict(state='failed',error=str(error)));raise


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=('prepare','run','assess'))
    parser.add_argument('--output',type=Path,default=Path('outputs/phase-carry-polarity'))
    args=parser.parse_args();result=prepare(args.output) if args.action=='prepare' else run(args.output) if args.action=='run' else assess(args.output)
    print({k:v for k,v in result.items() if k in ('completed','total','passed','estimated_new_moving_seconds','by_history')})
