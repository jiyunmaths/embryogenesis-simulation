"""Matched moving dilution/contact controls conditional on carry nonformation."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import fcntl
import multiprocessing
from pathlib import Path
import shutil
import time

import numpy as np
import torch

from .feedback_long import digest, save_checkpoint, restore_checkpoint
from .geometry_precision import checked_history
from .gpu_initiation_controls import ARMS, InitiationControlSimulation, checkpoint_fields
from .gpu_precision_control import PrecisionSimulation, library
from .neighbor_context import read, validate_graph
from .parameter_robustness import verify
from .parameter_robustness_moving import initialize, retain, endpoint_assay
from .phase_carry_polarity import current_evidence, trajectory_summary, positive_reference
from .polarity_robustness import refinement_comparison
from .resolution import _steps, write_json


def job_design(sources):
    return [dict(key=f'seed-{seed}_{level}_{arm}', seed=seed, level=level, dt=dt,
                 arm=arm, family='uniform', point=dict(key='polarity-0', ratio=27.5, chi=0.),
                 source=sources[seed][level], chemical_file=sources[seed]['chemical'],
                 initial_transport=sources[seed]['transport'])
            for seed in (7, 8) for arm in ARMS
            for level, dt in (('coarse', .00375), ('fine', .001875))]


def rescue_interpretation(forms, qualified):
    if not qualified:
        return 'unresolved_numerical_or_physical_checks'
    if forms['baseline']:
        return 'baseline_forms_not_a_failure_rescue'
    d, g, both = [forms[k] for k in ('dilution-off', 'fixed-conductances', 'combined')]
    if not both and (d or g):
        return 'single_arm_rescue_but_combined_nonrescue_investigate_nonmonotonicity'
    if not both:
        return 'dilution_and_contact_ablation_insufficient'
    if d and g:
        return 'each_single_intervention_restores_formation'
    if d:
        return 'removing_dilution_restores_formation'
    if g:
        return 'fixing_contact_conductances_restores_formation'
    return 'joint_removal_required_in_tested_context'


def cpu_chemical_step(state, delta, dt, beta, da, db):
    """Independent NumPy evaluation of the declared SSP-RK2 chemical update."""
    state = np.asarray(state, dtype=np.float64).copy()
    rate = max(1+da*max(-np.diag(delta)), beta+db*max(-np.diag(delta)))
    count = max(1, int(np.ceil(dt*rate/.2)))
    step = dt/count
    def reaction_transport(y):
        a, b = y
        return np.array([a*a/b-a+da*(delta@a), beta*(a*a-b)+db*(delta@b)])
    for _ in range(count):
        first = state+step*reaction_transport(state)
        state = .5*state+.5*(first+step*reaction_transport(first))
    return state


def exact_state(a, b, controls=False):
    keys = ['phi', 'activator', 'inhibitor', 'polarity', 'geometry', 'phase_carry', 'rounding']
    if controls:
        keys += ['initial_conductance', 'initial_volume', 'amount_source', 'last_amount_source',
                 'last_relative_volume_change']
    if any(not torch.equal(getattr(a, k), getattr(b, k)) for k in keys):
        raise ValueError('Control baseline/restart state differs')
    if (a.time != b.time or a.step_number != b.step_number or
            not torch.equal(a.matrices()[1], b.matrices()[1])):
        raise ValueError('Control baseline/restart clock or geometry differs')


def implementation_gate(folder, source, chemical, config, initial_transport):
    folder.mkdir(parents=True)
    make = lambda: initialize(source, chemical, 'uniform', dict(ratio=27.5, chi=0.), config)
    baseline = InitiationControlSimulation(make())
    original = PrecisionSimulation(make(), 'phase_carry')
    with np.load(initial_transport) as z:
        if (not np.array_equal(baseline.initial_conductance.cpu().numpy(), z['conductance']) or
                not np.array_equal(baseline.initial_volume.cpu().numpy(), z['volumes']) or
                not np.array_equal(baseline.ids, z['ids'])):
            raise ValueError('Changed actual initial transport')
    for _ in range(16):
        baseline.step(); original.step(); exact_state(baseline, original)
    del baseline, original
    rows = []
    for arm in ARMS:
        sim = InitiationControlSimulation(make(), arm)
        error, accounting = 0., 0.
        for _ in range(8):
            _, geometric, _, g = sim.matrices()
            delta, _ = sim.chemical_transport(geometric, g)
            before = np.array([sim.activator.cpu().numpy(), sim.inhibitor.cpu().numpy()])
            old = sim.geometry[:, 0].cpu().numpy().copy()
            expected = cpu_chemical_step(before, delta.cpu().numpy(), sim.config.dt, 2., .02, .55)
            sim.step()
            new = sim.geometry[:, 0].cpu().numpy().copy()
            source_term = expected*(new-old)[None, :] if not sim.dilution else np.zeros_like(expected)
            if sim.dilution:
                expected *= old/new
            actual = np.array([sim.activator.cpu().numpy(), sim.inhibitor.cpu().numpy()])
            error = max(error, float(abs(np.log(actual/expected)).max()))
            np.testing.assert_allclose(sim.last_amount_source.cpu().numpy(), source_term, rtol=1e-10, atol=1e-14)
            accounting = max(accounting, sim.audit()['volume_conversion_amount_error'])
            validate_graph(delta.cpu().numpy(), old)
        if error > 1e-11 or accounting > 2e-14:
            raise ValueError('Independent chemical/accounting gate failed')
        path = folder/f'{arm}-restart.npz'; sim.checkpoint(path)
        restart = InitiationControlSimulation.restore(path)
        exact_state(sim, restart, True)
        for _ in range(2):
            sim.step(); restart.step(); exact_state(sim, restart, True)
            if not torch.equal(sim._conversion_error, restart._conversion_error):
                raise ValueError('Accounting state differs on restart')
        rows.append(dict(arm=arm, passed=True, cpu_chemical_log_max=error,
                         amount_accounting_error_max=accounting, restart_exact_steps=2,
                         checkpoint=str(path.resolve()), checkpoint_sha256=digest(path)))
        del sim, restart
    torch.cuda.empty_cache()
    return dict(passed=True, baseline_exact_steps=16, controls=rows)


def prepare(root, parent=Path('outputs/phase-carry-polarity'), review=Path('docs/phase_carry_polarity_assessment.json')):
    root, parent, review = [Path(x).resolve() for x in (root, parent, review)]
    if root.exists():
        raise FileExistsError(root)
    old = read(parent/'protocol.json'); verify(old)
    summary = read(parent/'summary.json'); evidence = read(review)
    if (read(parent/'status.json')['state'] != 'completed' or not summary['passed'] or
            summary['protocol_sha256'] != digest(parent/'protocol.json') or
            evidence['protocol_sha256'] != digest(parent/'protocol.json') or
            evidence['original_summary_sha256'] != digest(parent/'summary.json') or
            evidence['followup_eligible_histories'] != [7, 8]):
        raise ValueError('Both histories must have qualified moving failure and matched frozen formation')
    if digest(evidence['full_review']) != evidence['full_review_sha256']:
        raise ValueError('Changed reviewed evidence')
    actual = read(evidence['full_review'])
    for row in actual['by_history']:
        if row['seed'] in (7, 8):
            ref = row['frozen_reference']
            if (not row['qualified'] or not row['moving_vs_frozen_failure_without_directional_tension'] or
                    not ref['persistent_contrast'] or ref['solver_log_error'] > 1e-8 or
                    digest(ref['saved_path']) != ref['saved_path_sha256']):
                raise ValueError('Unqualified actual-input frozen formation reference')
    positive_reference(parent, old)
    if (torch.cuda.get_device_name(0) != old['device']['name'] or
            torch.__version__ != old['device']['torch'] or torch.version.cuda != old['device']['torch_cuda'] or
            digest(library()._name) != old['device']['experimental_binary_sha256']):
        raise ValueError('Changed qualified GPU/runtime/kernel')
    root.mkdir(parents=True)
    sources, reused, gates = {}, {}, []
    inputs = [parent/'protocol.json', parent/'summary.json', parent/'status.json', review,
              Path(evidence['full_review']), Path('docs/phase_carry_polarity_reporting.json'),
              Path('docs/moving_geometry_initiation_controls.md')]
    for seed in (7, 8):
        folder = root/f'seed-{seed}'; folder.mkdir()
        sources[seed] = {}
        for level, name in (('coarse','coarse-source.npz'), ('fine','fine-source.npz'), ('chemical','initial_states.npz')):
            original, copied = parent/f'seed-{seed}'/name, folder/name
            shutil.copy2(original, copied); inputs += [original, copied]
            sources[seed][level] = str(copied)
        transport = folder/'initial-transport.npz'
        made = []
        for level in ('coarse', 'fine'):
            job = next(j for j in old['jobs'] if j['seed'] == seed and j['level'] == level and j['point']['chi'] == 0.)
            r, h = current_evidence(parent, old, job)
            if trajectory_summary(h, old)['persistent_contrast']:
                raise ValueError('Baseline forms: this is not a failure rescue experiment')
            reused[f'seed-{seed}_{level}_baseline'] = dict(root=str(parent), job=job)
            inputs += [parent/job['key']/name for name in ('result.json','history.json','latest_state.npz',
                         'context-gate.json','endpoint/protocol.json','endpoint/assay/result.json','endpoint/assay/paths.npz')]
            sim = InitiationControlSimulation(initialize(sources[seed][level], sources[seed]['chemical'], 'uniform',
                       job['point'], old['accepted_configs'][level]))
            _, delta, _, g = sim.matrices()
            np.testing.assert_array_equal(delta.cpu().numpy(), h[0]['delta'])
            np.testing.assert_array_equal(sim.geometry[:, 0].cpu().numpy(), h[0]['volumes'])
            np.testing.assert_array_equal([sim.activator.cpu().numpy(), sim.inhibitor.cpu().numpy()], h[0]['chemistry'])
            made.append(g.cpu().numpy().copy())
            if level == 'coarse':
                np.savez_compressed(transport, conductance=made[-1], volumes=h[0]['volumes'], delta=h[0]['delta'], ids=sim.ids)
            del sim
        np.testing.assert_array_equal(*made)
        sources[seed]['transport'] = str(transport); inputs.append(transport)
        for level in ('coarse','fine'):
            gate_folder = folder/f'gate-{level}'
            gates.append(dict(seed=seed, level=level, **implementation_gate(gate_folder, sources[seed][level],
                         sources[seed]['chemical'], old['accepted_configs'][level], transport)))
            inputs += [Path(t['checkpoint']) for t in gates[-1]['controls']]
    gate_file = root/'implementation-verification.json'
    write_json(gate_file, dict(passed=True, contexts=gates)); inputs.append(gate_file)
    jobs = job_design(sources)
    source_hashes = {**old['source_sha256'], **{str(Path(__file__).with_name(name).resolve()):
        digest(Path(__file__).with_name(name)) for name in ('moving_initiation_controls.py','gpu_initiation_controls.py')}}
    prior_durations = [read(parent/item['job']['key']/'result.json')['audit']['wall_seconds'] for item in reused.values()]
    p = dict(jobs=jobs, reused=reused, parent=str(parent), initial_transport_files={str(k):v['transport'] for k,v in sources.items()},
        accepted_configs=old['accepted_configs'], device=old['device'], histories=[7,8], arms=list(ARMS),
        start=150., duration=240., pilot_duration=60., late_window=24., pilot_late_window=12.,
        interval=.15, checkpoint_interval=3., independent_histories=2, new_histories=0,
        new_moving_jobs=12, reused_trajectories=4, endpoint_workers=2,
        criteria=old['criteria'], refinement_criteria=old['refinement_criteria'],
        endpoint_horizons=old['endpoint_horizons'], endpoint_interval=old['endpoint_interval'], endpoint_criteria=old['endpoint_criteria'],
        chemical_reference_log_max=1e-11, amount_accounting_max=2e-14,
        estimated_new_moving_seconds=3*sum(prior_durations),
        source_sha256=source_hashes, input_sha256={**old['input_sha256'], **{str(f.resolve()):digest(f) for f in inputs}},
        design='Two existing eligible histories 7/8, mature t=150 near-uniform starts, chi=0; beta=2, Da=.02, Db=.55, activity tension/adhesion .25/.35. Four factorial arms, two timesteps, 240 elapsed units. Reuse four qualified baseline paths; twelve new moving paths. Mechanics, polarity, volume and chemical feedback evolve in every arm.',
        interventions='Freeze actual initial symmetric conductances G0, retain changing masses: Delta=-M(t)^(-1)K0. Never freeze Delta0. Dilution off omits concentration conversion and records source per cell (Vnew-Vold)c after reaction/transport. On arms conserve conversion amounts; off arms match the explicit source to 2e-14. Chemical and geometric operators saved separately.',
        interpretation='Primary endpoint: sustained SD(log activator)>.1 throughout elapsed216-240, with original raw refinement and physical checks. Failed pilots block only their own pair. Spectra and frozen endpoints use the actual chemical operator, not geometric transport in fixed-contact arms. Different rescues imply different contributions; identical rescue patterns support a common tested contribution, not a unique universal cause. No new history, fresh zygote, inheritance/autonomy, spatial/closure/backend promotion. Capacity clamp and reconditioned geometry are subsequent diagnostics, not automatically launched.',
        history9_positive_benchmark=dict(root=str(parent), keys=[j['key'] for j in old['jobs'] if j['seed']==9 and j['point']['chi']==0.]))
    write_json(root/'protocol.json', p)
    write_json(root/'status.json', dict(state='prepared', completed=0, total=12, reused_trajectories=4,
                                      estimated_new_moving_seconds=p['estimated_new_moving_seconds']))
    return p


def validate_control_history(h, p, job, horizon):
    checked_history(h, p, job, horizon)
    with np.load(job['initial_transport']) as z:
        initial = z['conductance']
    if not np.allclose(np.array(h[0]['chemical_conductance']), initial, rtol=1e-14, atol=1e-14):
        raise ValueError('Changed initial conductances')
    for row in h:
        if row['initiation_arm'] != job['arm']:
            raise ValueError('History intervention changed')
        volume = np.asarray(row['volumes'])
        for key, gkey in (('delta','chemical_conductance'), ('geometric_delta','geometric_conductance')):
            delta, _ = validate_graph(row[key], volume)
            g = np.array(row[gkey]); expected = (g-np.diag(g.sum(1)))/volume[:,None]
            if not np.allclose(delta, expected, rtol=1e-12, atol=1e-13):
                raise ValueError('Conductance/operator mismatch')
        if ARMS[job['arm']][1] and not np.array_equal(row['chemical_conductance'], initial):
            raise ValueError('Fixed conductances changed')
        if ARMS[job['arm']][0] and np.any(row['cumulative_volume_amount_source']):
            raise ValueError('Unexpected dilution-on source')
    return h


def baseline_history(p, job, horizon):
    item = p['reused'][job['key']]; root = Path(item['root']); old = read(root/'protocol.json')
    _, h = current_evidence(root, old, item['job'])
    result=[]
    for raw in h:
        if raw['elapsed'] > horizon+1e-9: break
        row=dict(raw); g=np.array(raw['delta'])*np.array(raw['volumes'])[:,None]; np.fill_diagonal(g,0)
        row.update(geometric_delta=raw['delta'], chemical_conductance=g.tolist(), geometric_conductance=g.tolist(),
                   cumulative_volume_amount_source=np.zeros((2,len(raw['ids']))).tolist(), initiation_arm='baseline')
        result.append(row)
    return validate_control_history(result,p,job,horizon)


def completed_evidence(root,p,job):
    folder=Path(root)/job['key']; ph=digest(Path(root)/'protocol.json')
    r,h=read(folder/'result.json'),read(folder/'history.json')
    if (r['job']!=job or r['protocol_sha256']!=ph or not r['quality_pass'] or
            r['history_sha256']!=digest(folder/'history.json') or r['checkpoint_sha256']!=digest(folder/'latest_state.npz')):
        raise ValueError('Changed completed control evidence')
    host,audit,saved=restore_checkpoint(folder/'latest_state.npz',ph,job)
    if (saved!=h or audit!=r['audit'] or abs(host.time-p['start']-p['duration'])>1e-9 or
            not np.array_equal([host.activator,host.inhibitor],h[-1]['chemistry'])):
        raise ValueError('Control checkpoint/history mismatch')
    with np.load(folder/'latest_state.npz') as z:
        checkpoint_fields(z,host.phi.shape)
        if str(z['initiation_arm'])!=job['arm'] or not np.array_equal(z['amount_source'],h[-1]['cumulative_volume_amount_source']):
            raise ValueError('Control accounting/arm mismatch')
        with np.load(job['initial_transport']) as original:
            if not np.array_equal(z['initial_conductance'],original['conductance']):
                raise ValueError('Initial conductances lost on restart')
    validate_control_history(h,p,job,p['duration'])
    return r,h


def advance(root,job,p,horizon):
    if horizon not in (p['pilot_duration'],p['duration']): raise ValueError('Undeclared horizon')
    if job['key'] in p['reused']: return baseline_history(p,job,horizon)
    root=Path(root);folder=root/job['key'];folder.mkdir(exist_ok=True);ph=digest(root/'protocol.json')
    cp=folder/'latest_state.npz'
    if (folder/'result.json').exists():
        _,h=completed_evidence(root,p,job)
        return validate_control_history([r for r in h if r['elapsed']<=horizon+1e-9],p,job,horizon)
    if cp.exists():
        host,audit,history=restore_checkpoint(cp,ph,job);sim=InitiationControlSimulation.restore(cp)
        if sim.initiation_arm!=job['arm']: raise ValueError('Restart intervention changed')
        validate_control_history(history,p,job,sim.time-p['start'])
        if not np.array_equal(history[-1]['chemistry'],[host.activator,host.inhibitor]): raise ValueError('Restart chemistry mismatch')
    else:
        sim=InitiationControlSimulation(initialize(job['source'],job['chemical_file'],'uniform',job['point'],p['accepted_configs'][job['level']]),job['arm'])
        history=[];audit=dict(max_volume_error=0.,min_radius=1e100,max_clipping=0.,boundary_max=0.,dilution_error_max=0.,wall_seconds=0.)
    with np.load(job['initial_transport']) as z:
        np.testing.assert_array_equal(sim.initial_conductance.cpu().numpy(),z['conductance'])
    start=_steps(p['start'],job['dt']);stop=_steps(p['start']+horizon,job['dt'])
    every=_steps(p['interval'],job['dt']);save_every=_steps(p['checkpoint_interval'],job['dt'])
    if not start<=sim.step_number<=_steps(p['start']+p['duration'],job['dt']) or abs(sim.time-sim.step_number*job['dt'])>1e-9:
        raise ValueError('Invalid control clock')
    if sim.step_number>stop:
        del sim;torch.cuda.empty_cache()
        return validate_control_history([r for r in history if r['elapsed']<=horizon+1e-9],p,job,horizon)
    clock=time.perf_counter();base=audit['wall_seconds']
    while sim.step_number<=stop:
        quality=sim.audit()
        for key in ('max_volume_error','max_clipping'): audit[key]=max(audit[key],quality[key])
        audit['min_radius']=min(audit['min_radius'],quality['min_radius'])
        audit['dilution_error_max']=max(audit['dilution_error_max'],quality['volume_conversion_amount_error'])
        if audit['dilution_error_max']>p['amount_accounting_max']: raise RuntimeError('Volume-conversion accounting failure')
        offset=sim.step_number-start
        if offset%every==0:
            elapsed=offset*job['dt']
            if not history or abs(history[-1]['elapsed']-elapsed)>1e-10:
                row=retain(sim.observe(elapsed),job);row['precision']=sim.precision_diagnostics();history.append(row)
            audit['boundary_max']=max(audit['boundary_max'],history[-1]['boundary_occupancy'])
            if audit['boundary_max']>=p['criteria']['boundary_max']: raise RuntimeError('Boundary failure')
            audit['wall_seconds']=base+time.perf_counter()-clock
            write_json(folder/'status.json',dict(state='running',elapsed=elapsed,target_horizon=horizon,audit=audit,
                       log_activator_sd=history[-1]['log_activator_sd']))
            # Preserve every 0.15-unit sample, publish the full history only at
            # checkpoints; avoid repeated growing JSON serialization each sample.
            if offset%save_every==0 or sim.step_number==stop:
                save_checkpoint(sim,cp,audit,history,ph,job);write_json(folder/'history.json',history)
                print(f'{job["key"]}: elapsed={elapsed:g}, target={horizon:g}, wall={audit["wall_seconds"]:.1f}s',flush=True)
        if sim.step_number==stop: break
        sim.step()
    validate_control_history(history,p,job,horizon)
    if horizon==p['duration']:
        write_json(folder/'result.json',dict(job=job,protocol_sha256=ph,history_sha256=digest(folder/'history.json'),
                   checkpoint_sha256=digest(cp),quality_pass=True,audit=audit,**trajectory_summary(history,p)))
        write_json(folder/'status.json',dict(state='completed',elapsed=horizon))
    else: write_json(folder/'status.json',dict(state='awaiting_pair_refinement',elapsed=horizon,audit=audit))
    del sim;torch.cuda.empty_cache();return history


def pair_report(root,p,seed,arm,horizon,histories):
    folder=Path(root)/'pairs'/f'seed-{seed}_{arm}';folder.mkdir(parents=True,exist_ok=True)
    label='pilot' if horizon==p['pilot_duration'] else 'full';file=folder/f'{label}-refinement.json'
    ph=digest(Path(root)/'protocol.json')
    result=refinement_comparison(histories['coarse'],histories['fine'],p,horizon)
    geometric=max(float(np.linalg.norm(np.array(a['geometric_delta'])-b['geometric_delta'])/
        max(np.linalg.norm(b['geometric_delta']),1e-30)) for a,b in zip(histories['coarse'],histories['fine']))
    result['geometric_relative_transport_max']=geometric
    result['passed'] &= geometric<=p['refinement_criteria']['relative_transport_max']
    result['interpretation']=f'Same-start raw timestep halving in history {seed}, {arm}; spectra are frozen diagnostics.'
    if file.exists():
        saved=read(file)
        if saved['protocol_sha256']!=ph: raise ValueError('Changed comparison protocol')
        for f,h in saved['history_sha256'].items():
            if digest(f)!=h: raise ValueError('Changed immutable comparison histories')
        for key,value in result.items():
            if saved[key]!=value: raise ValueError('Changed comparison metrics')
        return saved
    files=[]
    for level in ('coarse','fine'):
        target=folder/f'{label}-{level}.json';write_json(target,histories[level]);files.append(target)
    result.update(seed=seed,arm=arm,protocol_sha256=ph,history_sha256={str(f.resolve()):digest(f) for f in files})
    write_json(file,result);return result


def run_pairs(root,p,executor):
    done,pending,failed=0,[],[]
    for seed in p['histories']:
        for arm in p['arms']:
            jobs=[j for j in p['jobs'] if j['seed']==seed and j['arm']==arm]
            h={}
            for job in jobs:
                write_json(Path(root)/'status.json',dict(state='running',stage='matched_pilot',current_job=job['key'],completed=done,total=12))
                h[job['level']]=advance(root,job,p,p['pilot_duration'])
            if not pair_report(root,p,seed,arm,p['pilot_duration'],h)['passed']:
                failed.append(dict(seed=seed,arm=arm));continue
            for job in jobs:
                write_json(Path(root)/'status.json',dict(state='running',stage='long_moving_controls',current_job=job['key'],completed=done,total=12,failed_pilot_contexts=failed))
                h[job['level']]=advance(root,job,p,p['duration'])
                if job['key'] not in p['reused']:
                    done+=1;pending.append(executor.submit(endpoint_assay,root,job,p,h[job['level']][-1]))
            pair_report(root,p,seed,arm,p['duration'],h)
    return done,pending,failed


def assess(root,p):
    verify(p)
    rows,endpoints,missing=[],[],[]
    for job in p['jobs']:
        if job['key'] in p['reused']:
            h=baseline_history(p,job,p['duration']);item=p['reused'][job['key']]
            e=read(Path(item['root'])/item['job']['key']/'endpoint/assay/result.json')
        elif not (Path(root)/job['key']/'result.json').exists():
            missing.append(job['key']);continue
        else:
            _,h=completed_evidence(root,p,job);e=endpoint_assay(root,job,p,h[-1])
        rows.append(dict(job=job,reused=job['key'] in p['reused'],**trajectory_summary(h,p)))
        endpoints.append(dict(job=job['key'],numerical_pass=e['numerical_pass'],all_trials_settled=e['all_trials_settled'],phase=e['phase'],local_bistability=e['local_bistability_supported']))
    by_history=[]
    for seed in p['histories']:
        history_endpoints=[e for e in endpoints if next(j for j in p['jobs'] if j['key']==e['job'])['seed']==seed]
        qualified=len(history_endpoints)==8 and all(e['numerical_pass'] and e['all_trials_settled'] for e in history_endpoints)
        forms={}
        for arm in p['arms']:
            file=Path(root)/'pairs'/f'seed-{seed}_{arm}'/'full-refinement.json'
            if not file.exists(): qualified=False;continue
            pair=read(file)
            if pair['protocol_sha256']!=digest(Path(root)/'protocol.json'):
                raise ValueError('Changed full comparison protocol')
            for f,h in pair['history_sha256'].items():
                if digest(f)!=h: raise ValueError('Changed full comparison evidence')
            histories={level:read(file.parent/f'full-{level}.json') for level in ('coarse','fine')}
            pair_report(root,p,seed,arm,p['duration'],histories)
            qualified &= pair['passed']
            matching=[r for r in rows if r['job']['seed']==seed and r['job']['arm']==arm]
            if len(matching)!=2: qualified=False;continue
            forms[arm]=next(r['persistent_contrast'] for r in matching if r['job']['level']=='fine')
        by_history.append(dict(seed=seed,numerical_pass=bool(qualified),formation= forms,
                               conclusion=rescue_interpretation(forms,qualified)))
    passed=not missing and all(r['numerical_pass'] for r in by_history) and all(e['numerical_pass'] and e['all_trials_settled'] for e in endpoints)
    result=dict(protocol_sha256=digest(Path(root)/'protocol.json'),passed=bool(passed),completed=sum(not r['reused'] for r in rows),total=12,
                independent_histories=2,new_histories=0,reused_trajectories=4,results=rows,missing=missing,by_history=by_history,endpoints=endpoints,scope=p['interpretation'])
    write_json(Path(root)/'summary.json',result);return result


def run(root):
    root=Path(root).resolve()
    with (root/'coordinator.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        p=read(root/'protocol.json');verify(p);torch.set_num_threads(1)
        if (torch.cuda.get_device_name(0)!=p['device']['name'] or torch.__version__!=p['device']['torch'] or
                torch.version.cuda!=p['device']['torch_cuda'] or digest(library()._name)!=p['device']['experimental_binary_sha256']):
            raise ValueError('Changed qualified GPU/runtime/kernel')
        try:
            with ProcessPoolExecutor(max_workers=p['endpoint_workers'],mp_context=multiprocessing.get_context('spawn')) as executor:
                done,pending,failed=run_pairs(root,p,executor)
                write_json(root/'status.json',dict(state='running',stage='endpoint_assays',completed=done,total=12,failed_pilot_contexts=failed))
                for future in pending: future.result()
            result=assess(root,p)
            write_json(root/'status.json',dict(state='completed' if result['passed'] else 'completed_with_unresolved_checks',
                                              completed=result['completed'],total=12,passed=result['passed'],failed_pilot_contexts=failed))
            return result
        except Exception as error:
            write_json(root/'status.json',dict(state='failed',error=str(error)));raise


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=('prepare','run','assess'))
    parser.add_argument('--output',type=Path,default=Path('outputs/moving-initiation-controls'))
    args=parser.parse_args()
    result=prepare(args.output) if args.action=='prepare' else run(args.output) if args.action=='run' else assess(args.output,read(args.output/'protocol.json'))
    print({k:v for k,v in result.items() if k in ('passed','completed','total','estimated_new_moving_seconds','by_history')})
