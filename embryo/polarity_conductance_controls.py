"""Matched polarity/contact intervention, reusing qualified carry trajectories.

The backend and all previously pinned scientific sources remain unchanged.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import fcntl
import multiprocessing
from pathlib import Path
import shutil
import time

import numpy as np
import torch

from . import moving_initiation_controls as contact
from .feedback_long import digest, save_checkpoint, restore_checkpoint
from .gpu_initiation_controls import InitiationControlSimulation
from .gpu_precision_control import PrecisionSimulation, library
from .neighbor_context import read, validate_graph
from .parameter_robustness import verify
from .parameter_robustness_moving import initialize, retain, endpoint_assay
from .phase_carry_polarity import current_evidence, trajectory_summary, initiation_comparison
from .resolution import _steps, write_json

HISTORIES = (9, 7, 8)  # First finish the original polarity-suppression history.
CONTRASTS = (0., .35)
ARMS = ('baseline', 'fixed-conductances')


def job_design(sources):
    return [dict(key=f'seed-{seed}_{level}_chi-{chi:g}_{arm}', seed=seed,
                 level=level, dt=dt, arm=arm, family='uniform',
                 point=dict(key=f'polarity-{chi:g}', ratio=27.5, chi=chi),
                 source=sources[seed][level], chemical_file=sources[seed]['chemical'],
                 initial_transport=sources[seed]['transport'])
            for seed in HISTORIES for chi in CONTRASTS for arm in ARMS
            for level, dt in (('coarse', .00375), ('fine', .001875))]


def interpretation(forms, qualified):
    """Report intervention outcomes separately from original suppression evidence."""
    if not qualified:
        return dict(qualified=False, conclusion='unresolved_numerical_or_physical_checks')
    e0, ep = forms['chi-0_baseline'], forms['chi-0.35_baseline']
    f0, fp = forms['chi-0_fixed-conductances'], forms['chi-0.35_fixed-conductances']
    original = initiation_comparison(e0, ep, True)
    fixed = initiation_comparison(f0, fp, True)
    if not f0:
        conclusion = ('positive_contrast_only_forms_with_fixed_conductances' if fp else
                      'neither_fixed_conductance_branch_forms_no_successful_control')
    elif not fp:
        conclusion = 'polarity_suppression_persists_with_fixed_conductances'
    elif ep:
        conclusion = 'both_fixed_branches_form_positive_baseline_already_forms'
    elif e0:
        conclusion = 'conductance_preservation_removes_observed_polarity_suppression'
    else:
        conclusion = 'fixed_conductances_enable_both_contrasts_without_original_polarity_specific_loss'
    return dict(qualified=True, original_polarity_comparison=original,
                fixed_conductance_polarity_comparison=fixed,
                positive_contrast_conductance_rescue=bool(not ep and fp), conclusion=conclusion)


def context_key(job):
    return f'chi-{job["point"]["chi"]:g}_{job["arm"]}'


def check_reuse_context(job, original, kind):
    if (kind not in ('phase', 'contact') or
            any(job[k] != original[k] for k in ('seed', 'level', 'dt', 'family')) or
            any(job['point'][k] != original['point'][k] for k in ('ratio', 'chi')) or
            (kind == 'phase' and (job['arm'] != 'baseline' or original['arm'] != 'phase_carry')) or
            (kind == 'contact' and job['arm'] != original['arm'])):
        raise ValueError('Reused trajectory is not the declared matched intervention')


def check_config(config, job, p):
    expected = dict(p['accepted_configs'][job['level']])
    expected.update(seed=job['seed'], signal_dh=.02*job['point']['ratio'], polarity_tension=job['point']['chi'])
    if config != expected or config['dt'] != job['dt']:
        raise ValueError('Checkpoint physical parameters differ from declared job')


def validate_history(history, p, job, horizon):
    # The older helper requires a positive duration. A saved initial checkpoint
    # is also a valid restart, provided its entire initial observation is checked.
    if horizon == 0.:
        if len(history) != 1 or history[0]['elapsed'] != 0. or history[0]['time'] != p['start']:
            raise ValueError('Invalid initial checkpoint clock')
        with np.load(job['chemical_file']) as z:
            if not np.array_equal(history[0]['chemistry'], z['uniform']) or not np.array_equal(history[0]['ids'], z['ids']):
                raise ValueError('Changed initial chemical state or cell order')
        for key in ('delta', 'geometric_delta'):
            validate_graph(history[0][key], history[0]['volumes'])
        with np.load(job['initial_transport']) as z:
            for key, field in (('chemical_conductance', 'conductance'), ('volumes', 'volumes')):
                np.testing.assert_allclose(history[0][key], z[field], rtol=1e-14, atol=1e-14)
            expected = (z['conductance']-np.diag(z['conductance'].sum(1)))/z['volumes'][:,None]
            for key in ('delta','geometric_delta'):
                np.testing.assert_allclose(history[0][key],expected,rtol=1e-12,atol=1e-13)
            np.testing.assert_allclose(history[0]['geometric_conductance'],z['conductance'],rtol=1e-14,atol=1e-14)
        if history[0]['initiation_arm'] != job['arm'] or np.any(history[0]['cumulative_volume_amount_source']):
            raise ValueError('Invalid initial intervention/accounting')
        return history
    return contact.validate_control_history(history, p, job, horizon)


def reused_evidence(p, job):
    item = p['reused'][job['key']]
    check_reuse_context(job, item['job'], item['kind'])
    root = Path(item['root']); old = read(root/'protocol.json')
    if digest(root/'protocol.json') != item['protocol_sha256']:
        raise ValueError('Changed reused protocol')
    r, raw = (current_evidence(root, old, item['job']) if item['kind'] == 'phase' else
              contact.completed_evidence(root, old, item['job']))
    history = []
    for record in raw:
        row = dict(record)
        if item['kind'] == 'phase':
            g = np.array(row['delta'])*np.array(row['volumes'])[:, None]
            np.fill_diagonal(g, 0.)
            row.update(geometric_delta=row['delta'], chemical_conductance=g.tolist(),
                       geometric_conductance=g.tolist(), initiation_arm='baseline',
                       cumulative_volume_amount_source=np.zeros((2,len(row['ids']))).tolist())
        history.append(row)
    validate_history(history, p, job, p['duration'])
    return r, history


def completed_evidence(root, p, job):
    result, history = contact.completed_evidence(root, p, job)
    import json
    with np.load(Path(root)/job['key']/'latest_state.npz') as z:
        check_config(json.loads(str(z['metadata']))['config'], job, p)
    return result, history


def actual_context_gate(folder, source, chemical, point, config, transport):
    folder.mkdir(parents=True)
    make = lambda: initialize(source, chemical, 'uniform', point, config)
    original = PrecisionSimulation(make(), 'phase_carry')
    baseline = InitiationControlSimulation(make())
    for _ in range(16):
        original.step(); baseline.step(); contact.exact_state(original, baseline)
    del original, baseline
    rows = []
    for arm in ARMS:
        sim = InitiationControlSimulation(make(), arm)
        with np.load(transport) as z:
            np.testing.assert_array_equal(sim.initial_conductance.cpu().numpy(), z['conductance'])
            np.testing.assert_array_equal(sim.initial_volume.cpu().numpy(), z['volumes'])
            np.testing.assert_array_equal(sim.ids, z['ids'])
        error, accounting = 0., 0.
        for _ in range(8):
            _, geometric, _, g = sim.matrices()
            delta, _ = sim.chemical_transport(geometric, g)
            state = np.array([sim.activator.cpu().numpy(), sim.inhibitor.cpu().numpy()])
            old = sim.geometry[:,0].cpu().numpy().copy()
            expected = contact.cpu_chemical_step(state, delta.cpu().numpy(), sim.config.dt, 2., .02, .55)
            validate_graph(delta.cpu().numpy(), old)
            sim.step(); new = sim.geometry[:,0].cpu().numpy().copy()
            expected *= old/new
            actual = np.array([sim.activator.cpu().numpy(), sim.inhibitor.cpu().numpy()])
            error = max(error, float(abs(np.log(actual/expected)).max()))
            accounting = max(accounting, sim.audit()['volume_conversion_amount_error'])
            if torch.any(sim.amount_source) or torch.any(sim.last_amount_source):
                raise ValueError('Dilution-retaining control creates a volume amount source')
        if error > 1e-11 or accounting > 2e-14:
            raise ValueError('Actual-context independent chemistry/accounting gate failed')
        checkpoint = folder/f'{arm}-restart.npz'; sim.checkpoint(checkpoint)
        restart = InitiationControlSimulation.restore(checkpoint)
        contact.exact_state(sim, restart, True)
        for _ in range(2):
            sim.step(); restart.step(); contact.exact_state(sim, restart, True)
        rows.append(dict(arm=arm, passed=True, cpu_chemical_log_max=error,
                         accounting_max=accounting, restart_exact_steps=2,
                         checkpoint=str(checkpoint.resolve()), checkpoint_sha256=digest(checkpoint)))
        del sim, restart
    torch.cuda.empty_cache()
    return dict(passed=True, baseline_exact_steps=16, controls=rows)


def prepare(root, phase=Path('outputs/phase-carry-polarity'), contacts=Path('outputs/moving-initiation-controls')):
    root, phase, contacts = [Path(x).resolve() for x in (root, phase, contacts)]
    if root.exists():
        raise FileExistsError(root)
    old, cp = read(phase/'protocol.json'), read(contacts/'protocol.json')
    verify(old); verify(cp)
    reviews = [Path('docs/phase_carry_polarity_assessment.json'), Path('docs/moving_initiation_controls_assessment.json')]
    inputs = [Path('docs/polarity_conductance_controls.md')]
    for parent, protocol, review in zip((phase,contacts), (old,cp), reviews):
        summary, evidence = read(parent/'summary.json'), read(review)
        if (read(parent/'status.json')['state'] != 'completed' or not summary['passed'] or
                summary['protocol_sha256'] != digest(parent/'protocol.json') or
                evidence['protocol_sha256'] != digest(parent/'protocol.json') or
                evidence['original_summary_sha256'] != digest(parent/'summary.json') or
                digest(evidence['full_review']) != evidence['full_review_sha256']):
            raise ValueError('Requires both completed and reviewed parent studies')
        inputs += [parent/name for name in ('protocol.json','summary.json','status.json')]
        inputs += [review, Path(evidence['full_review'])]
    if any(row['conclusion'] != 'fixing_contact_conductances_restores_formation' for row in read(contacts/'summary.json')['by_history']):
        raise ValueError('Parent does not qualify the proposed conductance-rescue follow-up')
    if (torch.cuda.get_device_name(0) != old['device']['name'] or torch.__version__ != old['device']['torch'] or
            torch.version.cuda != old['device']['torch_cuda'] or digest(library()._name) != old['device']['experimental_binary_sha256']):
        raise ValueError('Changed qualified GPU/runtime/kernel')
    root.mkdir(parents=True)
    sources, reused, gates = {}, {}, []
    for seed in HISTORIES:
        folder = root/f'seed-{seed}'; folder.mkdir(); sources[seed] = {}
        for level, name in (('coarse','coarse-source.npz'),('fine','fine-source.npz'),('chemical','initial_states.npz')):
            original, copied = phase/f'seed-{seed}'/name, folder/name
            shutil.copy2(original, copied); inputs += [original,copied]; sources[seed][level] = str(copied)
        sources[seed]['transport'] = str(folder/'initial-transport.npz')
    jobs = job_design(sources)
    for job in jobs:
        chi = job['point']['chi']
        if job['arm'] == 'baseline':
            original = next(j for j in old['jobs'] if j['seed']==job['seed'] and j['level']==job['level'] and j['point']['chi']==chi)
            item = old['reused'].get(original['key'], dict(root=str(phase),job=original))
            reused[job['key']] = dict(**item,kind='phase',protocol_sha256=digest(Path(item['root'])/'protocol.json'))
        elif chi == 0. and job['seed'] in (7,8):
            original = next(j for j in cp['jobs'] if j['seed']==job['seed'] and j['level']==job['level'] and j['arm']==job['arm'])
            reused[job['key']] = dict(root=str(contacts),job=original,kind='contact',protocol_sha256=digest(contacts/'protocol.json'))
    assert len(jobs)==24 and len(reused)==16
    p = dict(jobs=jobs,reused=reused,accepted_configs=old['accepted_configs'],device=old['device'],
             histories=list(HISTORIES),contrasts=list(CONTRASTS),arms=list(ARMS),
             start=150.,duration=240.,pilot_duration=60.,late_window=24.,pilot_late_window=12.,
             interval=.15,checkpoint_interval=3.,independent_histories=3,new_histories=0,
             new_moving_jobs=8,reused_trajectories=16,endpoint_workers=2,
             criteria=old['criteria'],refinement_criteria=old['refinement_criteria'],
             endpoint_horizons=old['endpoint_horizons'],endpoint_interval=old['endpoint_interval'],endpoint_criteria=old['endpoint_criteria'],
             chemical_reference_log_max=1e-11,amount_accounting_max=2e-14)
    for seed in HISTORIES:
        matrices=[]
        for chi in CONTRASTS:
            for level in ('coarse','fine'):
                job=next(j for j in jobs if j['seed']==seed and j['level']==level and j['point']['chi']==chi and j['arm']=='baseline')
                item=reused[job['key']];original_p=read(Path(item['root'])/'protocol.json')
                r,h=current_evidence(item['root'],original_p,item['job'])
                sim=InitiationControlSimulation(initialize(job['source'],job['chemical_file'],'uniform',job['point'],p['accepted_configs'][level]))
                _,delta,_,g=sim.matrices();matrices.append(g.cpu().numpy().copy())
                np.testing.assert_array_equal(delta.cpu().numpy(),h[0]['delta'])
                np.testing.assert_array_equal(sim.geometry[:,0].cpu().numpy(),h[0]['volumes'])
                np.testing.assert_array_equal([sim.activator.cpu().numpy(),sim.inhibitor.cpu().numpy()],h[0]['chemistry'])
                np.testing.assert_array_equal(sim.polarity.cpu().numpy(),h[0]['polarity'])
                if chi==0. and level=='coarse':
                    np.savez_compressed(job['initial_transport'],conductance=matrices[-1],volumes=h[0]['volumes'],delta=h[0]['delta'],ids=sim.ids)
                    inputs.append(Path(job['initial_transport']))
                del sim
                if chi==.35 or seed==9:
                    gate_folder=root/f'seed-{seed}'/f'gate-{level}-chi-{chi:g}'
                    gate=actual_context_gate(gate_folder,job['source'],job['chemical_file'],job['point'],p['accepted_configs'][level],job['initial_transport'])
                    gates.append(dict(seed=seed,chi=chi,level=level,**gate))
                    inputs += [Path(r['checkpoint']) for r in gate['controls']]
                print(f'initial context verified: history={seed}, chi={chi:g}, level={level}',flush=True)
        for matrix in matrices[1:]:np.testing.assert_array_equal(matrix,matrices[0])
    # Verify every reused intervention and endpoint, not only baseline starts.
    for job in jobs:
        if job['key'] not in reused:continue
        r,h=reused_evidence(p,job)
        item=reused[job['key']];src=Path(item['root'])/item['job']['key']
        ep,e=read(src/'endpoint/protocol.json'),read(src/'endpoint/assay/result.json');verify(ep)
        if (not e['numerical_pass'] or not e['all_trials_settled'] or e['protocol_sha256']!=digest(src/'endpoint/protocol.json') or
                e['paths_sha256']!=digest(src/'endpoint/assay/paths.npz') or ep['moving_protocol_sha256']!=item['protocol_sha256'] or
                ep['moving_history_sha256']!=digest(src/'history.json')):
            raise ValueError('Unqualified reused endpoint')
        inputs += [Path(item['root'])/'protocol.json']+[src/name for name in ('result.json','history.json','latest_state.npz','endpoint/protocol.json','endpoint/source.npz','endpoint/assay/result.json','endpoint/assay/paths.npz')]
    gate_file=root/'implementation-verification.json';write_json(gate_file,dict(passed=True,new_contexts=gates,reused_zero_contrast_gate=str(contacts/'implementation-verification.json')))
    inputs += [gate_file,contacts/'implementation-verification.json']
    source_files=[Path(__file__),Path(__file__).parents[1]/'tests/test_polarity_conductance_controls.py']
    p['source_sha256']={**old['source_sha256'],**cp['source_sha256'],**{str(f.resolve()):digest(f) for f in source_files}}
    p['input_sha256']={**old['input_sha256'],**cp['input_sha256'],**{str(f.resolve()):digest(f) for f in inputs}}
    medians={level:float(np.median([read(contacts/j['key']/'result.json')['audit']['wall_seconds'] for j in cp['jobs'] if j['level']==level and j['key'] not in cp['reused']])) for level in ('coarse','fine')}
    p['estimated_new_moving_seconds']=sum(medians[j['level']] for j in jobs if j['key'] not in reused)
    p['design']='Three existing histories 7/8/9; chi=0/.35 crossed with evolving/fixed initial symmetric chemical conductances. Both timesteps and mature near-uniform t=150 starts, 240 elapsed units. Twenty-four numerical paths: eight new, sixteen qualified read-only reuses. History 9 scheduled first; no new developmental histories.'
    p['interpretation']='Require full raw refinement, physical and endpoint checks. Only history 9 currently demonstrates original polarity-specific suppression; histories 7/8 have neither original branch forming. Fixed-conductance rescue at positive contrast identifies a tested conductance contribution, not universal or complete causal mediation. A forming fixed chi=0 control is necessary to interpret residual fixed-conductance polarity suppression. Dilution retained, actual M(t) retained: Delta=-M(t)^(-1)K0. No new biology/identity/inheritance/zygote/spatial/closure/backend qualification.'
    write_json(root/'protocol.json',p)
    # Freeze reused paired evidence before launch; later reads verify it exactly.
    for seed in HISTORIES:
        for chi in CONTRASTS:
            for arm in ARMS:
                group=[j for j in jobs if j['seed']==seed and j['point']['chi']==chi and j['arm']==arm]
                if all(j['key'] in reused for j in group):
                    h={j['level']:reused_evidence(p,j)[1] for j in group}
                    for horizon in (p['pilot_duration'],p['duration']):
                        cut={level:[r for r in history if r['elapsed']<=horizon+1e-9] for level,history in h.items()}
                        if not pair_report(root,p,group[0],horizon,cut)['passed']:
                            raise ValueError('Reused paired timestep evidence fails the original gate')
    write_json(root/'status.json',dict(state='prepared',completed=0,total=8,reused_trajectories=16,estimated_new_moving_seconds=p['estimated_new_moving_seconds']))
    return p


def pair_report(root,p,job,horizon,histories):
    return contact.pair_report(root,p,job['seed'],context_key(job),horizon,histories)


def advance(root,job,p,horizon):
    if horizon not in (p['pilot_duration'],p['duration']):raise ValueError('Undeclared horizon')
    if job['key'] in p['reused']:
        _,h=reused_evidence(p,job)
        return validate_history([r for r in h if r['elapsed']<=horizon+1e-9],p,job,horizon)
    root=Path(root);folder=root/job['key'];folder.mkdir(exist_ok=True);cp=folder/'latest_state.npz';ph=digest(root/'protocol.json')
    if (folder/'result.json').exists():
        _,h=completed_evidence(root,p,job)
        return validate_history([r for r in h if r['elapsed']<=horizon+1e-9],p,job,horizon)
    if cp.exists():
        from dataclasses import asdict
        host,audit,history=restore_checkpoint(cp,ph,job);check_config(asdict(host.config),job,p)
        sim=InitiationControlSimulation.restore(cp)
        if sim.initiation_arm!=job['arm']:raise ValueError('Restart intervention changed')
        validate_history(history,p,job,sim.time-p['start'])
        np.testing.assert_array_equal(history[-1]['chemistry'],[host.activator,host.inhibitor])
    else:
        sim=InitiationControlSimulation(initialize(job['source'],job['chemical_file'],'uniform',job['point'],p['accepted_configs'][job['level']]),job['arm'])
        history=[];audit=dict(max_volume_error=0.,min_radius=1e100,max_clipping=0.,boundary_max=0.,dilution_error_max=0.,wall_seconds=0.)
    with np.load(job['initial_transport']) as z:
        np.testing.assert_array_equal(sim.initial_conductance.cpu().numpy(),z['conductance'])
    start,stop=_steps(p['start'],job['dt']),_steps(p['start']+horizon,job['dt'])
    every,save_every=_steps(p['interval'],job['dt']),_steps(p['checkpoint_interval'],job['dt'])
    if not start<=sim.step_number<=_steps(p['start']+p['duration'],job['dt']) or abs(sim.time-sim.step_number*job['dt'])>1e-9:
        raise ValueError('Invalid continuation clock')
    if sim.step_number>stop:
        del sim;torch.cuda.empty_cache()
        return validate_history([r for r in history if r['elapsed']<=horizon+1e-9],p,job,horizon)
    clock=time.perf_counter();base=audit['wall_seconds']
    while sim.step_number<=stop:
        quality=sim.audit()
        for key in ('max_volume_error','max_clipping'):audit[key]=max(audit[key],quality[key])
        audit['min_radius']=min(audit['min_radius'],quality['min_radius'])
        audit['dilution_error_max']=max(audit['dilution_error_max'],quality['volume_conversion_amount_error'])
        if audit['dilution_error_max']>p['amount_accounting_max']:raise RuntimeError('Amount accounting failure')
        offset=sim.step_number-start
        if offset%every==0:
            elapsed=offset*job['dt']
            if not history or abs(history[-1]['elapsed']-elapsed)>1e-10:
                row=retain(sim.observe(elapsed),job);row['precision']=sim.precision_diagnostics();history.append(row)
            audit['boundary_max']=max(audit['boundary_max'],history[-1]['boundary_occupancy'])
            if audit['boundary_max']>=p['criteria']['boundary_max']:raise RuntimeError('Boundary failure')
            audit['wall_seconds']=base+time.perf_counter()-clock
            write_json(folder/'status.json',dict(state='running',elapsed=elapsed,target_horizon=horizon,audit=audit,log_activator_sd=history[-1]['log_activator_sd']))
            if offset%save_every==0 or sim.step_number==stop:
                save_checkpoint(sim,cp,audit,history,ph,job);write_json(folder/'history.json',history)
                print(f'{job["key"]}: elapsed={elapsed:g}, target={horizon:g}, wall={audit["wall_seconds"]:.1f}s',flush=True)
        if sim.step_number==stop:break
        sim.step()
    validate_history(history,p,job,horizon)
    if horizon==p['duration']:
        write_json(folder/'result.json',dict(job=job,protocol_sha256=ph,history_sha256=digest(folder/'history.json'),checkpoint_sha256=digest(cp),quality_pass=True,audit=audit,**trajectory_summary(history,p)))
        write_json(folder/'status.json',dict(state='completed',elapsed=horizon))
    else:write_json(folder/'status.json',dict(state='awaiting_pair_refinement',elapsed=horizon,audit=audit))
    del sim;torch.cuda.empty_cache();return history


def run_pairs(root,p,executor):
    done,pending,failed=0,[],[]
    for seed in p['histories']:
        for chi in p['contrasts']:
            for arm in p['arms']:
                jobs=[j for j in p['jobs'] if j['seed']==seed and j['point']['chi']==chi and j['arm']==arm]
                h={}
                for job in jobs:
                    write_json(Path(root)/'status.json',dict(state='running',stage='matched_pilot',current_job=job['key'],completed=done,total=p['new_moving_jobs']))
                    h[job['level']]=advance(root,job,p,p['pilot_duration'])
                if not pair_report(root,p,jobs[0],p['pilot_duration'],h)['passed']:
                    failed.append(dict(seed=seed,chi=chi,arm=arm));continue
                for job in jobs:
                    write_json(Path(root)/'status.json',dict(state='running',stage='long_moving_controls',current_job=job['key'],completed=done,total=p['new_moving_jobs'],failed_pilot_contexts=failed))
                    h[job['level']]=advance(root,job,p,p['duration'])
                    if job['key'] not in p['reused']:
                        done+=1;pending.append(executor.submit(endpoint_assay,root,job,p,h[job['level']][-1]))
                pair_report(root,p,jobs[0],p['duration'],h)
    return done,pending,failed


def assess(root,p):
    verify(p);rows,endpoints,missing=[],[],[]
    for job in p['jobs']:
        if job['key'] in p['reused']:
            _,h=reused_evidence(p,job);item=p['reused'][job['key']]
            src=Path(item['root'])/item['job']['key']
            e=read(src/'endpoint/assay/result.json');ep=read(src/'endpoint/protocol.json');verify(ep)
            if e['protocol_sha256']!=digest(src/'endpoint/protocol.json') or e['paths_sha256']!=digest(src/'endpoint/assay/paths.npz'):
                raise ValueError('Changed reused endpoint')
        elif not (Path(root)/job['key']/'result.json').exists():
            missing.append(job['key']);continue
        else:
            _,h=completed_evidence(root,p,job);e=endpoint_assay(root,job,p,h[-1])
        rows.append(dict(job=job,reused=job['key'] in p['reused'],**trajectory_summary(h,p)))
        endpoints.append(dict(job=job['key'],numerical_pass=e['numerical_pass'],all_trials_settled=e['all_trials_settled'],phase=e['phase'],local_bistability=e['local_bistability_supported']))
    by_history=[]
    for seed in p['histories']:
        selected=[e for e in endpoints if next(j for j in p['jobs'] if j['key']==e['job'])['seed']==seed]
        qualified=len(selected)==8 and all(e['numerical_pass'] and e['all_trials_settled'] for e in selected)
        forms={}
        for chi in p['contrasts']:
            for arm in p['arms']:
                job=next(j for j in p['jobs'] if j['seed']==seed and j['point']['chi']==chi and j['arm']==arm)
                file=Path(root)/'pairs'/f'seed-{seed}_{context_key(job)}'/'full-refinement.json'
                if not file.exists():qualified=False;continue
                histories={level:read(file.parent/f'full-{level}.json') for level in ('coarse','fine')}
                qualified &= pair_report(root,p,job,p['duration'],histories)['passed']
                matching=[r for r in rows if r['job']['seed']==seed and r['job']['point']['chi']==chi and r['job']['arm']==arm]
                if len(matching)!=2:qualified=False;continue
                forms[context_key(job)]=next(r['persistent_contrast'] for r in matching if r['job']['level']=='fine')
        by_history.append(dict(seed=seed,formation=forms,**interpretation(forms,qualified)))
    passed=not missing and all(r['qualified'] for r in by_history)
    result=dict(protocol_sha256=digest(Path(root)/'protocol.json'),passed=bool(passed),completed=sum(not r['reused'] for r in rows),total=8,
                independent_histories=3,new_histories=0,reused_trajectories=16,results=rows,missing=missing,by_history=by_history,endpoints=endpoints,scope=p['interpretation'])
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
                write_json(root/'status.json',dict(state='running',stage='endpoint_assays',completed=done,total=8,failed_pilot_contexts=failed))
                for future in pending:future.result()
            result=assess(root,p)
            write_json(root/'status.json',dict(state='completed' if result['passed'] else 'completed_with_unresolved_checks',completed=result['completed'],total=8,passed=result['passed'],failed_pilot_contexts=failed))
            return result
        except Exception as error:
            write_json(root/'status.json',dict(state='failed',error=str(error)));raise


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=('prepare','run','assess'))
    parser.add_argument('--output',type=Path,default=Path('outputs/polarity-conductance-controls'))
    args=parser.parse_args()
    result=prepare(args.output) if args.action=='prepare' else run(args.output) if args.action=='run' else assess(args.output,read(args.output/'protocol.json'))
    print({k:v for k,v in result.items() if k in ('passed','completed','total','estimated_new_moving_seconds','by_history')})
