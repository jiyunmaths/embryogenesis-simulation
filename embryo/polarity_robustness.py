"""Fixed-ratio polarity test with long matched horizons and targeted time refinement.

Keep previous protocols and backend code immutable. Gate long continuations on
same-state 60-unit coarse/fine agreement in the prespecified initiation history.
"""
import argparse
from dataclasses import asdict
import fcntl
from pathlib import Path
import shutil
import time

import numpy as np

from .attribute_development import AttributeSimulation
from .attribute_persistence import rhs
from .cell_response_exchange import checked_solve
from .feedback_endpoint_bistability import chemical_jacobian
from .feedback_long import digest, save_checkpoint, restore_checkpoint
from .feedback_survival_validation import retime
from .gpu_response_runner import require_validation, compatible, PARAMETERS
from .neighbor_context import read, validate_graph
from .parameter_robustness import verify, CRITERIA as ENDPOINT_CRITERIA
from .parameter_robustness_moving import initialize, prefix, retain, endpoint_assay
from .resolution import write_json, _steps
from .validate_gpu_backend import CRITERIA as GPU_CRITERIA, discrepancies


CONTRASTS = (0., .35, .7)
REFINEMENT_CRITERIA = dict(chemical_log_max=.01, polarity_abs_max=.01,
    relative_axis_max=.01, relative_volume_max=.005, relative_transport_max=.01,
    growth_abs_max=.001, crossing_time_abs_max=.30, onset_time_abs_max=.30,
    same_late_contrast_required=True, same_crossing_presence_required=True)


def job_view(p, job):
    """Per-job timestep/config is declared in the master protocol, not inferred."""
    view = dict(p, dt=job['dt'], accepted_config=p['accepted_configs'][job['level']])
    return view


def verify_completed(parent, p):
    if read(parent/'status.json').get('state') != 'completed' or len(p['jobs']) != 18:
        raise ValueError('Requires the completed three-point parent study')
    files = []
    ph = digest(parent/'protocol.json')
    for job in p['jobs']:
        folder = parent/job['key']; r = read(folder/'result.json')
        e = read(folder/'endpoint/assay/result.json'); ep = read(folder/'endpoint/protocol.json'); verify(ep)
        cp = read(folder/'prefix/comparison.json')
        if (not r['quality_pass'] or not e['numerical_pass'] or not cp['passed'] or
                r['protocol_sha256'] != ph or cp['protocol_sha256'] != ph or
                r['checkpoint_sha256'] != digest(folder/'latest_state.npz') or
                r['history_sha256'] != digest(folder/'history.json') or
                e['protocol_sha256'] != digest(folder/'endpoint/protocol.json') or
                e['paths_sha256'] != digest(folder/'endpoint/assay/paths.npz')):
            raise ValueError('Unaccepted parent evidence: '+job['key'])
        for f,h in cp['evidence_sha256'].items():
            if digest(f) != h: raise ValueError('Changed parent context gate')
        files.extend(folder/f for f in ('result.json', 'history.json', 'latest_state.npz',
            'prefix/comparison.json', 'endpoint/protocol.json', 'endpoint/assay/result.json', 'endpoint/assay/paths.npz'))
    return files


def assert_same_physical_start(coarse, fine):
    if coarse.time != fine.time or coarse.ids.tolist() != fine.ids.tolist():
        raise ValueError('Refinement changed the starting clock/cell order')
    for name in ('phi', 'polarity', 'target', 'activator', 'inhibitor', 'parents', 'due'):
        if not np.array_equal(getattr(coarse, name), getattr(fine, name)):
            raise ValueError('Refinement changed physical state: '+name)
    for name in ('rng', 'signal_rng', 'fate_rng'):
        if getattr(coarse, name).bit_generator.state != getattr(fine, name).bit_generator.state:
            raise ValueError('Refinement changed a random stream')


def prepare(root, parent=Path('outputs/parameter-robustness-moving'),
            fine_validation=Path('outputs/exchange-response-histories-refined/gpu-validation')):
    root, parent, fine_validation = [Path(x).resolve() for x in (root, parent, fine_validation)]
    if root.exists(): raise FileExistsError(root)
    old = read(parent/'protocol.json'); verify(old)
    evidence = verify_completed(parent, old)
    accepted_coarse, coarse_evidence = require_validation(old['gpu_validation'])
    accepted_fine, fine_evidence = require_validation(fine_validation)
    if accepted_fine.dt != accepted_coarse.dt/2:
        raise ValueError('Expected exactly one timestep halving')
    fine_check = type(accepted_fine)(**asdict(accepted_fine)); fine_check.dt = accepted_coarse.dt
    compatible(fine_check, accepted_coarse)
    root.mkdir(parents=True)
    inputs = [parent/'protocol.json', parent/'status.json', *evidence, *coarse_evidence, *fine_evidence]
    jobs, reused = [], []
    for seed in (7, 8, 9):
        original = parent/f'seed-{seed}'
        folder = root/f'seed-{seed}'; folder.mkdir()
        coarse_source, fine_source, chemical = [folder/f for f in ('coarse-source.npz', 'fine-source.npz', 'initial_states.npz')]
        shutil.copy2(original/'source.npz', coarse_source)
        shutil.copy2(original/'initial_states.npz', chemical)
        host = AttributeSimulation.restore(coarse_source); compatible(host.config, accepted_coarse)
        refined = retime(coarse_source, fine_source, accepted_fine.dt); compatible(refined.config, accepted_fine)
        assert_same_physical_start(host, refined)
        inputs.extend((original/'source.npz', original/'initial_states.npz', coarse_source, fine_source, chemical))
        for chi in CONTRASTS:
            point = dict(key=f'polarity-{chi:g}', ratio=27.5, chi=chi)
            for family in ('uniform', 'pattern'):
                job = dict(key=f'seed-{seed}_coarse_{point["key"]}_{family}', seed=seed, level='coarse',
                           dt=accepted_coarse.dt, point=point, family=family,
                           source=str(coarse_source), chemical_file=str(chemical))
                jobs.append(job)
                if chi == .7:
                    old_job = next(j for j in old['jobs'] if j['seed'] == seed and j['family'] == family and
                                   j['point']['ratio'] == 27.5 and j['point']['chi'] == chi)
                    reused.append(dict(job_key=job['key'], parent_job=old_job, parent_folder=str(parent/old_job['key'])))
            if seed == 9:
                jobs.append(dict(key=f'seed-{seed}_fine_{point["key"]}_uniform', seed=seed, level='fine',
                    dt=accepted_fine.dt, point=point, family='uniform', source=str(fine_source), chemical_file=str(chemical)))
    sources = {**old['source_sha256'], str(Path(__file__).resolve()):digest(__file__)}
    p = dict(parent=str(parent), parent_protocol_sha256=digest(parent/'protocol.json'),
        validation_roots=dict(coarse=old['gpu_validation'], fine=str(fine_validation)),
        accepted_configs=dict(coarse=asdict(accepted_coarse), fine=asdict(accepted_fine)),
        jobs=jobs, reused=reused, independent_histories=3, total_jobs=len(jobs), ratio=27.5,
        contrasts=list(CONTRASTS), start=150., duration=240., interval=.15, checkpoint_interval=3.,
        late_window=24., pilot_duration=60., pilot_late_window=12., refinement_seed=9,
        prefix_duration=.6, prefix_native_threads=4, prefix_criteria=GPU_CRITERIA,
        criteria=dict(boundary_max=.01, dilution_error_max=GPU_CRITERIA['dilution_amount_error'], late_log_sd_min=.1),
        refinement_criteria=REFINEMENT_CRITERIA, endpoint_criteria=ENDPOINT_CRITERIA,
        endpoint_horizons=[240., 960.], endpoint_interval=2.,
        stages=['six unique matched frozen references to 240',
                'three coarse/fine near-uniform history-9 pilots to 60',
                'require the three matched 60-unit comparisons to pass before long continuation',
                'continue all eighteen coarse and three history-9 near-uniform fine jobs to 240',
                'endpoint chemical assays and full-horizon history-9 refinement assessment'],
        design='Fix beta=2, D_a=.02, D_b=.55 and activity-tension/adhesion laws. Vary only directional polarity contrast across 0,.35,.7. Use identical t=150 geometry/polarity and the exact paired chemical starts from the prior study within each history. No new noise, fate, preconditioning or zygote histories. Long horizon 240; coarse dt=.00375; near-uniform seed-9 checks at dt=.001875.',
        refinement_selection='History 9 is prespecified because it had the largest initial uniform growth and the strongest frozen near-uniform contrast in the completed parent; all three polarity settings are checked, not selected by new outcomes. Fine near-uniform runs continue through the complete 240-unit window. Patterned starts are not timestep-refined in this targeted study.',
        transfer='Reuse all six unchanged coarse chi=.7 first-60-unit trajectories and physical checkpoints. Validate parent protocol, complete observations, state, geometry, chemistry and material parameters before re-encoding checkpoint metadata for the new declared job. Never reset chemistry or polarity at t=210. Exact-context native/GPU gates are rechecked for each new job.',
        interpretation='Late contrast above .1 over the last 24 units, instantaneous spectra, chemical onset and local frozen endpoint stability are distinct outcomes. Missing contrast does not establish absence of other attractors. Fine/coarse disagreement is a numerical limitation, not an alternate mechanism.',
        limits='Three shared histories; no fresh zygote initiation, spatial convergence, GPU cleavage, arbitrary-parameter full-horizon native/GPU equivalence, or autonomous/irreversible identity. Targeted time refinement covers near-uniform seed 9 only. Area/conductance proxies retain the existing geometric closure.',
        source_sha256=sources, input_sha256={str(f.resolve()):digest(f) for f in inputs})
    write_json(root/'protocol.json', p)
    install_reused(root, p)
    write_json(root/'status.json', dict(state='prepared', stage='frozen_and_timestep_gate', completed=0, total=len(jobs),
                                      independent_histories=3, reused_prefixes=6))
    return p


def install_reused(root, p):
    """Transfer verified physical continuations, never reinitialize old chemistry."""
    root = Path(root); ph = digest(root/'protocol.json')
    for reuse in p['reused']:
        job = next(j for j in p['jobs'] if j['key'] == reuse['job_key'])
        folder = root/job['key']; folder.mkdir(exist_ok=True)
        if (folder/'latest_state.npz').exists(): continue
        parent = Path(reuse['parent_folder']); r = read(parent/'result.json')
        host, audit, history = restore_checkpoint(parent/'latest_state.npz', p['parent_protocol_sha256'], reuse['parent_job'])
        expected = initialize(job['source'], job['chemical_file'], job['family'], job['point'], p['accepted_configs']['coarse'])
        if (abs(host.time-(p['start']+p['pilot_duration'])) > 1e-9 or not np.array_equal(host.ids, expected.ids) or
                any(getattr(host.config,k)!=getattr(expected.config,k) for k in PARAMETERS) or
                history != read(parent/'history.json') or audit != r['audit']):
            raise ValueError('Invalid reuse of physical prefix')
        parent_start = AttributeSimulation.restore(reuse['parent_job']['source'])
        for name in ('phi', 'polarity', 'target'):
            # Compare the recorded parent start to the new start, not to its evolved endpoint.
            if not np.array_equal(getattr(parent_start, name), getattr(expected, name)):
                raise ValueError('Changed geometry/polarity before transfer')
        with np.load(job['chemical_file']) as z:
            if not np.array_equal(np.array(history[0]['chemistry']), z[job['family']]):
                raise ValueError('Changed chemical start before transfer')
        checked_history(history, p, job, p['pilot_duration'])
        if not np.array_equal(np.array([host.activator,host.inhibitor]),history[-1]['chemistry']):
            raise ValueError('Transferred chemical endpoint differs from observations')
        save_checkpoint(host, folder/'latest_state.npz', audit, history, ph, job)
        write_json(folder/'history.json', history)
        write_json(folder/'transfer.json', dict(parent_checkpoint=str(parent/'latest_state.npz'),
            parent_checkpoint_sha256=digest(parent/'latest_state.npz'), new_protocol_sha256=ph,
            preserved_chemistry=True, preserved_geometry_polarity=True, elapsed=p['pilot_duration']))
        write_json(folder/'status.json', dict(state='awaiting_timestep_gate', elapsed=p['pilot_duration'], audit=audit))


def checked_history(history, p, job, duration):
    expected = np.arange(_steps(duration, p['interval'])+1)*p['interval']
    if len(history) != len(expected) or not np.allclose([r['elapsed'] for r in history], expected, rtol=0, atol=1e-9):
        raise ValueError('Incomplete or misaligned observation clocks')
    with np.load(job['chemical_file']) as z:
        if any(not np.array_equal(r['ids'], z['ids']) for r in history) or not np.array_equal(history[0]['chemistry'], z[job['family']]):
            raise ValueError('Cell order/initial chemistry changed')
    for row in history:
        if abs(row['time']-p['start']-row['elapsed']) > 1e-9:
            raise ValueError('Physical and elapsed clocks disagree')
        for name in ('chemistry','volumes','delta'):
            if not np.isfinite(np.asarray(row[name],float)).all():
                raise ValueError('Nonfinite observation: '+name)
        if np.any(np.asarray(row['chemistry'])<=0) or np.any(np.asarray(row['volumes'])<=0):
            raise ValueError('Nonpositive chemical state/capacity')
    return history


def frozen_references(root, p):
    root = Path(root); folder = root/'frozen'; folder.mkdir(exist_ok=True)
    file = folder/'results.json'; ph = digest(root/'protocol.json')
    if file.exists():
        result = read(file)
        if result['protocol_sha256'] != ph: raise ValueError('Changed frozen reference protocol')
        for f,h in result['paths_sha256'].items():
            if digest(f) != h: raise ValueError('Changed long frozen reference')
        return result
    times = np.arange(_steps(p['duration'], p['interval'])+1)*p['interval']
    rows, paths = [], {}
    for seed in (7,8,9):
        host = AttributeSimulation.restore(root/f'seed-{seed}/coarse-source.npz'); graph = host.signaling_graph()
        delta, masses = validate_graph(graph.delta, graph.masses)
        fun = rhs(delta, 2., .02, .55)
        jac = lambda t,y:chemical_jacobian(y.reshape(2,-1), delta, 2., .02, .55)
        with np.load(root/f'seed-{seed}/initial_states.npz') as z:
            for family in ('uniform','pattern'):
                initial = z[family].copy(); trajectory, error = checked_solve(fun,jac,initial,times)
                spread = np.std(np.log(trajectory[:,0]),axis=1)
                target = folder/f'seed-{seed}_{family}.npz'
                np.savez_compressed(target,times=times,trajectory=trajectory,delta=delta,masses=masses,ids=z['ids'])
                paths[str(target.resolve())]=digest(target)
                late_min=float(spread[times>=p['duration']-p['late_window']-1e-9].min())
                rows.append(dict(seed=seed,family=family,solver_log_error=error,final_log_sd=float(spread[-1]),
                    late_min_log_sd=late_min,persistent_contrast=bool(late_min>p['criteria']['late_log_sd_min']),
                    first_contrast_crossing=first_crossing(times,spread,p['criteria']['late_log_sd_min'],'up')))
    result=dict(protocol_sha256=ph,independent_histories=3,unique_trajectories=6,results=rows,paths_sha256=paths,
                interpretation='One frozen ODE reference per history/start. Polarity and mechanics timestep are absent; share the same reference across those interventions.')
    write_json(file,result);return result


def first_crossing(times, values, threshold=0., direction='down'):
    times, values = np.asarray(times,float), np.asarray(values,float)
    if direction not in ('up','down') or times.ndim!=1 or values.shape!=times.shape or len(times)<2 or not np.isfinite(times).all() or not np.isfinite(values).all() or np.any(np.diff(times)<=0):
        raise ValueError('Finite ordered matching observations and up/down direction required')
    mask=(values[:-1]<=threshold)&(values[1:]>threshold) if direction=='up' else (values[:-1]>threshold)&(values[1:]<=threshold)
    indices=np.flatnonzero(mask)
    if not len(indices):return None
    i=int(indices[0]);return float(times[i]+(threshold-values[i])*(times[i+1]-times[i])/(values[i+1]-values[i]))


def refinement_comparison(coarse, fine, p, horizon):
    """No floor-normalized relative contrast at near-uniform states."""
    times=np.arange(_steps(horizon,p['interval'])+1)*p['interval']
    if len(coarse)!=len(fine) or len(coarse)!=len(times) or not np.allclose([r['elapsed'] for r in coarse],times,rtol=0,atol=1e-9) or not np.allclose([r['elapsed'] for r in fine],times,rtol=0,atol=1e-9):
        raise ValueError('Coarse/fine times are not aligned')
    if not np.array_equal(coarse[0]['chemistry'],fine[0]['chemistry']) or not np.array_equal(coarse[0]['polarity'],fine[0]['polarity']):
        raise ValueError('Timestep comparison changed chemical/polarity starts')
    errors=dict.fromkeys(('chemical_log_max','polarity_abs_max','relative_axis_max','relative_volume_max','relative_transport_max'),0.)
    for a,b in zip(coarse,fine):
        for key,value in discrepancies(a,b).items():errors[key]=max(errors[key],value)
    cg=np.array([r['uniform_growth_max'] for r in coarse]);fg=np.array([r['uniform_growth_max'] for r in fine])
    cs=np.array([r['log_activator_sd'] for r in coarse]);fs=np.array([r['log_activator_sd'] for r in fine])
    errors['growth_abs_max']=float(abs(cg-fg).max())
    crit=p['refinement_criteria']
    def crossing_error(a,b):return 0. if a is None and b is None else None if a is None or b is None else abs(a-b)
    crossings=[first_crossing(times,x) for x in (cg,fg)]
    onsets=[first_crossing(times,x,p['criteria']['late_log_sd_min'],'up') for x in (cs,fs)]
    ce=crossing_error(*crossings);oe=crossing_error(*onsets)
    late=p['pilot_late_window'] if horizon==p['pilot_duration'] else p['late_window']
    outcome=[bool(x[times>=horizon-late-1e-9].min()>p['criteria']['late_log_sd_min']) for x in (cs,fs)]
    passed=all(v<=crit[k] for k,v in errors.items()) and ce is not None and ce<=crit['crossing_time_abs_max'] and oe is not None and oe<=crit['onset_time_abs_max'] and outcome[0]==outcome[1]
    return dict(passed=bool(passed),errors=errors,crossing_times=crossings,crossing_time_error=ce,
        onset_times=onsets,onset_time_error=oe,persistent_contrast=outcome,horizon=horizon,
        interpretation='Same-state timestep halving in prespecified near-uniform history 9. Instantaneous modal growth is a frozen-snapshot diagnostic, not full moving-system stability.')


def advance(root, job, p, horizon):
    """Checkpoint at the declared pilot stop without publishing a final result."""
    from .gpu_backend import GpuSimulation
    root=Path(root);folder=root/job['key'];folder.mkdir(exist_ok=True)
    view=job_view(p,job);ph=digest(root/'protocol.json')
    if horizon not in (p['pilot_duration'],p['duration']):raise ValueError('Only declared stage horizons are allowed')
    prefix(root,job,view)
    cp=folder/'latest_state.npz'
    if (folder/'result.json').exists():
        result=read(folder/'result.json');history=read(folder/'history.json')
        if (result['protocol_sha256']!=ph or not result['quality_pass'] or
                result['history_sha256']!=digest(folder/'history.json') or result['checkpoint_sha256']!=digest(cp)):
            raise ValueError('Changed completed long evidence')
        checked_history(history,p,job,p['duration'])
        endpoint_assay(root,job,view,history[-1])
        return checked_history([r for r in history if r['elapsed']<=horizon+1e-9],p,job,horizon)
    if cp.exists():
        host,audit,history=restore_checkpoint(cp,ph,job)
        sim=GpuSimulation(host)
    else:
        sim=GpuSimulation(initialize(job['source'],job['chemical_file'],job['family'],job['point'],view['accepted_config']))
        history=[];audit=dict(max_volume_error=0.,min_radius=1e100,max_clipping=0.,boundary_max=0.,dilution_error_max=0.,wall_seconds=0.)
    start=_steps(p['start'],job['dt']);stop=_steps(p['start']+horizon,job['dt']);every=_steps(p['interval'],job['dt']);checkpoint_every=_steps(p['checkpoint_interval'],job['dt'])
    if sim.step_number<start or sim.step_number>_steps(p['start']+p['duration'],job['dt']) or abs(sim.time-sim.step_number*job['dt'])>1e-9:
        raise ValueError('Invalid continuation clock')
    if history:
        elapsed=(sim.step_number-start)*job['dt']
        checked_history(history,p,job,elapsed)
        if not np.array_equal(history[-1]['chemistry'],np.array([host.activator,host.inhibitor])):
            raise ValueError('Checkpoint chemistry differs from observations')
    if sim.step_number>stop:
        return checked_history([r for r in history if r['elapsed']<=horizon+1e-9],p,job,horizon)
    clock=time.perf_counter();base_wall=audit['wall_seconds']
    try:
        while sim.step_number<=stop:
            quality=sim.audit()
            for key in ('max_volume_error','max_clipping'):audit[key]=max(audit[key],quality[key])
            audit['min_radius']=min(audit['min_radius'],quality['min_radius'])
            audit['dilution_error_max']=max(audit['dilution_error_max'],quality['dilution_amount_error'])
            if audit['dilution_error_max']>p['criteria']['dilution_error_max']:raise RuntimeError('Dilution amount failure')
            offset=sim.step_number-start
            if offset%every==0:
                elapsed=offset*job['dt']
                if not history or abs(history[-1]['elapsed']-elapsed)>1e-10:history.append(retain(sim.observe(elapsed),job))
                row=history[-1];audit['boundary_max']=max(audit['boundary_max'],row['boundary_occupancy'])
                if audit['boundary_max']>=p['criteria']['boundary_max']:raise RuntimeError('Boundary screen failed')
                audit['wall_seconds']=base_wall+time.perf_counter()-clock
                write_json(folder/'status.json',dict(state='running',target_horizon=horizon,elapsed=elapsed,audit=audit,
                    log_activator_sd=row['log_activator_sd'],uniform_growth_max=row['uniform_growth_max']))
                if offset%checkpoint_every==0 or sim.step_number==stop:
                    save_checkpoint(sim.to_cpu(),cp,audit,history,ph,job);write_json(folder/'history.json',history)
                    print(f"{job['key']} elapsed={elapsed:.2f}/{horizon:g} SD={row['log_activator_sd']:.5f} growth={row['uniform_growth_max']:.5f}",flush=True)
            if sim.step_number==stop:break
            sim.step()
        checked_history(history,p,job,horizon)
        if horizon==p['duration']:
            late=[r['log_activator_sd'] for r in history if r['elapsed']>=horizon-p['late_window']-1e-9]
            result=dict(job=job,protocol_sha256=ph,quality_pass=True,audit=audit,history_sha256=digest(folder/'history.json'),
                checkpoint_sha256=digest(cp),late_min_log_activator_sd=min(late),persistent_contrast=bool(min(late)>p['criteria']['late_log_sd_min']))
            write_json(folder/'result.json',result);e=endpoint_assay(root,job,view,history[-1])
            write_json(folder/'status.json',dict(state='completed',elapsed=horizon,endpoint_numerical_pass=e['numerical_pass']))
        else:
            write_json(folder/'status.json',dict(state='awaiting_timestep_gate',elapsed=horizon,audit=audit))
        return history
    except Exception as error:
        write_json(folder/'status.json',dict(state='failed',time=sim.time,error=str(error),audit=audit));raise


def compare_level(root,p,horizon):
    root=Path(root);ph=digest(root/'protocol.json')
    pilot=horizon==p['pilot_duration']
    target=root/('pilot-refinement.json' if pilot else 'long-refinement.json')
    if target.exists():
        result=read(target)
        if result['protocol_sha256']!=ph or result['horizon']!=horizon:raise ValueError('Changed refinement protocol')
        for f,h in result['source_histories_sha256'].items():
            if digest(f)!=h:raise ValueError('Changed refinement evidence')
        return result
    evidence={}
    rows=[]
    for chi in p['contrasts']:
        jobs=[next(j for j in p['jobs'] if j['seed']==p['refinement_seed'] and j['family']=='uniform' and j['point']['chi']==chi and j['level']==level) for level in ('coarse','fine')]
        histories=[]
        for job in jobs:
            folder=Path(root)/job['key'];h=read(folder/'history.json')
            part=[r for r in h if r['elapsed']<=horizon+1e-9];checked_history(part,p,job,horizon)
            host,audit,saved=restore_checkpoint(folder/'latest_state.npz',digest(Path(root)/'protocol.json'),job)
            if saved!=h or host.time<p['start']+horizon-1e-9:raise ValueError('Timestep history/checkpoint mismatch')
            # Preserve the pilot observations before later histories are appended.
            source=folder/'pilot-history.json' if pilot else folder/'history.json'
            if pilot:
                if source.exists() and read(source)!=part:raise ValueError('Changed immutable pilot history')
                if not source.exists():write_json(source,part)
            evidence[str(source.resolve())]=digest(source)
            histories.append(part)
        rows.append(dict(chi=chi,**refinement_comparison(*histories,p,horizon)))
    result=dict(passed=all(r['passed'] for r in rows),seed=p['refinement_seed'],horizon=horizon,comparisons=rows,
        protocol_sha256=ph,source_histories_sha256=evidence)
    write_json(target,result);return result


def assess(root):
    root=Path(root);p=read(root/'protocol.json');verify(p);rows=[];failures=[]
    refs=frozen_references(root,p);byref={(r['seed'],r['family']):r for r in refs['results']}
    for job in p['jobs']:
        folder=root/job['key']
        if (folder/'status.json').exists() and read(folder/'status.json')['state']=='failed':failures.append(dict(key=job['key'],error=read(folder/'status.json')['error']))
        if not (folder/'result.json').exists():continue
        result=read(folder/'result.json');h=read(folder/'history.json');checked_history(h,p,job,p['duration'])
        if not result['quality_pass'] or result['protocol_sha256']!=digest(root/'protocol.json') or result['history_sha256']!=digest(folder/'history.json') or result['checkpoint_sha256']!=digest(folder/'latest_state.npz'):raise ValueError('Changed completed long result')
        endpoint=endpoint_assay(root,job,job_view(p,job),h[-1])
        t=np.array([r['elapsed'] for r in h]);g=np.array([r['uniform_growth_max'] for r in h]);sd=np.array([r['log_activator_sd'] for r in h])
        rows.append(dict(job=job,persistent_contrast=result['persistent_contrast'],late_min_log_sd=result['late_min_log_activator_sd'],
            final_log_sd=float(sd[-1]),uniform_growth_final=float(g[-1]),growth_crossing_time=first_crossing(t,g),
            contrast_onset_time=first_crossing(t,sd,p['criteria']['late_log_sd_min'],'up'),frozen=byref[(job['seed'],job['family'])],
            endpoint_phase=endpoint['phase'],endpoint_numerical_pass=endpoint['numerical_pass'],endpoint_bistability=endpoint['local_bistability_supported']))
    report=dict(completed=len(rows),total=len(p['jobs']),independent_histories=3,results=rows,failures=failures,
                endpoint_numerical_pass=bool(len(rows)==len(p['jobs']) and all(r['endpoint_numerical_pass'] for r in rows)))
    write_json(root/'summary.json',report);return report


def run_stages(root,p):
    frozen_references(root,p)
    if not (root/'pilot-refinement.json').exists():
        for chi in p['contrasts']:
            for level in ('coarse','fine'):
                job=next(j for j in p['jobs'] if j['seed']==p['refinement_seed'] and j['family']=='uniform' and j['level']==level and j['point']['chi']==chi)
                write_json(root/'status.json',dict(state='running',stage='timestep_pilot',current_job=job['key'],total=len(p['jobs'])))
                advance(root,job,p,p['pilot_duration'])
    gate=compare_level(root,p,p['pilot_duration'])
    if not gate['passed']:
        write_json(root/'status.json',dict(state='stopped_at_numerical_gate',stage='timestep_pilot',passed=False,total=len(p['jobs'])))
        return gate
    # Complete the three refined pairs first, then the remaining history cohort.
    ordered=sorted(p['jobs'],key=lambda j:(j['seed']!=p['refinement_seed'],j['point']['chi'],j['level']=='fine',j['family']=='pattern',j['seed']))
    for job in ordered:
        write_json(root/'status.json',dict(state='running',stage='long_continuations',current_job=job['key'],
            completed=sum((root/j['key']/'result.json').exists() for j in p['jobs']),total=len(p['jobs']),pilot_refinement_pass=True))
        try:advance(root,job,p,p['duration'])
        except (RuntimeError,FloatingPointError) as error:
            folder=root/job['key'];folder.mkdir(exist_ok=True);write_json(folder/'status.json',dict(state='failed',error=str(error)))
            print(f"FAILED {job['key']}: {error}",flush=True)
        assess(root)
    refinement=compare_level(root,p,p['duration']) if all((root/j['key']/'result.json').exists() for j in p['jobs'] if j['seed']==p['refinement_seed'] and j['family']=='uniform') else dict(passed=False,incomplete=True)
    summary=assess(root)
    passed=summary['completed']==len(p['jobs']) and not summary['failures'] and summary['endpoint_numerical_pass'] and refinement['passed']
    write_json(root/'status.json',dict(state='completed' if passed else 'completed_with_unresolved_checks',completed=summary['completed'],total=len(p['jobs']),long_refinement_pass=refinement['passed']))
    return summary


def run(root):
    import torch
    root=Path(root);p=read(root/'protocol.json');verify(p);torch.set_num_threads(1)
    for level,validation in p['validation_roots'].items():
        accepted,_=require_validation(validation)
        if asdict(accepted)!=p['accepted_configs'][level]:raise ValueError('Changed accepted backend baseline')
    with (root/'coordinator.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        try:return run_stages(root,p)
        except Exception as error:
            prior=read(root/'status.json')
            write_json(root/'status.json',dict(**{k:v for k,v in prior.items() if k not in ('state','error')},state='failed',error=str(error)))
            raise


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=('prepare','run','assess'))
    parser.add_argument('--output',type=Path,default=Path('outputs/polarity-robustness'))
    parser.add_argument('--parent',type=Path,default=Path('outputs/parameter-robustness-moving'))
    args=parser.parse_args()
    result=prepare(args.output,args.parent) if args.action=='prepare' else run(args.output) if args.action=='run' else assess(args.output)
    print({k:v for k,v in result.items() if k not in ('jobs','results','source_sha256','input_sha256','accepted_configs')})
