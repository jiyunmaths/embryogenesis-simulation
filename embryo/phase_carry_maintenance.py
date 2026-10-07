"""Carry-corrected maintenance on evolving mature geometry, with matched initiation reuse.

All kernels and earlier protocols remain unchanged. The two chemical preparations
are fixed experimental inputs, not independent histories or cell identities.
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
import torch

from .attribute_development import AttributeSimulation
from .benchmark_gpu_backend import cpu_audit
from .cell_response_moving import observe
from .feedback_long import digest, save_checkpoint, restore_checkpoint
from .geometry_precision import exact_state as baseline_exact_state
from .gpu_backend import GpuSimulation
from .gpu_precision_control import PrecisionSimulation, library
from .mechanics_float64_reference import run as reference_run, CRITERIA as REFERENCE_CRITERIA
from .moving_initiation_controls import cpu_chemical_step, exact_state as carry_exact_state
from .native_mechanics import NativeSimulation
from .neighbor_context import read, validate_graph
from .parameter_robustness import verify
from .parameter_robustness_moving import initialize, retain, endpoint_assay
from .phase_carry_polarity import current_evidence as initiation_evidence, trajectory_summary
from .polarity_conductance_controls import check_config
from .polarity_robustness import assert_same_physical_start, first_crossing, refinement_comparison
from .resolution import _steps, write_json
from .validate_gpu_backend import CRITERIA as GPU_CRITERIA, discrepancies

HISTORIES = (9, 7, 8)
CONTRASTS = (0., .35)
FAMILIES = ('pattern', 'uniform')


def job_design(sources):
    return [dict(key=f'seed-{seed}_{level}_chi-{chi:g}_{family}', seed=seed,
        level=level, dt=dt, family=family, arm='phase_carry',
        point=dict(key=f'polarity-{chi:g}', ratio=27.5, chi=chi),
        source=sources[seed][level], chemical_file=sources[seed]['chemical'])
        for seed in HISTORIES for chi in CONTRASTS for family in FAMILIES
        for level, dt in (('coarse', .00375), ('fine', .001875))]


def context_key(job):
    return f'chi-{job["point"]["chi"]:g}_{job["family"]}'


def validate_history(history, p, job, horizon):
    if job['family'] not in FAMILIES or job['arm'] != 'phase_carry':
        raise ValueError('Changed chemical preparation or numerical arm')
    times = np.array([0.]) if horizon == 0. else np.arange(_steps(horizon,p['interval'])+1)*p['interval']
    if len(history)!=len(times) or not np.allclose([r['elapsed'] for r in history],times,rtol=0,atol=1e-9):
        raise ValueError('Incomplete or misaligned observation clocks')
    with np.load(job['chemical_file']) as z:
        if not np.array_equal(history[0]['chemistry'],z[job['family']]):
            raise ValueError('Changed chemical preparation')
        ids=z['ids']
    for row in history:
        if not np.array_equal(ids,row['ids']) or abs(row['time']-p['start']-row['elapsed'])>1e-9:
            raise ValueError('Changed cell IDs or physical clock')
        chem=np.array(row['chemistry'])
        if chem.shape!=(2,len(ids)) or not np.isfinite(chem).all() or np.any(chem<=0):
            raise ValueError('Invalid chemistry')
        if abs(float(np.std(np.log(chem[0])))-row['log_activator_sd'])>1e-12:
            raise ValueError('Chemical contrast does not match saved state')
        validate_graph(row['delta'],row['volumes'])
    return history


def maintenance_summary(history,p):
    result=trajectory_summary(history,p)
    times=np.array([r['elapsed'] for r in history]);sd=np.array([r['log_activator_sd'] for r in history])
    threshold=p['criteria']['late_log_sd_min']
    return dict(**result, initial_log_sd=float(sd[0]),whole_window_min_log_sd=float(sd.min()),
                whole_window_retained=bool(sd.min()>threshold),
                first_contrast_loss_elapsed=first_crossing(times,sd,threshold,'down'))


def maintenance_comparison(zero,positive,qualified):
    if not qualified:return 'unresolved_numerical_or_physical_checks'
    if zero and positive:return 'both_contrasts_maintain'
    if zero:return 'positive_contrast_loses_maintenance'
    if positive:return 'positive_contrast_only_maintains'
    return 'neither_contrast_maintains'


def reused_evidence(p,job):
    item=p['reused'][job['key']];old_job=item['job'];root=Path(item['root'])
    if (job['family']!='uniform' or any(job[k]!=old_job[k] for k in ('seed','level','dt','family','arm')) or
            job['point']!=old_job['point'] or digest(root/'protocol.json')!=item['protocol_sha256']):
        raise ValueError('Reused path is not the matched carry initiation control')
    old=read(root/'protocol.json');result,h=initiation_evidence(root,old,old_job)
    validate_history(h,p,job,p['duration'])
    with np.load(root/old_job['key']/'latest_state.npz') as z:
        import json
        check_config(json.loads(str(z['metadata']))['config'],job,p)
    return result,h


def current_evidence(root,p,job):
    folder=Path(root)/job['key'];ph=digest(Path(root)/'protocol.json')
    result,h=read(folder/'result.json'),read(folder/'history.json')
    if (result['job']!=job or result['protocol_sha256']!=ph or not result['quality_pass'] or
            result['history_sha256']!=digest(folder/'history.json') or
            result['checkpoint_sha256']!=digest(folder/'latest_state.npz')):
        raise ValueError('Changed completed maintenance evidence')
    host,audit,saved=restore_checkpoint(folder/'latest_state.npz',ph,job)
    validate_history(h,p,job,p['duration']);check_config(asdict(host.config),job,p)
    if (saved!=h or audit!=result['audit'] or abs(host.time-p['start']-p['duration'])>1e-9 or
            host.step_number!=_steps(p['start']+p['duration'],job['dt']) or
            not np.array_equal(h[-1]['chemistry'],[host.activator,host.inhibitor])):
        raise ValueError('Maintenance checkpoint/history mismatch')
    with np.load(folder/'latest_state.npz') as z:
        if (str(z['precision_arm'])!='phase_carry' or z['phase_carry'].dtype!=np.float64 or
                z['phase_carry'].shape!=host.phi.shape or not np.isfinite(z['phase_carry']).all()):
            raise ValueError('Invalid maintenance carry residual')
    return result,h


def actual_context_gate(folder,job,p):
    """Strict native/carry prefix, original arithmetic check, chemistry and restart."""
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
    config=p['accepted_configs'][job['level']]
    make=lambda:initialize(job['source'],job['chemical_file'],job['family'],job['point'],config)
    original=GpuSimulation(make());baseline=PrecisionSimulation(make());carried=PrecisionSimulation(make(),'phase_carry')
    original.step();baseline.step();carried.step()
    baseline_exact_state(original,baseline);baseline_exact_state(original,carried)
    for _ in range(15):original.step();baseline.step();baseline_exact_state(original,baseline)
    del original,baseline,carried
    sim=PrecisionSimulation(make(),'phase_carry');chemical_error=accounting=0.
    for _ in range(8):
        delta=sim.matrices()[1].cpu().numpy();old=sim.geometry[:,0].cpu().numpy().copy()
        before=np.array([sim.activator.cpu().numpy(),sim.inhibitor.cpu().numpy()])
        expected=cpu_chemical_step(before,delta,sim.config.dt,2.,.02,.55)
        sim.step();expected*=old/sim.geometry[:,0].cpu().numpy()
        actual=np.array([sim.activator.cpu().numpy(),sim.inhibitor.cpu().numpy()])
        chemical_error=max(chemical_error,float(abs(np.log(actual/expected)).max()))
        accounting=max(accounting,sim.audit()['dilution_amount_error'])
    if chemical_error>p['chemical_reference_log_max'] or accounting>p['criteria']['dilution_error_max']:
        raise ValueError('Developed-context chemistry/accounting failed')
    checkpoint=folder/'carry-restart.npz';sim.checkpoint(checkpoint);restart=PrecisionSimulation.restore(checkpoint)
    carry_exact_state(sim,restart)
    for _ in range(4):sim.step();restart.step();carry_exact_state(sim,restart)
    del sim,restart
    cpu=make();cpu.__class__=NativeSimulation;cpu.native_threads=p['prefix_native_threads']
    gpu=PrecisionSimulation(make(),'phase_carry');histories=[[],[]]
    errors=dict.fromkeys(('chemical_log_max','polarity_abs_max','relative_axis_max','relative_volume_max','relative_transport_max'),0.)
    steps=_steps(p['prefix_duration'],job['dt']);every=_steps(p['interval'],job['dt'])
    for step in range(steps+1):
        cpu_audit(cpu);quality=gpu.audit()
        if quality['dilution_amount_error']>p['criteria']['dilution_error_max']:
            raise ValueError('Prefix amount accounting failed')
        if step%every==0:
            a,b=gpu.observe(step*job['dt']),observe(cpu,step*job['dt'])
            if max(a['boundary_occupancy'],b['boundary_occupancy'])>=p['criteria']['boundary_max']:
                raise ValueError('Prefix boundary failed')
            for key,value in discrepancies(a,b).items():errors[key]=max(errors[key],value)
            histories[0].append(a);histories[1].append(b)
        if step<steps:cpu.step();gpu.step()
    phi_error=float(abs(gpu.phi.cpu().numpy()-cpu.phi).max())
    passed=all(value<=p['prefix_criteria'][key] for key,value in errors.items()) and phi_error<=p['prefix_criteria']['final_phi_abs_max']
    files=[folder/'gpu-history.json',folder/'cpu-history.json',folder/'gpu-end.npz',folder/'cpu-end.npz',checkpoint]
    write_json(files[0],histories[0]);write_json(files[1],histories[1]);gpu.checkpoint(files[2]);cpu.checkpoint(files[3])
    del gpu,cpu;torch.cuda.empty_cache()
    result=dict(passed=bool(passed),job=job,baseline_exact_steps=16,zero_carry_first_step_exact=True,
        independent_chemical_steps=8,cpu_chemical_log_max=chemical_error,amount_accounting_max=accounting,
        carry_restart_exact_steps=4,native_carry_prefix_duration=p['prefix_duration'],errors=errors,phi_abs_max=phi_error,
        evidence_sha256={str(f.resolve()):digest(f) for f in files},
        scope='Exact developed preparation/geometry/parameters; strict short native/carry coupling comparison, not full-horizon CPU/GPU equivalence.')
    write_json(folder/'result.json',result)
    if not passed:raise RuntimeError('Developed-state native/carry prefix fails original limits')
    return result


def context_gate(root,job,p):
    file=Path(root)/job['key']/'context-gate.json';r=read(file)
    if r['job']!=job or not r['passed'] or r['protocol_sha256']!=digest(Path(root)/'protocol.json'):
        raise ValueError('Changed prelaunch context gate')
    for f,h in r['evidence_sha256'].items():
        if digest(f)!=h:raise ValueError('Changed prelaunch actual-context evidence')
    return r


def pair_report(root,p,job,horizon,histories):
    folder=Path(root)/'pairs'/f'seed-{job["seed"]}_{context_key(job)}';folder.mkdir(parents=True,exist_ok=True)
    label='pilot' if horizon==p['pilot_duration'] else 'full';file=folder/f'{label}-refinement.json';ph=digest(Path(root)/'protocol.json')
    if file.exists():
        r=read(file)
        if r['protocol_sha256']!=ph:raise ValueError('Changed pair protocol')
        for f,h in r['history_sha256'].items():
            if digest(f)!=h:raise ValueError('Changed immutable pair history')
        for level in ('coarse','fine'):
            if read(folder/f'{label}-{level}.json')!=histories[level]:raise ValueError('Changed compared history')
        return r
    r=refinement_comparison(histories['coarse'],histories['fine'],p,horizon)
    if job['family']=='pattern':
        threshold=p['criteria']['late_log_sd_min']
        losses=[first_crossing([x['elapsed'] for x in histories[level]],[x['log_activator_sd'] for x in histories[level]],threshold,'down') for level in ('coarse','fine')]
        error=0. if losses==[None,None] else None if None in losses else abs(losses[0]-losses[1])
        retained=[all(x['log_activator_sd']>threshold for x in histories[level]) for level in ('coarse','fine')]
        r.update(loss_times=losses,loss_time_error=error,whole_window_retained=retained)
        r['passed']=bool(r['passed'] and error is not None and error<=p['loss_time_abs_max'] and retained[0]==retained[1])
    files=[]
    for level in ('coarse','fine'):
        target=folder/f'{label}-{level}.json';write_json(target,histories[level]);files.append(target)
    r.update(seed=job['seed'],chi=job['point']['chi'],family=job['family'],protocol_sha256=ph,
        history_sha256={str(f.resolve()):digest(f) for f in files},
        interpretation=f'Raw matched-time carry refinement in history {job["seed"]}, {job["family"]} start, chi={job["point"]["chi"]:g}. Frozen spectra do not describe full moving-system stability.')
    write_json(file,r);return r
def advance(root, job, p, horizon):
    """Stop/resume at declared horizons without publishing a pilot as a result."""
    if horizon not in (p['pilot_duration'],p['duration']): raise ValueError('Undeclared horizon')
    root=Path(root); reuse=p['reused'].get(job['key'])
    if reuse:
        _,h=reused_evidence(p,job)
        return validate_history([r for r in h if r['elapsed']<=horizon+1e-9],p,job,horizon)
    context_gate(root,job,p)
    folder=root/job['key']; cp=folder/'latest_state.npz'; ph=digest(root/'protocol.json')
    if (folder/'result.json').exists():
        _,h=current_evidence(root,p,job)
        return validate_history([r for r in h if r['elapsed']<=horizon+1e-9],p,job,horizon)
    if cp.exists():
        host,audit,history=restore_checkpoint(cp,ph,job); sim=PrecisionSimulation.restore(cp)
        check_config(asdict(host.config),job,p)
        if sim.precision_arm != 'phase_carry': raise ValueError('Carry arm changed')
        validate_history(history,p,job,sim.time-p['start'])
        if not np.array_equal(history[-1]['chemistry'],[host.activator,host.inhibitor]): raise ValueError('Checkpoint chemistry changed')
    else:
        sim=PrecisionSimulation(initialize(job['source'],job['chemical_file'],job['family'],job['point'],p['accepted_configs'][job['level']]),'phase_carry')
        history=[];audit=dict(max_volume_error=0.,min_radius=1e100,max_clipping=0.,boundary_max=0.,dilution_error_max=0.,wall_seconds=0.)
    start=_steps(p['start'],job['dt']); stop=_steps(p['start']+horizon,job['dt'])
    every=_steps(p['interval'],job['dt']); save_every=_steps(p['checkpoint_interval'],job['dt'])
    if not start<=sim.step_number<=_steps(p['start']+p['duration'],job['dt']) or abs(sim.time-sim.step_number*job['dt'])>1e-9:
        raise ValueError('Invalid carry checkpoint clock')
    if sim.step_number>stop:
        del sim;torch.cuda.empty_cache()
        return validate_history([r for r in history if r['elapsed']<=horizon+1e-9],p,job,horizon)
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
            write_json(folder/'status.json',dict(state='running',elapsed=elapsed,target_horizon=horizon,audit=audit,log_activator_sd=history[-1]['log_activator_sd']))
            if offset%save_every==0 or sim.step_number==stop:
                save_checkpoint(sim,cp,audit,history,ph,job)
                write_json(folder/'history.json',history)
                print(f'{job["key"]}: elapsed={elapsed:g}, target={horizon:g}, wall={audit["wall_seconds"]:.1f}s',flush=True)
        if sim.step_number==stop:break
        sim.step()
    validate_history(history,p,job,horizon)
    if horizon==p['duration']:
        write_json(folder/'result.json',dict(job=job,protocol_sha256=ph,history_sha256=digest(folder/'history.json'),
            checkpoint_sha256=digest(cp),quality_pass=True,audit=audit,**maintenance_summary(history,p)))
        write_json(folder/'status.json',dict(state='completed',elapsed=horizon))
    else:write_json(folder/'status.json',dict(state='awaiting_pair_refinement',elapsed=horizon,audit=audit))
    del sim;torch.cuda.empty_cache();return history



def prepare(root,parent=Path('outputs/polarity-conductance-controls'),phase=Path('outputs/phase-carry-polarity')):
    root,parent,phase=[Path(x).resolve() for x in (root,parent,phase)]
    if root.exists():raise FileExistsError(root)
    old,ip=read(parent/'protocol.json'),read(phase/'protocol.json');verify(old);verify(ip)
    review=Path('docs/polarity_conductance_controls_assessment.json');e=read(review)
    if (read(parent/'status.json')['state']!='completed' or not read(parent/'summary.json')['passed'] or
            e['protocol_sha256']!=digest(parent/'protocol.json') or
            e['original_summary_sha256']!=digest(parent/'summary.json') or
            digest(e['full_review'])!=e['full_review_sha256'] or
            not read(phase/'summary.json')['passed']):
        raise ValueError('Requires qualified completed and reviewed initiation controls')
    if (torch.cuda.get_device_name(0)!=old['device']['name'] or torch.__version__!=old['device']['torch'] or
            torch.version.cuda!=old['device']['torch_cuda'] or digest(library()._name)!=old['device']['experimental_binary_sha256']):
        raise ValueError('Changed qualified GPU/runtime/kernel')
    torch.set_num_threads(1);root.mkdir(parents=True)
    inputs=[Path('docs/phase_carry_maintenance.md'),review,Path(e['full_review'])]
    inputs += [parent/f for f in ('protocol.json','summary.json','status.json')]
    sources={}
    for seed in HISTORIES:
        folder=root/f'seed-{seed}';folder.mkdir();sources[seed]={}
        for label,name in (('coarse','coarse-source.npz'),('fine','fine-source.npz'),('chemical','initial_states.npz')):
            original,copied=phase/f'seed-{seed}'/name,folder/name
            shutil.copy2(original,copied);inputs += [original,copied];sources[seed][label]=str(copied)
        with np.load(sources[seed]['chemical']) as z:
            if float(np.std(np.log(z['pattern'][0])))<=old['criteria']['late_log_sd_min']:
                raise ValueError('Developed preparation lacks initial contrast')
        a,b=[initialize(sources[seed][level],sources[seed]['chemical'],'pattern',dict(ratio=27.5,chi=.35),old['accepted_configs'][level]) for level in ('coarse','fine')]
        assert_same_physical_start(a,b);del a,b
    jobs=job_design(sources);reused={}
    for job in jobs:
        if job['family']!='uniform':continue
        original=next(j for j in ip['jobs'] if j['seed']==job['seed'] and j['level']==job['level'] and j['point']==job['point'])
        item=ip['reused'].get(original['key'],dict(root=str(phase),job=original))
        reused[job['key']]=dict(**item,protocol_sha256=digest(Path(item['root'])/'protocol.json'))
    p=dict(jobs=jobs,reused=reused,accepted_configs=old['accepted_configs'],device=old['device'],
        histories=list(HISTORIES),contrasts=list(CONTRASTS),families=list(FAMILIES),
        start=150.,duration=240.,pilot_duration=60.,late_window=24.,pilot_late_window=12.,
        interval=.15,checkpoint_interval=3.,independent_histories=3,new_histories=0,
        new_moving_jobs=12,reused_trajectories=12,endpoint_workers=2,
        criteria=old['criteria'],refinement_criteria=old['refinement_criteria'],loss_time_abs_max=.3,
        prefix_duration=.6,prefix_native_threads=4,prefix_criteria=GPU_CRITERIA,
        chemical_reference_log_max=1e-11,
        endpoint_horizons=old['endpoint_horizons'],endpoint_interval=old['endpoint_interval'],endpoint_criteria=old['endpoint_criteria'],
        design='Same mature t=150 geometry/polarity, developed and near-uniform chemical preparations, histories 7/8/9. New developed starts co-evolve chemistry/mechanics with chi=0/.35 at dt=.00375/.001875 to elapsed240. Reuse twelve qualified uniform paths. Evolving conservative conductances and actual volumes/dilution in every path; no fixed transport, new noise, nutrient input, fate switch, chemical preconditioning or fresh zygote history.',
        interpretation='Three existing developmental histories; twenty-four nested paths, twelve new developed maintenance continuations and twelve reused near-uniform initiation controls. Report maintenance at each contrast separately from formation, continuous retention separately from late-window recovery, and frozen endpoint coexistence separately from moving contrast. Loss is an admissible scientific outcome; numerical failure remains unresolved. Pattern survival does not establish cell identity, autonomy, inheritance or full moving-system stability. Mature starts and supplied developed preparations inherit old-method development; no broader ledger, spatial/closure or backend promotion.',
        estimated_new_moving_seconds=1.5*sum(r['audit']['wall_seconds'] for r in read(e['full_review'])['completed_paths'] if not r['reused']))
    gates=[]
    try:
        for job in jobs:
            if job['key'] in reused:
                _,h=reused_evidence(p,job)
                item=reused[job['key']];src=Path(item['root'])/item['job']['key'];endpoint=src/'endpoint'
                ep=read(endpoint/'protocol.json');verify(ep);assay=read(endpoint/'assay/result.json')
                if (not assay['numerical_pass'] or not assay['all_trials_settled'] or
                        assay['protocol_sha256']!=digest(endpoint/'protocol.json') or
                        assay['paths_sha256']!=digest(endpoint/'assay/paths.npz') or
                        ep['moving_protocol_sha256']!=item['protocol_sha256'] or
                        ep['moving_history_sha256']!=digest(src/'history.json')):
                    raise ValueError('Changed reused initiation endpoint')
                inputs += [Path(item['root'])/'protocol.json']+[src/f for f in ('result.json','history.json','latest_state.npz','endpoint/protocol.json','endpoint/source.npz','endpoint/assay/result.json','endpoint/assay/paths.npz')]
                continue
            write_json(root/'status.json',dict(state='preparing',stage='actual_developed_context',current_job=job['key'],context_checks_completed=len(gates),total_context_checks=12))
            folder=root/job['key']/'prelaunch';g=actual_context_gate(folder,job,p);gates.append(g)
            inputs += [folder/'result.json',*[Path(f) for f in g['evidence_sha256']]]
            print('Developed-context gate passed:',job['key'],flush=True)
        # Independent float64 mechanics at the prespecified history-9 positive-chi
        # developed start; all twelve contexts additionally receive native/carry checks.
        ref_job=next(j for j in jobs if j['seed']==9 and j['point']['chi']==.35 and j['level']=='coarse' and j['family']=='pattern')
        ref_host=initialize(ref_job['source'],ref_job['chemical_file'],'pattern',ref_job['point'],p['accepted_configs']['coarse'])
        ref_source=root/'developed-reference-source.npz';ref_host.checkpoint(ref_source);inputs.append(ref_source)
        contexts=[dict(key='history-9-developed-positive-chi',source=str(ref_source))]
        write_json(root/'status.json',dict(state='preparing',stage='independent_float64_mechanics',context_checks_completed=12,total_context_checks=12))
        r=reference_run(root/'float64-reference',contexts,REFERENCE_CRITERIA,.15,(.00375,.001875,.0009375))
        r['scope']='One developed positive-chi state from existing history 9; three nested timestep comparisons of mechanics only. Independent float64 fields/geometry with concentrations and polarity held fixed. This is a short component check, not three histories, full coupled stability or full-horizon backend equivalence.'
        write_json(root/'float64-reference/result.json',r)
        if not r['passed']:raise ValueError('Developed-state independent float64 mechanics check failed')
        inputs += [root/'float64-reference/result.json',*[Path(f) for f in r['field_sha256']]]
        p.update(reference_contexts=contexts,reference_criteria=REFERENCE_CRITERIA,reference_horizon=.15,reference_dts=[.00375,.001875,.0009375])
        gate_file=root/'implementation-verification.json';write_json(gate_file,dict(passed=True,actual_contexts=gates,reference=dict(passed=True,result_sha256=digest(root/'float64-reference/result.json'))));inputs.append(gate_file)
        files=[Path(__file__),Path('tests/test_phase_carry_maintenance.py')]
        p['source_sha256']={**old['source_sha256'],**{str(f.resolve()):digest(f) for f in files}}
        p['input_sha256']={**old['input_sha256'],**{str(f.resolve()):digest(f) for f in inputs}}
        write_json(root/'protocol.json',p);ph=digest(root/'protocol.json')
        for g in gates:
            write_json(root/g['job']['key']/'context-gate.json',dict(**g,protocol_sha256=ph))
        for seed in HISTORIES:
            for chi in CONTRASTS:
                matched=[j for j in jobs if j['seed']==seed and j['point']['chi']==chi and j['family']=='uniform']
                h={j['level']:reused_evidence(p,j)[1] for j in matched}
                for horizon in (p['pilot_duration'],p['duration']):
                    cut={level:[row for row in history if row['elapsed']<=horizon+1e-9] for level,history in h.items()}
                    if not pair_report(root,p,matched[0],horizon,cut)['passed']:
                        raise ValueError('Matched reused initiation comparison fails original limits')
        write_json(root/'status.json',dict(state='prepared',completed=0,total=12,reused_trajectories=12,actual_context_checks=12,reference_pass=True,estimated_new_moving_seconds=p['estimated_new_moving_seconds']))
        return p
    except Exception as error:
        write_json(root/'status.json',dict(state='preparation_failed',error=str(error)));raise


def verify_reference(root,p):
    r=read(Path(root)/'float64-reference/result.json')
    if (not r['passed'] or r['contexts']!=p['reference_contexts'] or r['criteria']!=p['reference_criteria'] or
            r['horizon']!=p['reference_horizon'] or r['dts']!=p['reference_dts']):
        raise ValueError('Changed independent developed-context reference')
    for f,h in r['field_sha256'].items():
        if digest(f)!=h:raise ValueError('Changed independent reference field')
    return r


def run_pairs(root,p,executor):
    done,pending,failed=0,[],[]
    for seed in p['histories']:
        for chi in p['contrasts']:
            jobs=[j for j in p['jobs'] if j['seed']==seed and j['point']['chi']==chi and j['family']=='pattern']
            h={}
            for job in jobs:
                write_json(Path(root)/'status.json',dict(state='running',stage='matched_maintenance_pilot',current_job=job['key'],completed=done,total=12,failed_pilot_contexts=failed))
                h[job['level']]=advance(root,job,p,p['pilot_duration'])
            if not pair_report(root,p,jobs[0],p['pilot_duration'],h)['passed']:
                failed.append(dict(seed=seed,chi=chi,family='pattern'));continue
            for job in jobs:
                write_json(Path(root)/'status.json',dict(state='running',stage='long_moving_maintenance',current_job=job['key'],completed=done,total=12,failed_pilot_contexts=failed))
                h[job['level']]=advance(root,job,p,p['duration']);done+=1
                pending.append(executor.submit(endpoint_assay,root,job,p,h[job['level']][-1]))
            pair_report(root,p,jobs[0],p['duration'],h)
    return done,pending,failed


def assess(root,p):
    root=Path(root);verify(p);verify_reference(root,p)
    rows,endpoints,missing=[],[],[]
    for job in p['jobs']:
        if job['key'] in p['reused']:
            _,h=reused_evidence(p,job);item=p['reused'][job['key']];src=Path(item['root'])/item['job']['key']
            ep=read(src/'endpoint/protocol.json');verify(ep);e=read(src/'endpoint/assay/result.json')
            if e['protocol_sha256']!=digest(src/'endpoint/protocol.json') or e['paths_sha256']!=digest(src/'endpoint/assay/paths.npz'):
                raise ValueError('Changed reused endpoint')
        elif not (root/job['key']/'result.json').exists():missing.append(job['key']);continue
        else:
            _,h=current_evidence(root,p,job);e=endpoint_assay(root,job,p,h[-1])
        rows.append(dict(job=job,reused=job['key'] in p['reused'],**maintenance_summary(h,p)))
        endpoints.append(dict(job=job['key'],numerical_pass=e['numerical_pass'],all_trials_settled=e['all_trials_settled'],phase=e['phase'],local_bistability=e['local_bistability_supported']))
    by_history=[]
    for seed in p['histories']:
        selected=[e for e in endpoints if next(j for j in p['jobs'] if j['key']==e['job'])['seed']==seed]
        qualified=len(selected)==8 and all(e['numerical_pass'] and e['all_trials_settled'] for e in selected)
        branches={}
        for chi in p['contrasts']:
            for family in FAMILIES:
                job=next(j for j in p['jobs'] if j['seed']==seed and j['point']['chi']==chi and j['family']==family)
                file=root/'pairs'/f'seed-{seed}_{context_key(job)}'/'full-refinement.json'
                if not file.exists():qualified=False;continue
                h={level:read(file.parent/f'full-{level}.json') for level in ('coarse','fine')}
                qualified &= pair_report(root,p,job,p['duration'],h)['passed']
                matched=[r for r in rows if r['job']['seed']==seed and r['job']['point']['chi']==chi and r['job']['family']==family]
                if len(matched)!=2:qualified=False;continue
                branches[context_key(job)]=next(r['persistent_contrast'] for r in matched if r['job']['level']=='fine')
        by_history.append(dict(seed=seed,qualified=bool(qualified),late_contrast=branches,
            maintenance_conclusion=maintenance_comparison(branches.get('chi-0_pattern',False),branches.get('chi-0.35_pattern',False),qualified)))
    passed=not missing and all(r['qualified'] for r in by_history)
    result=dict(protocol_sha256=digest(root/'protocol.json'),passed=bool(passed),completed=sum(not r['reused'] for r in rows),total=12,
        independent_histories=3,new_histories=0,reused_trajectories=12,results=rows,endpoints=endpoints,missing=missing,by_history=by_history,scope=p['interpretation'])
    write_json(root/'summary.json',result);return result


def run(root):
    root=Path(root).resolve()
    with (root/'coordinator.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        p=read(root/'protocol.json');verify(p);verify_reference(root,p);torch.set_num_threads(1)
        if (torch.cuda.get_device_name(0)!=p['device']['name'] or torch.__version__!=p['device']['torch'] or
                torch.version.cuda!=p['device']['torch_cuda'] or digest(library()._name)!=p['device']['experimental_binary_sha256']):
            raise ValueError('Changed qualified GPU/runtime/kernel')
        try:
            with ProcessPoolExecutor(max_workers=p['endpoint_workers'],mp_context=multiprocessing.get_context('spawn')) as executor:
                done,pending,failed=run_pairs(root,p,executor)
                write_json(root/'status.json',dict(state='running',stage='endpoint_assays',completed=done,total=12,failed_pilot_contexts=failed))
                for future in pending:future.result()
            result=assess(root,p)
            write_json(root/'status.json',dict(state='completed' if result['passed'] else 'completed_with_unresolved_checks',completed=result['completed'],total=12,passed=result['passed'],failed_pilot_contexts=failed))
            return result
        except Exception as error:
            write_json(root/'status.json',dict(state='failed',error=str(error)));raise


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=('prepare','run','assess'))
    parser.add_argument('--output',type=Path,default=Path('outputs/phase-carry-maintenance'))
    args=parser.parse_args()
    result=prepare(args.output) if args.action=='prepare' else run(args.output) if args.action=='run' else assess(args.output,read(args.output/'protocol.json'))
    print({k:v for k,v in result.items() if k in ('passed','completed','total','estimated_new_moving_seconds','by_history')})
