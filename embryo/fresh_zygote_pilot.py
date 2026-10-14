"""Prospectively matched new histories using the unchanged developmental solver.

This coordinator executes only the t=0--18 developmental stage. Long-window
qualification is a separate, explicitly gated continuation of its saved states.
"""
import argparse
from dataclasses import asdict, replace
from pathlib import Path
import fcntl
import hashlib
import json
import os
import time

import numpy as np
import torch

from . import development_control_qualification as controls
from . import development_space_refinement as space
from . import development_time_diagnosis as td
from . import qualify_repeated_development as q
from .model import Config
from .qualify_development_gpu import digest, read, write, pin, save_cpu, reference, exact
from .resolution import _steps

HISTORIES = (1101, 1102, 1103)
BRANCHES = ('off', 'tension_adhesion', 'polarity', 'full')
LEVELS = (('development', .00046875), ('development_half', .000234375))
PARENT = Path('output/studies/development-control-qualification')
ASSESSMENT = Path('docs/development_control_qualification_assessment.json')
PREFIX_LIMITS = dict(phase_abs_max=2e-5, effective_phase_abs_max=2e-5,
    chemical_log_max=1e-5, polarity_abs_max=1e-5,
    relative_volume_max=1e-5, relative_transport_max=1e-5)


def design(base):
    jobs = []
    for seed in HISTORIES:
        for branch in BRANCHES:
            feedback, tension, adhesion, polarity = controls.COUPLINGS[branch]
            for level, dt in LEVELS:
                cfg = replace(Config(**base), seed=seed, grid=120, dt=dt,
                    steps=_steps(18., dt), save_every=_steps(.15, dt),
                    feedback=feedback, fate_tension=tension, fate_adhesion=adhesion,
                    polarity_tension=polarity)
                key = f'h{seed}_{branch}_{level}'
                jobs.append(dict(key=key, seed=seed, branch=branch, level=level,
                    config=asdict(cfg), reused=False))
    return jobs


def config_changes(base, candidate, branch):
    changed = {key for key, value in candidate.items() if value != base[key]}
    allowed = {'seed', 'dt', 'steps', 'save_every', 'feedback',
               'fate_tension', 'fate_adhesion', 'polarity_tension'}
    if not changed <= allowed:
        raise ValueError('Physical model drift: '+str(changed-allowed))
    if candidate['seed'] not in HISTORIES or candidate['grid'] != 120:
        raise ValueError('Unregistered history or grid')
    if candidate['dt'] not in dict(LEVELS).values():
        raise ValueError('Unregistered developmental timestep')
    controls.config_changes(base, dict(candidate, seed=base['seed']), branch)
    return sorted(changed)


def initial_signature(path, job):
    """Exact matching of state/RNG after allowing the declared configuration/mode.

    Array hashes include dtype and shape. Metadata other than config is included,
    so the first cycle time, lineage and all RNG states must match within history.
    """
    with np.load(path, allow_pickle=False) as z:
        metadata = json.loads(str(z['metadata']))
        if metadata['config'] != job['config']:
            raise ValueError('Initial checkpoint configuration differs from job')
        if (metadata['time'] != 0 or metadata['step_number'] != 0 or
                metadata['divisions'] or len(z['ids']) != 1 or
                np.any(z['fate'] != 0) or np.any(z['polarity'] != 0) or
                np.any(z['activator'] != 1.03) or np.any(z['inhibitor'] != .98) or
                np.any(z['phase_carry'] != 0) or np.any(z['chemical_correction'] != 0) or
                np.any(z['rounding'] != 0) or json.loads(str(z['event_audits']))):
            raise ValueError('Not the declared fresh, unlabeled zygote')
        expected_mode = 'no_feedback' if job['branch'] == 'off' else 'direct'
        if str(z['attribute_mode']) != expected_mode:
            raise ValueError('Initial mechanical mode differs from branch')
        if str(z['chemical_method']) != td.gpu.METHOD or int(z['development_schema']) != td.gpu.SCHEMA:
            raise ValueError('Initial correction schema differs from qualified backend')
        del metadata['config']
        signature = dict(metadata=metadata, arrays={})
        for name in sorted(set(z.files)-{'metadata', 'attribute_mode'}):
            x = np.ascontiguousarray(z[name])
            signature['arrays'][name] = dict(dtype=str(x.dtype), shape=list(x.shape),
                sha256=hashlib.sha256(x.tobytes()).hexdigest())
        return signature


def match_starts(jobs):
    signatures = {}; reference_signatures = {}
    for job in jobs:
        sig = initial_signature(job['source'], job)
        if job['seed'] in reference_signatures and sig != reference_signatures[job['seed']]:
            raise ValueError('Initial states/random streams are not matched within history')
        reference_signatures[job['seed']] = sig
        signatures[job['key']] = hashlib.sha256(json.dumps(sig, sort_keys=True).encode()).hexdigest()
    if len(set(signatures.values())) != len(HISTORIES):
        raise ValueError('Histories do not have distinct initial random states')
    return signatures


def verify_seed_inventory(inventory):
    if (inventory['proposed_histories'] != list(HISTORIES) or
            not inventory['all_proposed_unused'] or
            set(HISTORIES) & set(inventory['used_seed_or_history_ids'])):
        raise ValueError('Fresh history identifiers were already used or inventory drifted')
    for record in inventory['files']:
        if digest(record['path']) != record['sha256']:
            raise ValueError('Seed inventory input changed: '+record['path'])


def verify_preserved(p):
    q.verify(p)
    for name, expected in read(p['preserved_manifest'])['preserved_sha256'].items():
        if digest(name) != expected:
            raise ValueError('Preserved historical evidence changed: '+name)


def prepare(root, seed_inventory):
    root = Path(root).resolve()
    if root.exists():
        raise FileExistsError('Fresh pilot directory required')
    accepted = read(ASSESSMENT)
    review = Path(accepted['independent_review']['path'])
    if not accepted['review_completed'] or not accepted['scoped_qualification_accepted']:
        raise ValueError('Independently accepted matched developmental controls required')
    if (digest(review) != accepted['independent_review']['sha256'] or
            digest(PARENT/'protocol.json') != accepted['protocol_sha256']):
        raise ValueError('Accepted developmental provenance changed')
    parent = read(PARENT/'protocol.json'); controls.verify_preserved(parent)
    inventory = read(seed_inventory); verify_seed_inventory(inventory)
    baseline = next(j for j in parent['jobs'] if j['key'] == 'off120')['config']
    root.mkdir(parents=True)
    write(root/'seed-inventory.json', inventory)
    inputs = [ASSESSMENT, review, PARENT/'protocol.json', PARENT/'summary.json',
        Path(parent['preserved_manifest']), root/'seed-inventory.json',
        Path(parent['instrumentation_gate'])]
    inputs += [Path(r['path']) for r in inventory['files']]
    jobs = design(baseline)
    for job in jobs:
        folder = root/job['key']; folder.mkdir()
        host = q.initial(Config(**job['config']))
        config_changes(baseline, asdict(host.config), job['branch'])
        source = folder/'source.npz'; save_cpu(host, source)
        job['source'] = str(source); inputs.append(source)
    signatures = match_starts(jobs)
    write(root/'matched-starts.json', dict(passed=True, signatures=signatures,
        independent_histories=3, nested_initial_states=24, initial_cells=1,
        state_and_rng_matched_within_history=True, distinct_random_states_across_histories=True))
    inputs.append(root/'matched-starts.json')
    sources = dict(parent['source_sha256'])
    sources.update(pin([Path(__file__), Path('docs/fresh_zygote_pilot.md'),
        Path('tests/test_fresh_zygote_pilot.py'), Path('scripts/with_study_output_path_alias.py'),
        Path('scripts/with_output_path_alias.py')]))
    p = {key:parent[key] for key in ('duration', 'interval', 'checkpoint_interval',
        'diagnostic_interval', 'snapshots', 'time_limits', 'physical_limits', 'restart_steps',
        'restart_stages', 'minimum_mature_window', 'method', 'gpu_uuid', 'dense_time_spectrum_limit',
        'instrumentation_gate', 'instrumentation_protocol_sha256', 'runtime_binary_sha256')}
    p.update(jobs=jobs, histories=list(HISTORIES), branches=list(BRANCHES),
        levels=[dict(level=level, dt=dt) for level,dt in LEVELS], grid=120,
        independent_histories=3, new_scientific_histories=3, new_numerical_paths=24,
        baseline_config=baseline, matched_start_signatures=signatures,
        prefix_steps=8, prefix_limits=PREFIX_LIMITS,
        source_sha256=sources, input_sha256=pin(inputs),
        accepted_parent_assessment=str(ASSESSMENT.resolve()),
        accepted_parent_review=str(review.resolve()), parent_protocol_sha256=digest(PARENT/'protocol.json'),
        preserved_manifest=parent['preserved_manifest'],
        time_pairs=[dict(seed=seed, branch=branch,
            cases=[f'h{seed}_{branch}_{level}' for level,_ in LEVELS])
            for seed in HISTORIES for branch in BRANCHES],
        stop_policy='Stop on any physical, prefix, restart or timestep gate failure; preserve '
            'failure and all unused registered histories. Never replace histories or change thresholds.',
        post_cleavage_plan=dict(executed_by_this_coordinator=False, inherited_time=18.,
            pilot_elapsed=60., confirmatory_elapsed=240., late_window=24.,
            candidate_dts=[.00375, .001875, .0009375],
            restart_policy='Retain visible phase, phase carry, chemical corrections, event state and RNG. '
                'No chemistry reset, supplied pattern, further division or silent backend handoff.',
            gates='Independent acceptance of developmental data first; actual-start CPU/GPU prefixes, '
                'exact correction-aware restart and matched timestep tests at each new mechanical context. '
                'Qualify larger mature timesteps before committing to the 60/240-unit horizons.',
            endpoints='Frozen actual-endpoint spectra and local uniform/patterned coexistence; '
                'response/context perturbations are a later separately matched experiment.'),
        scope='Three newly seeded zygote histories, four matched mechanical controls, two developmental '
            'timesteps at grid120 through t18. This stage measures development and numerical reproducibility; '
            'no long-window initiation/maintenance, biological identity, 32-cell or continuum acceptance.')
    write(root/'protocol.json', p)
    write(root/'status.json', dict(state='prepared', protocol_sha256=digest(root/'protocol.json')))
    return p


def prefix(root, job, p):
    """Independent NumPy full-state prefix plus exact GPU correction-aware restart."""
    folder = Path(root)/job['key']; output = folder/'actual-start-prefix.json'
    if output.exists():
        raise FileExistsError('Never overwrite a recorded actual-start check')
    sim = td.gpu.DevelopmentGpu.restore(job['source']); cpu = reference(job['source'])
    clone = td.gpu.DevelopmentGpu.restore(job['source']); exact(sim, clone)
    maxima = dict.fromkeys(PREFIX_LIMITS, 0.); failure = None; started=time.perf_counter()
    try:
        for _ in range(p['prefix_steps']):
            sim.step(); cpu.step(); clone.step(); exact(sim, clone)
            if not np.array_equal(sim.ids, cpu.ids) or sim.host.divisions or cpu.divisions:
                raise ValueError('Actual-start prefix unexpectedly crosses cleavage')
            gp = sim.phi.cpu().numpy(); gc = sim.phase_carry.cpu().numpy()
            gu = torch.stack((sim.activator,sim.inhibitor)).cpu().numpy()
            cu = np.stack((cpu.activator,cpu.inhibitor))
            gv = sim.geometry[:,0].cpu().numpy(); cv = cpu.volumes()
            gd = sim.matrices()[1].cpu().numpy(); cd = cpu.signaling_graph(cpu.contacts()[0]).delta
            errors = dict(phase_abs_max=float(abs(gp-cpu.phi).max()),
                effective_phase_abs_max=float(abs(gp.astype(float)+gc-cpu.phi.astype(float)-cpu.phase_carry).max()),
                chemical_log_max=float(abs(np.log(gu/cu)).max()),
                polarity_abs_max=float(abs(sim.polarity.cpu().numpy()-cpu.polarity).max()),
                relative_volume_max=float(abs(gv/cv-1).max()),
                relative_transport_max=float(np.linalg.norm(gd-cd)/(np.linalg.norm(cd) or 1.)))
            for key,value in errors.items(): maxima[key]=max(maxima[key],value)
            if not all(value<=p['prefix_limits'][key] for key,value in maxima.items()):
                raise FloatingPointError('Actual-start full-state reference gate failed')
        sim.checkpoint(folder/'prefix-gpu.npz'); save_cpu(cpu,folder/'prefix-cpu.npz')
    except Exception as exc:
        failure=dict(type=type(exc).__name__,message=str(exc),time=sim.time)
        sim.checkpoint(folder/'prefix-failure-state.npz')
    result=dict(passed=failure is None, case=job['key'], failure=failure,
        steps=p['prefix_steps'], errors=maxima, exact_gpu_restart=failure is None,
        seconds=time.perf_counter()-started, protocol_sha256=digest(Path(root)/'protocol.json'),
        evidence_sha256=pin(folder.glob('prefix-*.npz')))
    write(output,result); del sim,cpu,clone; torch.cuda.empty_cache()
    return result


def assess_pair(root, pair, p):
    cases=pair['cases']; output=Path(root)/f"time-h{pair['seed']}-{pair['branch']}.json"
    if output.exists(): raise FileExistsError('Never overwrite a recorded timestep comparison')
    reports=[controls.case_report(root,next(j for j in p['jobs'] if j['key']==case),p) for case in cases]
    temporal=q.comparison(root,*cases,p['time_limits'])
    dense=space.dense_pair(root,*cases) if all(r['passed'] for r in reports) else None
    result=dict(pair, passed=bool(all(r['passed'] for r in reports) and temporal['passed']
        and dense is not None and dense['peak']['relative_spectrum']<=p['dense_time_spectrum_limit']),
        temporal=temporal,dense=dense,protocol_sha256=digest(Path(root)/'protocol.json'),
        evidence_sha256=pin(Path(root)/case/'gpu-result.json' for case in cases))
    write(output,result); return result


def complete_design(pairs):
    expected={(seed,branch) for seed in HISTORIES for branch in BRANCHES}
    return (len(pairs)==len(expected) and {(r['seed'],r['branch']) for r in pairs}==expected
        and all(r['passed'] for r in pairs))


def validate_prelaunch(root,p):
    verify_preserved(p); verify_seed_inventory(read(Path(root)/'seed-inventory.json'))
    if p['jobs'] != [dict(j,source=str(Path(root)/j['key']/'source.npz')) for j in design(p['baseline_config'])]:
        raise ValueError('Frozen factorial design drifted')
    if match_starts(p['jobs']) != p['matched_start_signatures']:
        raise ValueError('Matched initial states drifted')
    gate=read(p['instrumentation_gate'])
    if not gate['passed'] or gate['protocol_sha256'] != p['instrumentation_protocol_sha256']:
        raise ValueError('Unchanged developmental instrumentation required')
    for name,expected in gate['evidence_sha256'].items():
        if digest(name)!=expected: raise ValueError('Instrumentation evidence changed')
    if os.environ.get('CUDA_VISIBLE_DEVICES') != p['gpu_uuid']:
        raise ValueError('Explicit dedicated GPU required')
    if read(Path(root)/'status.json') != dict(state='prepared',protocol_sha256=digest(Path(root)/'protocol.json')):
        raise ValueError('Unchanged prepared status required')


def run(root):
    root=Path(root).resolve()
    with (root/'coordinator.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        p=read(root/'protocol.json'); validate_prelaunch(root,p)
        torch.set_num_threads(1); lib=td.gpu.library()
        if digest(lib.binary_path)!=p['runtime_binary_sha256']:
            raise ValueError('Qualified developmental runtime binary changed')
        write(root/'runtime.json',dict(gpu=torch.cuda.get_device_name(0),torch=torch.__version__,
            binary=lib.binary_path,binary_sha256=digest(lib.binary_path)))
        started=time.perf_counter(); paths=[]; pairs=[]; prefixes=[]
        try:
            for pair in p['time_pairs']:
                for case in pair['cases']:
                    job=next(j for j in p['jobs'] if j['key']==case)
                    write(root/'status.json',dict(state='running_prefix',active=case,completed=len(paths),
                        total=24,independent_histories=3,protocol_sha256=digest(root/'protocol.json'),
                        seconds=time.perf_counter()-started))
                    record=prefix(root,job,p); prefixes.append(record)
                    if not record['passed']: raise RuntimeError('Actual-start gate failed: '+case)
                    write(root/'status.json',dict(state='running_development',active=case,completed=len(paths),
                        total=24,independent_histories=3,protocol_sha256=digest(root/'protocol.json'),
                        seconds=time.perf_counter()-started))
                    result=td.run_path(root,job); paths.append(dict(case=case,passed=result['passed'],seconds=result['seconds']))
                    if not result['passed']: raise RuntimeError('Developmental path failed: '+case)
                result=assess_pair(root,pair,p); pairs.append(result)
                if not result['passed']: raise RuntimeError('Developmental timestep gate failed: '+str(pair))
            verify_preserved(p)
            accepted=complete_design(pairs)
            result=dict(execution_completed=True,developmental_screens_passed=accepted,
                qualification_accepted=False,independent_review_pending=True,
                post_cleavage_started=False,independent_histories=3,new_scientific_histories=3,
                new_numerical_paths=len(paths),paths=paths,pairs=pairs,prefixes=prefixes,
                seconds=time.perf_counter()-started,protocol_sha256=digest(root/'protocol.json'),scope=p['scope'])
            write(root/'summary.json',result)
            write(root/'status.json',dict(state='completed_provisional',completed=len(paths),total=24,
                independent_histories=3,developmental_screens_passed=accepted,
                independent_review_pending=True,post_cleavage_started=False,seconds=result['seconds'],
                protocol_sha256=digest(root/'protocol.json')))
            return result
        except Exception as exc:
            write(root/'status.json',dict(state='failed',completed=len(paths),total=24,
                independent_histories=3,seconds=time.perf_counter()-started,
                failure=dict(type=type(exc).__name__,message=str(exc)),paths=paths,pairs=pairs,
                post_cleavage_started=False,protocol_sha256=digest(root/'protocol.json')))
            raise


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=('prepare','run'))
    parser.add_argument('root',type=Path)
    parser.add_argument('--seed-inventory',type=Path)
    args=parser.parse_args()
    if args.action=='prepare':
        if args.seed_inventory is None: parser.error('prepare requires --seed-inventory')
        prepare(args.root,args.seed_inventory)
    else: run(args.root)


if __name__=='__main__': main()
