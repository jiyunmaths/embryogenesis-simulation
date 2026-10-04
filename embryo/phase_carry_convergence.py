"""Third carry timestep plus a separate independent float64 mechanics check.

Reuse the completed fine/finer formation evidence read-only. Add just one
coarse carry trajectory, preserving all original failures and thresholds.
"""
import argparse
from dataclasses import asdict
import fcntl
from pathlib import Path

import numpy as np
import scipy
import torch

from .attribute_development import AttributeSimulation
from .feedback_long import digest, restore_checkpoint
from .feedback_survival_validation import retime
from .geometry_precision import checked_history, implementation_gate, worker
from .gpu_backend import library as accepted_library
from .gpu_precision_control import library
from .mechanics_float64_reference import run as reference_run, CRITERIA, DTS, HORIZON
from .neighbor_context import read
from .parameter_robustness import verify
from .parameter_robustness_moving import endpoint_assay
from .polarity_robustness import refinement_comparison
from .polarity_robustness_refinement import check_halving
from .resolution import write_json


def completed_evidence(root, p, job):
    folder = Path(root)/job['key']; ph = digest(Path(root)/'protocol.json')
    r, h = read(folder/'result.json'), read(folder/'history.json')
    if (r['job'] != job or r['protocol_sha256'] != ph or not r['quality_pass'] or
            r['history_sha256'] != digest(folder/'history.json') or
            r['checkpoint_sha256'] != digest(folder/'latest_state.npz')):
        raise ValueError('Changed completed moving evidence')
    checked_history(h, p, job, p['duration'])
    host, audit, saved = restore_checkpoint(folder/'latest_state.npz', ph, job)
    if (saved != h or audit != r['audit'] or host.time != p['start']+p['duration'] or
            not np.array_equal([host.activator, host.inhibitor], h[-1]['chemistry'])):
        raise ValueError('Moving checkpoint and history disagree')
    with np.load(folder/'latest_state.npz') as z:
        if (str(z['precision_arm']) != 'phase_carry' or z['phase_carry'].dtype != np.float64 or
                z['phase_carry'].shape != host.phi.shape or not np.isfinite(z['phase_carry']).all()):
            raise ValueError('Invalid completed carry state')
    return r, h


def prepare(root, parent=Path('outputs/phase-carry-formation')):
    root, parent = Path(root).resolve(), Path(parent).resolve()
    if root.exists(): raise FileExistsError(root)
    p = read(parent/'protocol.json'); verify(p)
    summary, status = read(parent/'summary.json'), read(parent/'status.json')
    if (status['state'] != 'completed' or not summary['passed'] or
            summary['protocol_sha256'] != digest(parent/'protocol.json')):
        raise ValueError('Requires the completed accepted carry pair')
    references = {}; inputs = [parent/f for f in ('protocol.json', 'summary.json', 'status.json')]
    for j in p['jobs']:
        completed_evidence(parent, p, j)
        references[j['level']] = dict(folder=str(parent/j['key']), job=j)
        folder = parent/j['key']
        ep = read(folder/'endpoint/protocol.json'); verify(ep)
        er = read(folder/'endpoint/assay/result.json')
        if (not er['numerical_pass'] or not er['all_trials_settled'] or
                er['protocol_sha256'] != digest(folder/'endpoint/protocol.json') or
                er['paths_sha256'] != digest(folder/'endpoint/assay/paths.npz') or
                ep['moving_history_sha256'] != digest(folder/'history.json') or
                ep['moving_protocol_sha256'] != digest(parent/'protocol.json')):
            raise ValueError('Changed completed frozen endpoint')
        inputs += [folder/f for f in ('result.json', 'history.json', 'latest_state.npz',
            'endpoint/protocol.json', 'endpoint/assay/result.json', 'endpoint/assay/paths.npz')]
    root.mkdir(parents=True)
    original_job = references['fine']['job']
    source = root/'source.npz'
    coarse = retime(original_job['source'], source, DTS[0])
    original = AttributeSimulation.restore(original_job['source'])
    check_halving(coarse, original)  # coarse-to-fine direction
    job = dict(original_job, key='phase_carry_coarse', level='coarse', dt=DTS[0], source=str(source))
    config = asdict(coarse.config)
    torch.set_num_threads(1)
    gate = implementation_gate(root, source, Path(job['chemical_file']), job['point'], config)
    device = dict(name=torch.cuda.get_device_name(0), torch=torch.__version__, torch_cuda=torch.version.cuda,
        experimental_binary=str(library()._name), experimental_binary_sha256=digest(library()._name),
        accepted_binary=str(accepted_library()._name), accepted_binary_sha256=digest(accepted_library()._name))
    write_json(root/'implementation-gate.json', dict(**gate, device=device))
    contexts = [dict(key='mature-uniform', source=str(source), chemical_file=job['chemical_file']),
        dict(key='patterned-endpoint', source=str(Path(references['finer']['folder'])/'latest_state.npz'))]
    inputs += [source, root/'implementation-gate.json', root/'gate-checkpoint.npz',
               Path(device['experimental_binary']), Path(device['accepted_binary'])]
    new_sources = [Path(__file__), Path(__file__).with_name('mechanics_float64_reference.py'),
                   Path(__file__).with_name('phase_carry_formation.py')]
    sources = {**p['source_sha256'], **{str(f.resolve()):digest(f) for f in new_sources}}
    protocol = dict(p, jobs=[job], accepted_configs={'coarse':config}, parent=str(parent),
        reused_references=references, dts=list(DTS), reuse_duration=0., device=device,
        independent_histories=1, new_histories=0, existing_history=9,
        reference_contexts=contexts, reference_horizon=HORIZON, reference_criteria=CRITERIA,
        reference_runtime=dict(numpy=np.__version__, scipy=scipy.__version__, torch_threads=1),
        trend_criteria=dict(chemical_pair_error_ratio_min=1.5),
        estimated_moving_seconds=summary['results'][0]['audit']['wall_seconds']/2,
        design='One new dt=.00375 carry run, starting from exactly the same mature t=150 history-9 geometry, RNG/lineage/polarity and near-uniform chemistry as the completed dt=.001875/.0009375 pair. Compare all three through elapsed 240. Same 72^3 grid, beta=2, D_a=.02, D_b=.55, chi=0, activity tension/adhesion .25/.35. No new histories; no cleavage.',
        reference_design='BEFORE the new moving run, independently evaluate float64 NumPy/SciPy mechanics on the shared initial geometry and actual finest patterned endpoint. Six .15-unit mechanics-only runs at three dt values; freeze chemical concentrations and polarity, but recompute evolving geometry/volume forces/overlaps. Compare baseline/carry CUDA against the CPU fields every step. All component-assay residuals start at zero. This is separate from co-evolving moving comparisons.',
        interpretation='Require both adjacent moving pairs to pass the unchanged original raw tolerances, plus >=1.5 reduction in maximum chemical pair error under halving, all physical quality checks, and dual-solver endpoint stationarity. A decreasing three-point error trend and short independent mechanics check do not establish full-system high-precision equivalence, spatial convergence, history robustness, positive-chi behavior, initiation from zygotes or autonomous cell identity. Keep old no-carry failures unchanged.',
        source_sha256=sources, input_sha256={**p['input_sha256'], **{str(f.resolve()):digest(f) for f in inputs}})
    # Drop inherited labels that refer to the old continuation design.
    protocol.pop('estimated_new_seconds', None)
    write_json(root/'protocol.json', protocol)
    write_json(root/'status.json', dict(state='prepared', completed=0, total=1, reused_trajectories=2,
        independent_histories=1, new_histories=0, estimated_moving_seconds=protocol['estimated_moving_seconds']))
    return protocol


def error_trend(coarse_error, fine_error, minimum=1.5):
    if not np.isfinite([coarse_error, fine_error]).all() or min(coarse_error, fine_error) < 0:
        raise ValueError('Finite nonnegative adjacent errors required')
    ratio = coarse_error/fine_error if fine_error > 0 else None
    return dict(coarse_fine=coarse_error, fine_finer=fine_error, ratio=ratio,
        observed_pair_order=None if ratio is None or ratio <= 0 else float(np.log2(ratio)),
        passed=ratio is not None and ratio >= minimum,
        interpretation='Ratio of adjacent-pair maximum errors at common saved times, not proof of asymptotic order or an absolute solution-error estimate.')


def reference_evidence(root, p):
    root = Path(root); gate = read(root/'reference-verification.json')
    result_file = root/'float64-reference/result.json'
    if (gate['protocol_sha256'] != digest(root/'protocol.json') or not gate['passed'] or
            gate['result_sha256'] != digest(result_file)):
        raise ValueError('Changed independent reference gate')
    result = read(result_file)
    if (not result['passed'] or result['contexts'] != p['reference_contexts'] or
            result['criteria'] != p['reference_criteria'] or result['horizon'] != p['reference_horizon'] or
            result['dts'] != p['dts']): raise ValueError('Changed independent reference design')
    for f, h in result['field_sha256'].items():
        if digest(f) != h: raise ValueError('Changed independent reference fields')
    return result


def assess(root):
    root = Path(root); p = read(root/'protocol.json'); verify(p)
    ref = reference_evidence(root, p)
    old = read(Path(p['parent'])/'protocol.json')
    histories = {}; rows = []
    for level, reuse in p['reused_references'].items():
        r, h = completed_evidence(Path(reuse['folder']).parent, old, reuse['job'])
        histories[level] = h; rows.append(dict(level=level, reused=True, result=r))
    job = p['jobs'][0]; r, h = completed_evidence(root, p, job)
    histories['coarse'] = h; rows.append(dict(level='coarse', reused=False, result=r))
    e = endpoint_assay(root, job, p, h[-1])
    comparisons = {key:refinement_comparison(histories[a], histories[b], p, p['duration'])
        for key, a, b in (('coarse_fine', 'coarse', 'fine'), ('fine_finer', 'fine', 'finer'))}
    trend = error_trend(comparisons['coarse_fine']['errors']['chemical_log_max'],
        comparisons['fine_finer']['errors']['chemical_log_max'], p['trend_criteria']['chemical_pair_error_ratio_min'])
    endpoint_pass = bool(e['numerical_pass'] and e['all_trials_settled'])
    summary = dict(protocol_sha256=digest(root/'protocol.json'), completed=1, total=1,
        moving_trajectories=3, reused_trajectories=2, independent_histories=1, new_histories=0,
        reference_pass=ref['passed'], reference=ref, quality_pass=True, results=rows,
        full_refinement=comparisons, chemical_error_trend=trend,
        new_endpoint=dict(numerical_pass=e['numerical_pass'], all_trials_settled=e['all_trials_settled'],
            phase=e['phase'], local_bistability=e['local_bistability_supported']),
        passed=bool(all(c['passed'] for c in comparisons.values()) and trend['passed'] and endpoint_pass),
        scope=p['interpretation'])
    write_json(root/'summary.json', summary); return summary


def run(root):
    root = Path(root).resolve()
    with (root/'coordinator.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        p = read(root/'protocol.json'); verify(p); torch.set_num_threads(1)
        if p['reference_runtime'] != dict(numpy=np.__version__, scipy=scipy.__version__, torch_threads=1):
            raise ValueError('Changed independent reference runtime')
        device = p['device']
        if (device['name'] != torch.cuda.get_device_name(0) or device['torch'] != torch.__version__ or
                device['torch_cuda'] != torch.version.cuda or
                digest(library()._name) != device['experimental_binary_sha256'] or
                digest(accepted_library()._name) != device['accepted_binary_sha256']):
            raise ValueError('Changed GPU hardware/software')
        stage = 'independent_float64_reference'
        try:
            if not (root/'reference-verification.json').exists():
                write_json(root/'status.json', dict(state='running', stage=stage, completed=0, total=1))
                result = reference_run(root/'float64-reference', p['reference_contexts'],
                    p['reference_criteria'], p['reference_horizon'], p['dts'])
                write_json(root/'reference-verification.json', dict(protocol_sha256=digest(root/'protocol.json'),
                    passed=result['passed'], result_sha256=digest(root/'float64-reference/result.json')))
                if not result['passed']: raise RuntimeError('Independent mechanics reference failed; moving run not launched')
            reference_evidence(root, p)
            stage = 'moving_formation'
            write_json(root/'status.json', dict(state='running', stage=stage, completed=0, total=1,
                reference_pass=True, current_job=p['jobs'][0]['key'], reused_trajectories=2))
            worker(root, p['jobs'][0], p)
            stage = 'endpoint_and_three_timestep_assessment'
            write_json(root/'status.json', dict(state='running', stage=stage, completed=1, total=1))
            summary = assess(root)
            write_json(root/'status.json', dict(state='completed' if summary['passed'] else 'completed_with_unresolved_checks',
                completed=1, total=1, passed=summary['passed'], quality_pass=True, reference_pass=True))
            return summary
        except Exception as error:
            write_json(root/'status.json', dict(state='failed', stage=stage, error=str(error))); raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'run', 'assess'))
    parser.add_argument('--output', type=Path, default=Path('outputs/phase-carry-convergence'))
    args = parser.parse_args()
    result = prepare(args.output) if args.action == 'prepare' else run(args.output) if args.action == 'run' else assess(args.output)
    print({k:v for k,v in result.items() if k in ('completed', 'total', 'passed', 'estimated_moving_seconds', 'chemical_error_trend')})
