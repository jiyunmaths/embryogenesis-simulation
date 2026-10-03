"""Audited parallel native-backend continuation of the moving response study."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import fcntl
import json
import multiprocessing as mp
from pathlib import Path
from . import fast_response_resume as adapter
from . import cell_exchange_response_moving as study
from .native_mechanics import NativeSimulation
from .feedback_long import digest
from .resolution import write_json


def worker(args):
    root, job, threads = args
    # Each spawned process owns its module globals. Reuse the established
    # checkpoint/diagnostic loop, substituting only the simulation factory.
    class ConfiguredSimulation(NativeSimulation):
        native_threads = threads
    previous = adapter.FastAttributeSimulation
    adapter.FastAttributeSimulation = ConfiguredSimulation
    try:
        return adapter.worker((root, job))
    finally:
        adapter.FastAttributeSimulation = previous


def jobs(root, protocol):
    rows = []
    for background in protocol['backgrounds']:
        if not study.source_ready(protocol['sources'][background]):
            raise ValueError('Source not ready: '+background)
        study.materialize(root, protocol, background)
        child = root/background
        p = json.loads((child/'protocol.json').read_text())
        study.verify(p)
        for job in p['jobs']:
            folder = child/adapter.name(job)
            completed = (folder/'result.json').exists()
            if completed:
                result = json.loads((folder/'result.json').read_text())
                if not result['quality_pass'] or result['protocol_sha256'] != digest(child/'protocol.json'):
                    raise ValueError('Invalid saved result: '+str(folder))
            rows.append((child, job, folder, completed))
    return rows


def prepare(root, validation, workers=3, threads=2):
    root, validation = Path(root).resolve(), Path(validation).resolve()
    if not 1 <= workers <= 6 or not 1 <= threads <= 6 or workers*threads > 6:
        raise ValueError('Use at most six native threads across workers')
    target = root/'native-backend-transition.json'
    if target.exists():
        raise FileExistsError(target)
    vp = json.loads((validation/'protocol.json').read_text())
    study.verify(vp)
    if not json.loads((validation/'comparison.json').read_text())['passed']:
        raise ValueError('Native validation did not pass')
    evidence = [validation/'protocol.json', validation/'comparison.json']
    for job in vp['jobs']:
        f = validation/job['key']/'result.json'
        result = json.loads(f.read_text())
        if not result['passed'] or result['protocol_sha256'] != digest(validation/'protocol.json'):
            raise ValueError('Invalid validation result')
        evidence.append(f)
    protocol = json.loads((root/'protocol.json').read_text())
    study.verify(protocol)
    old_transition = root/'backend-transition.json'
    old = json.loads(old_transition.read_text())
    for f, h in {**old['source_sha256'], **old['validation_sha256']}.items():
        if digest(f) != h:
            raise ValueError('Previous backend provenance changed: '+f)
    rows = jobs(root, protocol)
    code = [Path(__file__), Path(adapter.__file__), Path(study.__file__)]
    manifest = dict(created_utc=datetime.now(timezone.utc).isoformat(), workers=workers,
        threads_per_worker=threads, original_protocol_sha256=digest(root/'protocol.json'),
        previous_transition_sha256=digest(old_transition),
        source_sha256={**vp['source_sha256'], **{str(f.resolve()): digest(f) for f in code}},
        validation_sha256={str(f): digest(f) for f in evidence},
        completed_results={str(folder/'result.json'): digest(folder/'result.json')
                           for _, _, folder, done in rows if done},
        pending=[dict(child=str(child), job=job, folder=str(folder),
                      checkpoint_sha256=digest(folder/'latest_state.npz') if (folder/'latest_state.npz').exists() else None)
                 for child, job, folder, done in rows if not done],
        scope='Execution-only transition; original scientific protocols, accepted outputs, timestep and checkpoints preserved. Full control/pulse validation passed with four threads; per-cell arithmetic is independent of thread scheduling. Three spawned workers use two native threads each, BLAS one thread.')
    write_json(target, manifest)
    return manifest


def run(root):
    root = Path(root).resolve()
    with (root/'native-run.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        manifest_path = root/'native-backend-transition.json'
        m = json.loads(manifest_path.read_text())
        for f, h in {**m['source_sha256'], **m['validation_sha256'], **m['completed_results']}.items():
            if digest(f) != h:
                raise ValueError('Transition evidence changed: '+f)
        if digest(root/'protocol.json') != m['original_protocol_sha256'] or digest(root/'backend-transition.json') != m['previous_transition_sha256']:
            raise ValueError('Original protocol or prior transition changed')
        p = json.loads((root/'protocol.json').read_text()); study.verify(p)
        rows = jobs(root, p)
        done = sum(row[3] for row in rows)
        pending = [row for row in rows if not row[3]]
        def status(state, **extra):
            write_json(root/'status.json', dict(state=state, completed=done, total=len(rows),
                workers=m['workers'], threads_per_worker=m['threads_per_worker'], backend='native',
                native_transition_sha256=digest(manifest_path), **extra))
        try:
            for child, job, folder, _ in pending:
                folder.mkdir(exist_ok=True)
                execution = folder/'native-backend-execution.json'
                if not execution.exists():
                    expected = next(x for x in m['pending'] if x['folder'] == str(folder))
                    ck = folder/'latest_state.npz'
                    actual = digest(ck) if ck.exists() else None
                    if actual != expected['checkpoint_sha256']:
                        raise ValueError('Resume checkpoint changed before migration')
                    write_json(execution, dict(native_transition_sha256=digest(manifest_path),
                        resume_checkpoint_sha256=actual, job=job))
                elif json.loads(execution.read_text())['native_transition_sha256'] != digest(manifest_path):
                    raise ValueError('Execution manifest mismatch')
            status('running', pending=[str(row[2].relative_to(root)) for row in pending])
            with ProcessPoolExecutor(max_workers=m['workers'], mp_context=mp.get_context('spawn')) as pool:
                futures = {pool.submit(worker, (str(child), job, m['threads_per_worker'])): folder
                           for child, job, folder, _ in pending}
                for future in as_completed(futures):
                    future.result(); done += 1
                    status('running', last_completed=str(futures[future].relative_to(root)))
            study.assess(root, p)
            status('completed')
        except Exception as exc:
            status('failed', error=str(exc)); raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('prepare', 'run'))
    parser.add_argument('--output', type=Path, default=Path('outputs/cell-exchange-response-moving'))
    parser.add_argument('--validation', type=Path, default=Path('outputs/native-mechanics-validation'))
    args = parser.parse_args()
    prepare(args.output, args.validation) if args.command == 'prepare' else run(args.output)
