"""Finish an existing serial pilot using independent concurrent control branches.

The model and frozen source protocol are unchanged. Worker protocols select one
arm from that protocol; their reports are combined under the original criteria.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from embryo import moving_causal
from embryo.resolution import write_json


def worker(output):
    # A single-arm worker cannot make between-arm comparisons. Defer those to
    # the parent after collecting all four complete reports.
    def defer(path, protocol):
        return {'all_quality_pass': None}
    moving_causal.compare = defer
    moving_causal.run(output)


def finish(output, serial_pid):
    output = Path(output).resolve()
    protocol_path = output / 'protocol.json'
    protocol = json.loads(protocol_path.read_text())
    digest = hashlib.sha256(protocol_path.read_bytes()).hexdigest()
    initial_status = json.loads((output / 'status.json').read_text())
    if initial_status['state'] != 'running' or initial_status['arm'] != 'full':
        raise ValueError('This handoff requires the serial full branch still running')
    expected = ['python', '-u', '-m', 'embryo.moving_causal', 'run']
    command = Path(f'/proc/{serial_pid}/cmdline').read_bytes().decode().split('\0')
    if command[:5] != expected:
        raise ValueError('Serial process identity mismatch')
    work = output / 'parallel-workers'
    work.mkdir()
    manifest = {'original_protocol_sha256': digest, 'serial_pid': serial_pid,
                'method': 'unchanged run function, separate process per remaining arm; comparisons deferred',
                'workers': {}, 'state': 'running'}
    write_json(work / 'execution.json', manifest)
    processes = []
    for arm in protocol['arms'][1:]:
        path = work / arm
        path.mkdir()
        selected = dict(protocol, arms=[arm], parent_protocol_sha256=digest)
        write_json(path / 'protocol.json', selected)
        write_json(path / 'status.json', {'state': 'prepared', 'protocol_sha256': hashlib.sha256((path / 'protocol.json').read_bytes()).hexdigest()})
        log = (path / 'run.log').open('w')
        process = subprocess.Popen([sys.executable, str(Path(__file__).resolve()), '--worker', str(path)], stdout=log, stderr=subprocess.STDOUT, env=dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1'))
        log.close()
        processes.append((arm, process, path))
        manifest['workers'][arm] = {'pid': process.pid, 'output': str(path)}
    write_json(work / 'execution.json', manifest)
    print('Concurrent control workers started', flush=True)
    # run() writes the full viewer last, before entering the next branch. Wait
    # for the next branch's initial status to ensure the full file is closed.
    while True:
        status = json.loads((output / 'status.json').read_text())
        if status.get('state') == 'failed':
            raise RuntimeError(status)
        if status.get('arm') == 'no_feedback' and (output / 'full' / 'viewer.html').exists():
            break
        time.sleep(5)
    os.kill(serial_pid, signal.SIGTERM)
    while Path(f'/proc/{serial_pid}').exists():
        time.sleep(1)
    # Preserve any redundant serial no-feedback prefix; never mix it with the
    # independent completed worker trajectory.
    write_json(work / 'serial_handoff_status.json', status)
    if (output / 'no_feedback').exists():
        shutil.move(str(output / 'no_feedback'), str(work / 'serial-no-feedback-prefix'))
    write_json(output / 'status.json', {'state': 'running', 'protocol_sha256': digest,
               'completed_arms': ['full'], 'execution': 'parallel controls', 'workers': manifest['workers']})
    print('Full branch retained; redundant serial continuation stopped', flush=True)
    for arm, process, path in processes:
        code = process.wait()
        if code != 0:
            raise RuntimeError(f'{arm} worker failed with exit {code}')
        if not (path / arm / 'final_state.npz').exists():
            raise RuntimeError(f'{arm} has no final checkpoint')
        shutil.copytree(path / arm, output / arm)
        print(f'Collected {arm}', flush=True)
    if hashlib.sha256(protocol_path.read_bytes()).hexdigest() != digest:
        raise RuntimeError('Original protocol changed during execution')
    result = moving_causal.compare(output, protocol)
    write_json(output / 'status.json', {'state': 'completed', 'protocol_sha256': digest,
               'completed_arms': protocol['arms'], 'all_quality_pass': result['all_quality_pass'],
               'execution_manifest': str(work / 'execution.json')})
    manifest['state'] = 'completed'
    write_json(work / 'execution.json', manifest)
    print(json.dumps(result['causal_effect_checks']), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--worker', type=Path)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--serial-pid', type=int)
    args = parser.parse_args()
    if args.worker:
        worker(args.worker)
    else:
        finish(args.output, args.serial_pid)
