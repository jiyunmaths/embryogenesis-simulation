"""One-shot local completion watcher for the explicitly scheduled CPU review."""
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
STUDY = REPO / 'output/studies/development-control-qualification'
REPLAY_CASES = ('off120', 'tension_adhesion_space120', 'polarity_space120', 'full_space120')


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''): digest.update(block)
    return digest.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def write_status(state, **details):
    target = HERE / 'schedule-status.json'; temporary = target.with_suffix('.tmp')
    temporary.write_text(json.dumps(dict(state=state, utc=datetime.now(timezone.utc).isoformat(),
                                        **details), indent=2) + '\n')
    temporary.replace(target)


def process_identity(pid):
    try:
        content = Path(f'/proc/{pid}/stat').read_text()
    except FileNotFoundError:
        return None
    fields = content[content.rfind(')') + 2:].split()
    return dict(state=fields[0], start_ticks=fields[19])


def readiness(status, summary, coordinator_alive, expected_protocol):
    """Never start on partial evidence, before process exit, or after failure."""
    if status.get('state') == 'failed': return 'failed', 'Simulation coordinator reported failure'
    if coordinator_alive: return 'waiting', 'Coordinator is still running'
    if status.get('state') not in ('completed_provisional', 'completed_with_unresolved_checks'):
        return 'failed', 'Coordinator exited before a terminal completed status'
    if summary is None: return 'failed', 'Completed status has no complete summary'
    required = {(case, stage) for case in REPLAY_CASES for stage in (2, 4, 8)}
    if not summary.get('execution_completed') or summary.get('protocol_sha256') != expected_protocol:
        return 'failed', 'Summary does not match the scheduled completed study'
    if summary.get('new_scientific_histories') != 0:
        return 'failed', 'Summary scientific-history scope changed'
    if len(summary.get('paths', [])) != 13 or len(summary.get('branch_rows', [])) != 3:
        return 'failed', 'Summary omits required paths or branches'
    rows = summary.get('replay_rows', [])
    if len(rows) != 12 or {(x['case'], x['stage']) for x in rows} != required:
        return 'failed', 'Summary omits or duplicates required replay cases'
    if status.get('replays_completed') != 12:
        return 'failed', 'Terminal status does not record twelve completed replays'
    return 'ready', 'Coordinator exited with all required evidence recorded'


def main():
    with (HERE / 'schedule.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        request = read(HERE / 'schedule-request.json')
        if (HERE / 'independent_review.json').exists():
            write_status('already_reviewed', report=str(HERE / 'independent_review.json'))
            return
        for file, expected in request['review_source_sha256'].items():
            if sha(file) != expected: raise ValueError('Scheduled review source changed: ' + file)
        deadline = time.monotonic() + request['timeout_seconds']
        while True:
            if sha(STUDY / 'protocol.json') != request['protocol_sha256']:
                raise ValueError('Scheduled study protocol changed')
            identity = process_identity(request['coordinator_pid'])
            alive = identity is not None and identity['state'] != 'Z' and identity['start_ticks'] == request['coordinator_start_ticks']
            try:
                status = read(STUDY / 'status.json')
                summary = read(STUDY / 'summary.json') if (STUDY / 'summary.json').exists() else None
            except json.JSONDecodeError:
                if not alive: raise ValueError('Incomplete JSON evidence after coordinator exit')
                status = dict(state='running'); summary = None
            state, reason = readiness(status, summary, alive, request['protocol_sha256'])
            if state == 'failed':
                write_status('not_started', reason=reason, study_status=status)
                return
            if state == 'ready': break
            write_status('waiting_for_completion', reason=reason, coordinator_pid=request['coordinator_pid'],
                         study_status=status, poll_seconds=request['poll_seconds'])
            if time.monotonic() >= deadline:
                write_status('timed_out', reason='Coordinator did not complete within twelve hours')
                return
            time.sleep(request['poll_seconds'])
        for file, expected in request['review_source_sha256'].items():
            if sha(file) != expected: raise ValueError('Scheduled review source changed: ' + file)
        if not (STUDY / 'timing.json').exists(): raise ValueError('Completed study timing is missing')
        write_status('reviewing', argv=request['review_argv'], study_status=status)
        env = dict(os.environ, CUDA_VISIBLE_DEVICES='', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
                   MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1', PYTHONUNBUFFERED='1')
        with (HERE / 'review.log').open('x') as log:
            completed = subprocess.run(request['review_argv'], cwd=REPO, env=env,
                                       stdout=log, stderr=subprocess.STDOUT)
        if completed.returncode != 0:
            write_status('review_failed', exit_code=completed.returncode, log=str(HERE / 'review.log'))
            return
        report = read(HERE / 'independent_review.json')
        if not report['review_completed'] or report['protocol_sha256'] != request['protocol_sha256']:
            raise ValueError('Reviewer did not produce the required matching report')
        write_status('completed', qualification_accepted=report['qualification_accepted'],
                     review_seconds=report['review_seconds'], report=str(HERE / 'independent_review.json'),
                     original_reports_unchanged=report['original_reports_unchanged'])


if __name__ == '__main__':
    try:
        main()
    except BlockingIOError:
        raise SystemExit('A completion watcher already holds this schedule lock')
    except Exception as exc:
        write_status('schedule_failed', error=dict(type=type(exc).__name__, message=str(exc)))
        raise
