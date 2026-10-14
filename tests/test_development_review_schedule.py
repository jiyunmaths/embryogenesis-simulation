"""The deferred saved-array review must not consume incomplete evidence."""
import importlib.util
from pathlib import Path

import pytest

PATH = Path(__file__).resolve().parents[1] / 'archive/study-reviews/2026-10-11-development-controls/after_completion.py'
spec = importlib.util.spec_from_file_location('development_review_schedule', PATH)
schedule = importlib.util.module_from_spec(spec); spec.loader.exec_module(schedule)


def complete_summary():
    return dict(execution_completed=True, protocol_sha256='expected', new_scientific_histories=0,
                paths=[{}] * 13, branch_rows=[{}] * 3,
                replay_rows=[dict(case=case, stage=stage)
                             for case in schedule.REPLAY_CASES for stage in (2, 4, 8)])


def test_waits_for_running_coordinator_even_if_summary_exists():
    state, _ = schedule.readiness(dict(state='completed_provisional', replays_completed=12),
                                 complete_summary(), True, 'expected')
    assert state == 'waiting'


@pytest.mark.parametrize('state', ['completed_provisional', 'completed_with_unresolved_checks'])
def test_completed_evidence_can_be_reviewed_after_exit(state):
    assert schedule.readiness(dict(state=state, replays_completed=12), complete_summary(), False, 'expected')[0] == 'ready'


def test_unfinished_or_failed_job_is_never_reviewed():
    assert schedule.readiness(dict(state='running'), None, False, 'expected')[0] == 'failed'
    assert schedule.readiness(dict(state='failed'), None, True, 'expected')[0] == 'failed'


def test_duplicate_replay_cannot_replace_missing_stage():
    summary = complete_summary(); summary['replay_rows'][-1] = summary['replay_rows'][0]
    assert schedule.readiness(dict(state='completed_provisional', replays_completed=12), summary, False, 'expected')[0] == 'failed'


def test_mismatched_protocol_or_missing_summary_is_rejected():
    status = dict(state='completed_provisional', replays_completed=12)
    assert schedule.readiness(status, complete_summary(), False, 'other')[0] == 'failed'
    assert schedule.readiness(status, None, False, 'expected')[0] == 'failed'
