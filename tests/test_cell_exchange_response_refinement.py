import json
import numpy as np
import pytest

from embryo.cell_exchange_response_refinement import assess, BACKGROUNDS
from embryo.cell_response_moving import name
from embryo.cell_response_moving_refinement import CRITERIA
from embryo.feedback_long import digest


def fixture_study(tmp_path):
    baseline, fine = tmp_path / 'baseline', tmp_path / 'fine'
    times = np.arange(0., 61., 3.)
    targets = [27, 20]
    protocol = dict(baseline=str(baseline), backgrounds=BACKGROUNDS,
                    targets=targets, factor=.9, duration=60., interval=3.,
                    refinement_criteria=CRITERIA,
                    classification_criteria=dict(reference_separation_min=.01),
                    scope='Synthetic matched continuations', limits='Test fixture',
                    source_sha256={}, input_sha256={})
    for root in (baseline, fine):
        root.mkdir()
        (root / 'protocol.json').write_text(json.dumps(protocol))
        for background in BACKGROUNDS:
            child = root / background
            child.mkdir()
            jobs = [dict(family=background, target=None, factor=1.)] + [
                dict(family=background, target=cell, factor=.9) for cell in targets]
            cp = dict(jobs=jobs, source_sha256={}, input_sha256={})
            (child / 'protocol.json').write_text(json.dumps(cp))
            np.savez(child / 'initial_states.npz', ids=targets, masses=[1., 2.])
            for job in jobs:
                folder = child / name(job)
                folder.mkdir()
                (folder / 'result.json').write_text(json.dumps(dict(
                    quality_pass=True, protocol_sha256=digest(child / 'protocol.json'))))
                path = np.ones((len(times), 2, 2))
                if job['target'] is not None:
                    index = targets.index(job['target'])
                    reference = index if background == 'unexchanged' else 1-index
                    tau = (1., 5.)[reference]
                    path[:, 0, index] = np.exp(np.log(.9)*np.exp(-times/tau))
                    path[:, 1, index] = np.exp(np.log(.9)*.2*(times/tau)*np.exp(-times/tau))
                history = [dict(elapsed=float(t), ids=targets, chemistry=x.tolist())
                           for t, x in zip(times, path)]
                (folder / 'history.json').write_text(json.dumps(history))
    return baseline, fine


def pulse_folder(root, background='fresh_exchange', cell=27):
    return root / background / name(dict(family=background, target=cell, factor=.9))


def test_matching_responses_and_transfer_pass(tmp_path):
    _, fine = fixture_study(tmp_path)
    report = assess(fine)
    assert report['passed']
    assert report['max_raw_chemical_log_error'] == 0
    assert len(report['trials']) == 4
    assert all(row['fine']['nearest_reference'] == 'donor' for row in report['classifications'])


def test_destination_response_replacement_detects_classification_flip(tmp_path):
    _, fine = fixture_study(tmp_path)
    source = pulse_folder(fine, 'unexchanged') / 'history.json'
    (pulse_folder(fine) / 'history.json').write_text(source.read_text())
    report = assess(fine)
    assert not report['passed']
    row = report['classifications'][0]
    assert row['coarse']['nearest_reference'] == 'donor'
    assert row['fine']['nearest_reference'] == 'destination'
    assert not row['passed']


def test_lost_recovery_is_not_treated_as_time_agreement(tmp_path):
    _, fine = fixture_study(tmp_path)
    path = pulse_folder(fine) / 'history.json'
    history = json.loads(path.read_text())
    for row in history[1:]:
        row['chemistry'][0][0] = float(np.exp(.2*np.log(.9)))
    path.write_text(json.dumps(history))
    report = assess(fine)
    row = next(r for r in report['trials'] if r['background'] == 'fresh_exchange' and r['cell'] == 27)
    assert not report['passed']
    assert row['coarse_metrics']['target_recovery_time'] is not None
    assert row['fine_metrics']['target_recovery_time'] is None
    assert row['recovery_time_errors']['target_recovery_time'] is None


@pytest.mark.parametrize('error', ['ids', 'times', 'quality', 'protocol', 'missing'])
def test_invalid_evidence_is_rejected(tmp_path, error):
    _, fine = fixture_study(tmp_path)
    folder = pulse_folder(fine)
    if error in ('ids', 'times'):
        path = folder / 'history.json'
        history = json.loads(path.read_text())
        history[1]['ids' if error == 'ids' else 'elapsed'] = [20, 27] if error == 'ids' else 3.1
        path.write_text(json.dumps(history))
    elif error == 'missing':
        (folder / 'result.json').unlink()
    else:
        path = folder / 'result.json'
        result = json.loads(path.read_text())
        result['quality_pass' if error == 'quality' else 'protocol_sha256'] = False if error == 'quality' else 'wrong'
        path.write_text(json.dumps(result))
    with pytest.raises((ValueError, FileNotFoundError)):
        assess(fine)
