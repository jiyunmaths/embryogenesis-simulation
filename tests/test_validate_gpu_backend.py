import json
import numpy as np
import pytest

pytest.importorskip('torch')
from embryo.validate_gpu_backend import assess, discrepancies, CRITERIA
from embryo.cell_response_moving import name
from embryo.feedback_long import digest


def fixture_study(tmp_path):
    baseline, gpu = tmp_path/'cpu', tmp_path/'gpu'
    baseline.mkdir(); gpu.mkdir()
    jobs = [dict(family=f, target=t, factor=1. if t is None else .9)
            for f in ('unexchanged', 'fresh_exchange') for t in (None, 27)]
    protocol = dict(baseline=str(baseline), jobs=jobs, duration=60., interval=3., criteria=CRITERIA,
                    scope='Synthetic', limits='Fixture', source_sha256={}, input_sha256={})
    (gpu/'protocol.json').write_text(json.dumps(protocol))
    times = np.arange(0., 61., 3.)
    for family in ('unexchanged', 'fresh_exchange'):
        (baseline/family).mkdir()
        np.savez(baseline/family/'initial_states.npz', ids=[27, 20], masses=[1., 1.])
    for job in jobs:
        for root in (baseline/job['family'], gpu):
            folder = root/name(job); folder.mkdir()
            tau = 1. if job['family'] == 'unexchanged' else 5.
            path = np.ones((len(times), 2, 2))
            if job['target'] is not None:
                path[:, 0, 0] = np.exp(np.log(.9)*np.exp(-times/tau))
            history = [dict(time=float(t), elapsed=float(t), ids=[27, 20], chemistry=x.tolist(),
                            polarity=[[.1, 0., 0.], [.1, 0., 0.]], axis_ratio=1.3,
                            volumes=[1., 1.], delta=[[-1., 1.], [1., -1.]]) for t, x in zip(times, path)]
            (folder/'history.json').write_text(json.dumps(history))
            np.savez(folder/'latest_state.npz', phi=np.ones((2, 3, 3, 3), np.float32))
            if root == gpu:
                audit = dict(max_volume_error=.01, min_radius=5., max_clipping=0.,
                             max_boundary=0., dilution_amount_error=0.)
                (folder/'result.json').write_text(json.dumps(dict(passed=True,
                    protocol_sha256=digest(gpu/'protocol.json'), audit=audit)))
    return gpu, jobs


def test_full_matching_paths_pass(tmp_path):
    gpu, _ = fixture_study(tmp_path)
    report = assess(gpu)
    assert report['passed'] and report['scientific_ready']
    assert len(report['runs']) == 4 and len(report['responses']) == 2


def test_altered_response_fails_even_if_worker_report_claims_pass(tmp_path):
    gpu, jobs = fixture_study(tmp_path)
    path = gpu/name(jobs[-1])/'history.json'
    history = json.loads(path.read_text()); history[1]['chemistry'][0][0] *= 1.01
    path.write_text(json.dumps(history))
    report = assess(gpu)
    assert not report['passed'] and not report['scientific_ready']
    assert not report['runs'][-1]['passed'] and not report['responses'][-1]['passed']


def test_endpoint_field_is_checked_independently(tmp_path):
    gpu, jobs = fixture_study(tmp_path)
    np.savez(gpu/name(jobs[-1])/'latest_state.npz', phi=np.full((2, 3, 3, 3), .9, np.float32))
    assert not assess(gpu)['passed']


@pytest.mark.parametrize('kind', ['ids', 'time', 'nonfinite', 'nonpositive', 'missing'])
def test_invalid_or_incomplete_observations_rejected(tmp_path, kind):
    gpu, jobs = fixture_study(tmp_path)
    path = gpu/name(jobs[-1])/'history.json'
    history = json.loads(path.read_text())
    if kind == 'ids': history[1]['ids'] = [20, 27]
    elif kind == 'time': history[1]['time'] += .1
    elif kind == 'nonfinite': history[1]['chemistry'][0][0] = float('nan')
    elif kind == 'nonpositive': history[1]['chemistry'][0][0] = 0.
    else: history.pop()
    path.write_text(json.dumps(history))
    with pytest.raises(ValueError): assess(gpu)
