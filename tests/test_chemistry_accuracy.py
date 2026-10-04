import numpy as np
import pytest

from embryo import chemistry_accuracy as study
from embryo.feedback_long import digest
from embryo.resolution import write_json


def small_protocol(tmp_path):
    delta = np.array([[-3., 3.], [1.5, -1.5]])
    masses = np.array([1., 2.]); initial = np.array([[.3, 1.5], [1.2, .8]])
    file = tmp_path/'context.npz'
    np.savez(file, delta=delta, masses=masses, initial=initial, ids=[7, 9])
    jobs = [dict(key=f'{backend}-{dt}', context='test', backend=backend, dt=dt)
            for backend in ('native_chemistry', 'torch_cpu_chemistry') for dt in (.00375, .001875, .0009375)]
    p = dict(contexts=[dict(key='test', duration=.3, source=str(file))], jobs=jobs,
        dts=[.00375, .001875, .0009375], beta=2., da=.02, db=.55, interval=.15,
        criteria=study.CRITERIA, source_sha256={}, input_sha256={str(file):digest(file)},
        interpretation='Synthetic scientific solver verification')
    write_json(tmp_path/'protocol.json', p)
    return p


def test_fixed_transport_preserves_measured_operator_and_conservative_nullspaces():
    delta = np.array([[-3., 3.], [1.5, -1.5]])
    masses = np.array([1., 2.]); transport = study.fixed_transport(delta, masses)
    np.testing.assert_array_equal(transport.delta.toarray(), delta)
    np.testing.assert_allclose(masses@transport.delta.toarray(), 0, atol=0)
    np.testing.assert_allclose(transport.delta@np.ones(2), 0, atol=0)
    with pytest.raises(ValueError): study.fixed_transport(delta, [1., 1.])


def test_actual_production_functions_converge_to_independent_nonlinear_reference(tmp_path):
    pytest.importorskip('torch')
    p = small_protocol(tmp_path); study.references(tmp_path, p)
    for job in p['jobs']:
        result = study.worker((tmp_path, job))
        assert result['accuracy_pass'] and result['minimum_concentration'] > 0
    report = study.assess(tmp_path)
    assert report['passed'] and report['cpu_torch_log_max'] < 1e-12
    for row in report['comparisons']:
        assert row['order_resolved']
        assert row['orders'] == pytest.approx([2., 2.], abs=.03)
    job = p['jobs'][0]; before = digest(tmp_path/job['key']/'paths.npz')
    study.worker((tmp_path, job))
    assert digest(tmp_path/job['key']/'paths.npz') == before


def test_even_hash_consistent_clock_changes_fail_scientific_alignment(tmp_path):
    p = small_protocol(tmp_path); study.references(tmp_path, p)
    for job in p['jobs']:
        study.worker((tmp_path, job))
    folder = tmp_path/p['jobs'][0]['key']; file = folder/'paths.npz'
    with np.load(file) as z:arrays={k:z[k].copy() for k in z.files}
    arrays['times'][1] += .01; np.savez_compressed(file, **arrays)
    result = study.read(folder/'result.json'); result['paths_sha256'] = digest(file)
    write_json(folder/'result.json', result)
    with pytest.raises(ValueError, match='physical clock'):
        study.assess(tmp_path)
