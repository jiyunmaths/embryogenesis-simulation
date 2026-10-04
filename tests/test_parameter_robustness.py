import numpy as np
import pytest

from embryo import parameter_robustness as study
from embryo.feedback_endpoint_bistability import chemical_jacobian
from embryo.resolution import write_json
from embryo.feedback_long import digest


def graph():
    masses = np.array([1., 2., 3.])
    conductance = np.array([[0., 1., 2.], [1., 0., 3.], [2., 3., 0.]])
    stiffness = np.diag(conductance.sum(axis=1))-conductance
    return -stiffness/masses[:, None], masses


@pytest.mark.parametrize('ratio', (8., 20., 40.))
def test_mass_symmetric_modal_spectrum_matches_full_chemical_jacobian(ratio):
    delta, masses = graph()
    result = study.frozen_spectrum(delta, masses, 2., .02, ratio)
    jac = chemical_jacobian(np.ones((2, 3)), delta, 2., .02, .02*ratio)
    assert result['uniform_jacobian_max_real'] == pytest.approx(np.linalg.eigvals(jac).real.max(), abs=1e-13)
    np.testing.assert_allclose(sorted(result['lambdas']), sorted(np.linalg.eigvals(-delta).real), atol=1e-14)
    assert study.validate_parameters(2., .02, ratio) == .02*ratio


def test_frozen_spectrum_has_no_mechanics_parameter_and_uniform_threshold_is_graph_specific():
    import inspect
    assert 'chi' not in inspect.signature(study.frozen_spectrum).parameters
    delta, masses = graph()
    threshold = study.uniform_threshold(delta, masses, 2., .02, upper=80.)
    assert threshold is not None
    assert abs(study.frozen_spectrum(delta, masses, 2., .02, threshold)['uniform_jacobian_max_real']) < 1e-10
    assert study.frozen_spectrum(delta, masses, 2., .02, threshold-.1)['uniform_jacobian_max_real'] < 0
    assert study.frozen_spectrum(delta, masses, 2., .02, threshold+.1)['uniform_jacobian_max_real'] > 0


@pytest.mark.parametrize('beta,da,ratio', [(1., .02, 20.), (2., 0., 20.), (2., .02, np.nan), (2., .02, -1.)])
def test_invalid_parameters_rejected(beta, da, ratio):
    with pytest.raises(ValueError): study.validate_parameters(beta, da, ratio)


def trial(kind, error=0.):
    return dict(endpoint_kind=kind, solver_pass=True, pattern_return_log_rms=error)


def test_local_coexistence_separates_initiation_from_maintenance():
    trials = [trial('patterned') for _ in range(3)]+[trial('uniform') for _ in range(2)]
    r = study.classify(trials, -.1, study.CRITERIA)
    assert r['local_bistability_supported'] and r['maintenance_supported']
    assert not r['initiation_supported'] and r['phase'] == 'coexistence'
    r = study.classify([trial('patterned') for _ in range(5)], .1, study.CRITERIA)
    assert r['initiation_supported'] and r['maintenance_supported'] and not r['local_bistability_supported']


def test_failed_solver_or_unsettled_state_never_becomes_a_negative_attractor_claim():
    trials = [trial('uniform') for _ in range(5)]
    trials[0]['solver_pass'] = False
    assert study.classify(trials, -.1, study.CRITERIA)['phase'] == 'unresolved'
    trials[0] = trial('unresolved')
    assert study.classify(trials, -.1, study.CRITERIA)['phase'] == 'unresolved'
    trials[0] = trial('uniform')
    assert study.classify(trials, -.1, study.CRITERIA)['phase'] == 'uniform_from_sampled_starts'


def test_different_pattern_endpoints_do_not_count_as_local_return():
    trials = [trial('patterned') for _ in range(5)]
    trials[1]['pattern_return_log_rms'] = .2
    r = study.classify(trials, -.1, study.CRITERIA)
    assert not r['local_bistability_supported'] and not r['maintenance_supported']


def test_worker_real_dual_solver_and_hash_checked_resume(tmp_path):
    delta, masses = graph()
    source = tmp_path/'source.npz'
    pattern = np.ones((2, 3))
    np.savez(source, delta=delta, volumes=masses, trajectories=pattern[None,None], ids=np.arange(3))
    p = dict(noise_seeds=[0, 1], pattern_log_noise=.01, uniform_log_noise=.001,
             horizons=[30., 120.], interval=2., criteria=study.CRITERIA)
    write_json(tmp_path/'protocol.json', p)
    job = dict(job_key='fixture', source=str(source), seed=7, branch='switch_on', ratio=8.,
               beta=2., da=.02, db=.16)
    result = study.worker((tmp_path, job, p))
    assert result['numerical_pass'] and result['phase'] == 'uniform_from_sampled_starts'
    assert len(result['trials']) == 5 and all(t['solver_log_error'] < 1e-8 for t in result['trials'])
    assert study.worker((tmp_path, job, p)) == result
    with np.load(tmp_path/'fixture'/'paths.npz') as z:
        np.testing.assert_allclose(z['initial_states'][1:]@masses, np.tile(pattern@masses, (4, 1)), atol=1e-14)
    with (tmp_path/'fixture'/'paths.npz').open('ab') as f: f.write(b'changed')
    with pytest.raises(ValueError, match='Changed'): study.worker((tmp_path, job, p))


def test_verify_rejects_dependency_change(tmp_path):
    f = tmp_path/'dependency'; f.write_text('old')
    p = dict(source_sha256={}, input_sha256={str(f):digest(f)})
    study.verify(p); f.write_text('new')
    with pytest.raises(ValueError, match='Changed'): study.verify(p)
