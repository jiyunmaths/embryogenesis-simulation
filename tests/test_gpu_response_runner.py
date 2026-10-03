import json
import os
import numpy as np
import pytest

pytest.importorskip('torch')
from embryo.gpu_response_runner import compatible, require_validation, prepare
from embryo.model import Config


@pytest.mark.parametrize('state,passed,completed', [
    ('prepared', False, 0), ('running', False, 3), ('failed', False, 0), ('completed', False, 4),
    ('completed', True, 2)])
def test_incomplete_or_failed_gate_cannot_prepare_scientific_run(tmp_path, state, passed, completed):
    gate = tmp_path/'gate'; gate.mkdir()
    (gate/'status.json').write_text(json.dumps(dict(state=state, passed=passed, completed=completed)))
    output = tmp_path/'study'
    with pytest.raises(ValueError, match='four-run'):
        prepare(output, tmp_path/'source', gate)
    assert not output.exists()


def test_missing_or_stale_acceptance_report_is_rejected(tmp_path):
    (tmp_path/'status.json').write_text(json.dumps(dict(state='completed', passed=True, completed=4)))
    (tmp_path/'protocol.json').write_text(json.dumps(dict(source_sha256={}, input_sha256={})))
    (tmp_path/'comparison.json').write_text(json.dumps(dict(scientific_ready=True, protocol_sha256='stale')))
    with pytest.raises(ValueError, match='invalid'):
        require_validation(tmp_path)


def test_other_histories_allowed_but_unvalidated_physical_parameters_rejected():
    accepted = Config(grid=72, extent=2.24, dt=.00375)
    other = Config(grid=72, extent=2.24, dt=.00375, seed=8, steps=80000)
    compatible(other, accepted)
    other.signal_dh *= 2
    with pytest.raises(ValueError, match='signal_dh'):
        compatible(other, accepted)
    other.signal_dh = accepted.signal_dh
    other.dt *= .5
    with pytest.raises(ValueError, match='dt'):
        compatible(other, accepted)


@pytest.mark.skipif(os.environ.get('EMBRYO_CUDA_TESTS') != '1', reason='GPU opt-in required')
def test_adapter_preserves_schema_and_resumes_without_reapplying_pulse(tmp_path, monkeypatch):
    from embryo import gpu_response_runner as runner
    from embryo.attribute_development import AttributeSimulation
    from embryo.feedback_long import digest
    from embryo.cell_response_moving import name
    source = tmp_path/'source'; source.mkdir()
    host = AttributeSimulation(Config(grid=24, extent=1.4, interface_width=.12,
                                      max_cells=1, dt=.00375), 'direct')
    host.polarity[:] = [.2, .1, 0.]
    checkpoint = source/'source.npz'; host.checkpoint(checkpoint)
    np.savez(source/'initial_states.npz', state=np.array([host.activator, host.inhibitor]),
             ids=host.ids, masses=host.volumes())
    jobs = [dict(family='state', target=None, factor=1.), dict(family='state', target=0, factor=.9)]
    p = dict(checkpoint=str(checkpoint), start=0., dt=.00375, duration=.015, interval=.00375,
             jobs=jobs, source_sha256={}, input_sha256={})
    (source/'protocol.json').write_text(json.dumps(p))
    # Mock only the acceptance evidence, not the backend, chemistry or worker.
    monkeypatch.setattr(runner, 'require_validation', lambda _: (host.config, []))
    gate = tmp_path/'gate'; gate.mkdir(); (gate/'protocol.json').write_text('{}')
    output = tmp_path/'gpu'
    runner.prepare(output, source, gate)
    original = digest(source/'protocol.json')
    runner.run(output)
    assert digest(source/'protocol.json') == original
    assert json.loads((output/'status.json').read_text())['completed'] == 2
    folder = output/name(jobs[1])
    result = json.loads((folder/'result.json').read_text())
    assert result['quality_pass'] and result['audit']['max_clipping'] == 0
    before = digest(folder/'latest_state.npz')
    history = json.loads((folder/'history.json').read_text())
    assert len(history) == 5 and history[0]['chemistry'][0][0] == pytest.approx(.9*host.activator[0])
    # Simulate loss of final result, retaining the accepted endpoint checkpoint.
    (folder/'result.json').unlink()
    runner.run(output)
    assert digest(folder/'latest_state.npz') == before
    assert json.loads((folder/'history.json').read_text()) == history
    restarted = AttributeSimulation.restore(folder/'latest_state.npz')
    assert restarted.time == .015 and restarted.ids.tolist() == [0]
    assert restarted.signal_rng.bit_generator.state == host.signal_rng.bit_generator.state
