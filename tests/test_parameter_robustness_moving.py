from dataclasses import asdict

import numpy as np
import pytest

from embryo import parameter_robustness_moving as study
from embryo.attribute_development import AttributeSimulation
from embryo.model import Config


def source(tmp_path):
    c = Config(grid=16, extent=.8, max_cells=2, dt=.00375, steps=40000,
               save_every=40, differentiation=False, signal_transport='conservative')
    host = AttributeSimulation(c, mode='direct')
    partition = .5*(1+np.tanh(host.xyz[0]/c.interface_width))
    field = host.phi[0]
    host.phi = np.array([field*partition, field*(1-partition)], dtype=np.float32)
    host.target = host.volumes()
    host.ids = np.array([3, 4]); host.parents = np.array([1, 1]); host.due = np.full(2, np.inf)
    host.fate = np.zeros(2); host.activator = np.array([.2, 1.8]); host.inhibitor = np.ones(2)
    host.polarity = np.array([[1., 0., 0.], [-1., 0., 0.]])
    host.step_number = 40000; host.time = 150.
    cp = tmp_path/'source.npz'; host.checkpoint(cp)
    chemicals = tmp_path/'chemicals.npz'
    np.savez(chemicals, ids=host.ids, pattern=np.array([host.activator, host.inhibitor]), uniform=np.ones((2, 2)))
    return host, cp, chemicals


def test_parameter_restart_preserves_geometry_lineage_clock_and_baseline_chemistry(tmp_path):
    host, cp, chemicals = source(tmp_path)
    point = dict(ratio=27.5, chi=.7)
    changed = study.initialize(cp, chemicals, 'pattern', point, asdict(host.config))
    assert changed.config.signal_da == .02 and changed.config.signal_dh == .55
    assert changed.config.polarity_tension == .7
    assert changed.time == host.time and changed.step_number == host.step_number
    for key in ('phi', 'polarity', 'ids', 'parents', 'target', 'activator', 'inhibitor'):
        np.testing.assert_array_equal(getattr(changed, key), getattr(host, key))
    assert changed.rng.bit_generator.state == host.rng.bit_generator.state
    assert changed.config.fate_tension == host.config.fate_tension
    assert changed.config.fate_adhesion == host.config.fate_adhesion
    assert AttributeSimulation.restore(cp).config.signal_dh == .4
    zero = study.initialize(cp, chemicals, 'uniform', dict(ratio=20., chi=0.), asdict(host.config))
    assert zero.config.polarity_enabled and zero.attribute_mode == 'direct'
    np.testing.assert_array_equal(zero.activator, np.ones(2))


def test_baseline_gate_is_enforced_before_new_parameter_changes(tmp_path):
    host, cp, chemicals = source(tmp_path)
    wrong = asdict(host.config); wrong['surface_tension'] *= 1.1
    with pytest.raises(ValueError, match='outside validated'):
        study.initialize(cp, chemicals, 'pattern', dict(ratio=27.5, chi=.7), wrong)
    np.savez(chemicals, ids=[4, 3], pattern=np.ones((2, 2)), uniform=np.ones((2, 2)))
    with pytest.raises(ValueError, match='cell order'):
        study.initialize(cp, chemicals, 'pattern', dict(ratio=20., chi=.35), asdict(host.config))


@pytest.mark.parametrize('ratio,chi', [(0., .35), (20., -1.), (20., 1.), (20., np.nan), (np.inf, .35)])
def test_invalid_or_negative_tension_regime_rejected(ratio, chi):
    with pytest.raises(ValueError): study.validate_point(dict(ratio=ratio, chi=chi))


def test_instantaneous_spectrum_uses_actual_graph_not_a_chi_scaling():
    masses = np.array([1., 2.]); delta = np.array([[-3., 3.], [1.5, -1.5]])
    row = dict(delta=delta.tolist(), volumes=masses.tolist(), time=150., chemistry=[[1, 1], [1, 1]])
    first = study.retain(row, dict(point=dict(ratio=20., chi=0.)))
    second = study.retain(row, dict(point=dict(ratio=20., chi=.7)))
    assert first == second
    third = study.retain(row, dict(point=dict(ratio=40., chi=.7)))
    assert third['uniform_growth_max'] > first['uniform_growth_max']


def test_failed_full_gate_or_incomplete_frozen_map_cannot_prepare(tmp_path, monkeypatch):
    from embryo.resolution import write_json
    frozen = tmp_path/'frozen'; frozen.mkdir()
    write_json(frozen/'protocol.json', dict(source_sha256={}, input_sha256={}))
    write_json(frozen/'status.json', dict(state='running'))
    write_json(frozen/'summary.json', dict(numerical_pass=False))
    with pytest.raises(ValueError, match='frozen map first'):
        study.prepare(tmp_path/'new', frozen, tmp_path/'validation')
    assert not (tmp_path/'new').exists()


def test_interrupted_moving_run_restores_chemistry_clock_and_observations_without_reinitializing(tmp_path, monkeypatch):
    from embryo import gpu_backend
    from embryo.resolution import write_json
    from embryo.neighbor_context import read
    host, cp, chemical = source(tmp_path)
    job = dict(key='pilot', seed=7, source=str(cp), chemical_file=str(chemical), family='pattern',
               point=dict(key='challenge', ratio=27.5, chi=.7))
    p = dict(accepted_config=asdict(host.config), start=150., duration=.3, dt=.00375, interval=.15,
             checkpoint_interval=.15, late_window=.15,
             criteria=dict(boundary_max=.01, dilution_error_max=2e-14, late_log_sd_min=.1))
    write_json(tmp_path/'protocol.json', p)
    fail = [True]
    class SmallCheckpointBackend:
        def __init__(self, host): self.host = host
        def __getattr__(self, key): return getattr(self.host, key)
        def to_cpu(self): return self.host
        def audit(self):
            return dict(max_volume_error=0., min_radius=5., max_clipping=0., dilution_amount_error=0.)
        def step(self):
            if fail[0] and self.host.step_number == 40045:
                fail[0] = False
                raise RuntimeError('simulated interruption after checkpoint')
            self.host.activator += 1e-5
            self.host.step_number += 1; self.host.time = self.host.step_number*self.host.config.dt
        def observe(self, elapsed):
            graph = self.host.signaling_graph()
            return dict(elapsed=elapsed, time=self.host.time, ids=self.host.ids.tolist(),
                chemistry=np.array([self.host.activator, self.host.inhibitor]).tolist(),
                volumes=graph.masses.tolist(), delta=graph.delta.tolist(), boundary_occupancy=0.,
                log_activator_sd=float(np.std(np.log(self.host.activator))))
    monkeypatch.setattr(gpu_backend, 'GpuSimulation', SmallCheckpointBackend)
    monkeypatch.setattr(study, 'prefix', lambda *args:dict(passed=True))
    monkeypatch.setattr(study, 'endpoint_assay', lambda *args:dict(numerical_pass=True))
    with pytest.raises(RuntimeError, match='simulated interruption'):
        study.continuation(tmp_path, job, p)
    saved = AttributeSimulation.restore(tmp_path/'pilot'/'latest_state.npz')
    assert saved.step_number == 40040 and saved.time == pytest.approx(150.15)
    r = study.continuation(tmp_path, job, p)
    final = AttributeSimulation.restore(tmp_path/'pilot'/'latest_state.npz')
    assert final.step_number == 40080 and final.time == pytest.approx(150.3)
    np.testing.assert_allclose(final.activator, host.activator+80e-5, rtol=1e-12, atol=1e-14)
    history = read(tmp_path/'pilot'/'history.json')
    assert len(history) == 3
    np.testing.assert_allclose([r['elapsed'] for r in history], [0., .15, .3], atol=1e-14)
    assert study.continuation(tmp_path, job, p) == r
