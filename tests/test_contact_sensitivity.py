import numpy as np
import pytest

from embryo.model import Config
from embryo.resolution import manufactured
from embryo.contact_sensitivity import scan, assess, run


def test_cutoff_scan_preserves_state_and_monotonically_removes_edges():
    sim = manufactured(Config(grid=32, extent=2.24))
    original = sim.phi.copy()
    rows = scan(sim, [0., .01, .02, .04, .08])
    assert [r['edges'] for r in rows] == sorted([r['edges'] for r in rows], reverse=True)
    baseline = rows[2]
    assert baseline['conductance_relative_change'] == 0
    assert baseline['spectrum_relative_l2_change'] == 0
    assert baseline['maximum_growth_rate_change'] == 0
    assert all(r['amount_conservation_residual'] < 1e-10 for r in rows)
    np.testing.assert_array_equal(sim.phi, original)
    assert sim.config.graph_contact_cutoff == .02


def test_report_preserves_checkpoint_and_rejects_overwrite(tmp_path):
    sim = manufactured(Config(grid=32, extent=2.24))
    checkpoint = tmp_path/'state.npz'
    sim.checkpoint(checkpoint)
    before = checkpoint.read_bytes()
    result = run([checkpoint], tmp_path/'report')
    assert result['sources_unchanged']
    assert checkpoint.read_bytes() == before
    assert len(result['cases'][0]['rows']) == 6
    with pytest.raises(FileExistsError):
        run([checkpoint], tmp_path/'report')


def test_local_gate_reports_growth_sensitivity_without_stress_case_failure():
    row = {'cutoff': .02, 'component_count_unchanged': True,
           'unstable_mode_count_unchanged': True, 'conductance_relative_change': 0.,
           'spectrum_relative_l2_change': 0., 'maximum_growth_rate_change': 0.,
           'constant_residual': 0., 'amount_conservation_residual': 0.}
    criteria = {'conductance_relative_max': .05, 'spectrum_relative_max': .05, 'growth_absolute_max': .01}
    stress = {**row, 'cutoff': .08, 'maximum_growth_rate_change': 1.}
    assert all(assess([row, stress], criteria).values())
    local = {**row, 'cutoff': .04, 'maximum_growth_rate_change': .02}
    assert not assess([row, local], criteria)['local_growth_change_below_0_01']
