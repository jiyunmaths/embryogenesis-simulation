import json
import pytest

from embryo.resolution import prepare as prepare_source, run as run_source
from embryo.refinement_confirmation import prepare, run, evaluate


def test_confirmation_preserves_source_and_records_every_check(tmp_path):
    source, output = tmp_path/'source', tmp_path/'confirmation'
    prepare_source(source, grids=(32,36,40), duration=.015)
    run_source(source)
    before = (source/'comparison.json').read_bytes()
    prepare(source, output, grid=44, probes=(36,40,44))
    frozen_protocol = (output/'protocol.json').read_bytes()
    result = run(output)
    assert (source/'comparison.json').read_bytes() == before
    assert (output/'protocol.json').read_bytes() == frozen_protocol
    assert len(result['pairs']) == 12
    assert len(result['calibration']) == 18
    assert result['checks']['archived_sources_unchanged']
    assert result['passed'] == all(result['checks'].values())
    assert json.loads((output/'status.json').read_text())['state'] == 'completed'
    with pytest.raises(ValueError):
        run(output)
    with pytest.raises(FileExistsError):
        prepare(source, output)


def test_quality_gates_reject_failed_variant():
    names = ['a','b','c']
    protocol = {'space_cases': names, 'criteria': {
        'field_relative_l2_max': .01, 'initial_analytic_error_max': .0025,
        'diagnostic_spread_max': .001, 'signal_relative_l2_max': .01,
        'axis_ratio_relative_difference_max': .01, 'boundary_occupancy_max': .01,
        'cell_volume_error_max': .05, 'minimum_radius_grid_cells': 4.}}
    rows = [{'cases': pair, 'probe_grid': probe, 'order': order, 'final_relative_l2': value}
            for probe in (88,112,128) for order in (3,5)
            for pair, value in [(names[:2],.004),(names[1:],.002)]]
    calibration = [{'relative_l2': .001}]
    reports = [{'max_boundary_occupancy': 0., 'max_cell_volume_error': .01,
                'max_clipped_fraction': 0., 'min_radius_grid_cells': 6.}]
    def check():
        return evaluate(protocol, rows, calibration, reports, {'activator': .001, 'inhibitor': .001}, .001, True)
    assert all(check().values())
    rows[-1]['final_relative_l2'] = .02
    checks = check()
    assert not checks['all_six_field_comparisons_below_1_percent']
    assert not checks['field_difference_decreases_for_every_diagnostic']
    assert not checks['diagnostic_spread_below_0_1_percentage_point']
    calibration[0]['relative_l2'] = .003
    assert not check()['initial_reconstruction_below_0_25_percent']
