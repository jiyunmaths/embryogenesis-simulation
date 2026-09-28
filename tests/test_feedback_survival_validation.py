from dataclasses import asdict
import json
import numpy as np
import pytest
from embryo.attribute_development import AttributeSimulation
from embryo.model import Config
from embryo.feedback_survival_validation import retime, compare_refinement, development
from embryo.resolution import write_json


def test_retime_preserves_state_streams_and_next_physical_time(tmp_path):
    sim = AttributeSimulation(Config(grid=24, max_cells=1, interface_width=.12,
                                    dt=.0075, steps=80, save_every=8), 'no_feedback')
    sim.step(); sim.step()
    source = tmp_path/'source.npz'; sim.checkpoint(source)
    fine = retime(source, tmp_path/'fine.npz', .00375)
    for key in ('phi','activator','inhibitor','polarity','ids','parents','due','target','fate'):
        np.testing.assert_array_equal(getattr(sim,key), getattr(fine,key))
    for key in ('rng','fate_rng','signal_rng'):
        assert getattr(sim,key).bit_generator.state == getattr(fine,key).bit_generator.state
    assert fine.time == sim.time and fine.step_number == 2*sim.step_number
    before, after = asdict(sim.config), asdict(fine.config)
    assert {k for k in before if before[k] != after[k]} == {'dt','steps','save_every'}
    sim.step(); fine.step(); fine.step()
    assert fine.time == sim.time
    with pytest.raises(ValueError): retime(source, tmp_path/'bad.npz', .004)


def test_comparison_checks_entire_trajectory_and_cell_ids(tmp_path):
    baseline = tmp_path/'baseline'; refined = tmp_path/'refined'
    p = dict(baseline=str(baseline), criteria=dict(max_chemical_log_rms=.02,
                                                 max_relative_axis_ratio_error=.01))
    row = dict(metrics=dict(time=90., axis_ratio=1.2), cell_order=[1],
               cells=[dict(attributes=dict(log_activator=0.,log_inhibitor=0.),context=dict(volume=1.))])
    result = dict(contrast_survives=True, initial_ordering_retained=True, quality_pass=True)
    for root in (baseline, refined):
        for arm in ('switch_on','keep_off'):
            (root/arm).mkdir(parents=True)
            write_json(root/arm/'history.json', [row]); write_json(root/arm/'result.json', result)
    assert all(r['passed'] for r in compare_refinement(tmp_path,p).values())
    row['cells'][0]['attributes']['log_activator'] = .1
    write_json(refined/'switch_on'/'history.json',[row])
    assert not compare_refinement(tmp_path,p)['switch_on']['passed']
    row['cell_order'] = [2]; write_json(refined/'switch_on'/'history.json',[row])
    with pytest.raises(ValueError, match='identities'): compare_refinement(tmp_path,p)


def test_development_checkpoint_resumes_without_reseeding(tmp_path):
    # Short one-cell fixture exercises fresh development, checkpoint and resume;
    # it is deliberately ineligible for the sixteen-cell scientific assay.
    cfg = Config(grid=24, extent=1.5, interface_width=.12, dt=.0075,
                 max_cells=1, steps=2, save_every=1)
    p = dict(source_config=asdict(cfg),start=.015,interval=.0075)
    write_json(tmp_path/'protocol.json',p)
    first = development(tmp_path,8,p)
    checkpoint = first['checkpoint']
    a = AttributeSimulation.restore(checkpoint)
    assert a.config.seed == 8 and a.time == .015 and not first['eligible_for_survival']
    (tmp_path/'seed-8'/'development'/'result.json').unlink()
    second = development(tmp_path,8,p)
    b = AttributeSimulation.restore(second['checkpoint'])
    np.testing.assert_array_equal(a.phi,b.phi)
    assert len(json.loads((tmp_path/'seed-8'/'development'/'history.json').read_text())) == 3
    assert second == first
