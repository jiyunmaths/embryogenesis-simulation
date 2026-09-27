"""Axis diagnostics, persistence decisions, and state-preserving interventions."""

from copy import deepcopy
from dataclasses import asdict
import json

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from embryo.model import Config, Simulation
from embryo.shape import (shape_tensor, axial_angle, observe, summarize, apply_mode,
                          run, compare)


def points():
    return np.array([[3,0,0],[-3,0,0],[0,2,0],[0,-2,0],[0,0,1],[0,0,-1]]).T


def test_tensor_matches_independent_covariance_and_rotates_with_shape():
    cloud=points();weights=np.ones(6)
    reference=shape_tensor(weights,cloud)
    np.testing.assert_allclose(reference['covariance'],np.diag([3.,4/3,1/3]),atol=1e-14)
    assert reference['tensor_axis_ratio']==pytest.approx(3.)
    rotation=Rotation.from_rotvec([.37,-.2,.56]).as_matrix()
    shifted=rotation@cloud+np.array([4,-8,2])[:,None]
    actual=shape_tensor(weights,shifted)
    np.testing.assert_allclose(actual['covariance'],rotation@np.diag([3.,4/3,1/3])@rotation.T,atol=1e-14)
    assert actual['tensor_axis_ratio']==pytest.approx(3.)
    assert axial_angle(actual['principal_axis'],rotation@np.array([1,0,0]))<1e-5
    np.testing.assert_allclose(actual['centroid'],[4,-8,2],atol=1e-14)


def test_spherical_and_oblate_spectra_do_not_manufacture_a_unique_long_axis():
    sphere=np.c_[np.eye(3),-np.eye(3)]
    assert shape_tensor(np.ones(6),sphere)['principal_axis'] is None
    oblate=np.diag([2.,2.,1.])@sphere
    result=shape_tensor(np.ones(6),oblate)
    assert result['tensor_axis_ratio']==pytest.approx(2.)
    assert result['principal_axis'] is None


def test_axes_are_unoriented_and_zero_vectors_are_undefined():
    assert axial_angle([1,0,0],[-1,0,0])==0
    assert axial_angle([1,0,0],[0,1,0])==pytest.approx(90)
    assert axial_angle(None,[1,0,0]) is None
    assert axial_angle([0,0,0],[1,0,0]) is None


def history():
    return [dict(time=float(t),axis_ratio=1.3,principal_axis=[(-1.)**t,0,0],
                 dividing_cells=0,activator_std=.5,angle_to_first_cleavage_degrees=0.,
                 max_cell_volume_error=.01,boundary_occupancy=.001,clipped_fraction=0.,
                 min_radius_grid_cells=5.,max_cell_components_phi_0_5=1,contact_components=1)
            for t in range(11)]


def test_persistence_is_separate_from_signal_contrast_and_numerical_quality():
    rows=history();result=summarize(rows,5.,1.)
    assert result['finite_window_shape_persistence']
    assert result['numerical_screen_passed']
    assert result['maximum_late_axis_rotation_degrees']==0
    for row in rows:
        row['activator_std']=0.;row['min_radius_grid_cells']=3.9
    result=summarize(rows,5.,1.)
    assert result['finite_window_shape_persistence']
    assert not result['checks']['late_signaling_contrast']
    assert not result['numerical_screen_passed']


@pytest.mark.parametrize('change',[{'principal_axis':None},{'axis_ratio':1.01},{'dividing_cells':1},
                                    {'principal_axis':[0.,1.,0.]}])
def test_transient_or_unidentifiable_shape_cannot_pass_persistence(change):
    rows=history();rows[-1].update(change)
    assert not summarize(rows,5.,1.)['finite_window_shape_persistence']


def test_missing_late_window_or_recent_cleavage_is_not_persistent():
    assert not summarize(history()[8:],5.,1.)['finite_window_shape_persistence']
    assert not summarize(history(),5.,8.)['finite_window_shape_persistence']


def small():
    return Simulation(Config(grid=12,interface_width=.18,max_cells=1,competence_cells=1))


@pytest.mark.parametrize('mode',['full','no_feedback','no_polarity_tension'])
def test_interventions_preserve_physical_state_and_random_streams(mode):
    sim=small();sim.step();before=deepcopy(sim)
    apply_mode(sim,mode)
    for name in ('phi','target','fate','activator','inhibitor','polarity','ids','due'):
        np.testing.assert_array_equal(getattr(sim,name),getattr(before,name))
    for name in ('rng','fate_rng','signal_rng'):
        assert getattr(sim,name).bit_generator.state==getattr(before,name).bit_generator.state
    assert sim.time==before.time and sim.step_number==before.step_number
    assert sim.config.feedback==(mode!='no_feedback')
    assert sim.config.polarity_tension==(0 if mode=='no_polarity_tension' else before.config.polarity_tension)


def test_observation_does_not_change_fields_or_random_streams():
    sim=small();before=deepcopy(sim)
    row=observe(sim)
    assert row['principal_axis'] is None
    assert row['max_cell_components_phi_0_5']==1
    assert row['contact_components']==1
    for name in ('rng','fate_rng','signal_rng'):
        assert getattr(sim,name).bit_generator.state==getattr(before,name).bit_generator.state
    sim.step();before.step()
    np.testing.assert_array_equal(sim.phi,before.phi)


def test_short_checkpoint_run_exports_actual_frames_and_failed_persistence(tmp_path):
    sim=small();checkpoint=tmp_path/'start.npz';sim.checkpoint(checkpoint)
    directory=tmp_path/'full'
    result=run(directory,checkpoint=checkpoint,duration=.06,late_duration=.03,sample_interval=.015)
    assert not result['summary']['finite_window_shape_persistence']
    assert result['provenance']['origin']=='mature_checkpoint'
    assert len(result['history'])==5
    assert json.loads((directory/'analysis.json').read_text())==result
    saved=Simulation.restore(directory/'final_state.npz')
    for _ in range(4):sim.step()
    np.testing.assert_array_equal(saved.phi,sim.phi)
    for filename in ('viewer.html','trajectory.json','lineage.json','protocol.json'):
        assert (directory/filename).stat().st_size>0
    original=(directory/'analysis.json').read_bytes()
    with pytest.raises(FileExistsError):run(directory)
    assert (directory/'analysis.json').read_bytes()==original


def test_compare_retains_negative_causal_result_and_developmental_origin(tmp_path):
    directories=[]
    for i,(mode,origin,ratio) in enumerate([('full','mature_checkpoint',1.3),('no_feedback','mature_checkpoint',1.32),('no_feedback','zygote',1.4)]):
        rows=history()
        for r in rows:r['axis_ratio']=ratio
        provenance={'origin':origin,'sha256':'same' if origin=='mature_checkpoint' else None}
        parameters=asdict(Config(feedback=mode=='full'))
        report={'mode':mode,'provenance':provenance,'parameters':parameters,'history':rows,'summary':summarize(rows,5.,1.)}
        path=tmp_path/f'run{i}';path.mkdir();(path/'analysis.json').write_text(json.dumps(report));directories.append(path)
    result=compare(directories,tmp_path/'comparison')
    assert len(result['matched_branch_comparisons'])==1
    assert not result['matched_branch_comparisons'][0]['excess_above_0_05_throughout_late_window']
    assert result['matched_branch_comparisons'][0]['mean_axis_ratio_excess']==pytest.approx(-.02)
    assert len(result['developmental_comparisons'])==1
    assert not result['developmental_comparisons'][0]['matched_checkpoint']
    assert (tmp_path/'comparison'/'persistence.png').stat().st_size>0


def test_checkpoint_with_incomplete_population_is_rejected(tmp_path):
    sim=Simulation(Config(grid=12,interface_width=.18,max_cells=2));sim.checkpoint(tmp_path/'bad.npz')
    with pytest.raises(ValueError,match='finished all cleavages'):
        run(tmp_path/'out',checkpoint=tmp_path/'bad.npz')
    assert not (tmp_path/'out').exists()


def test_extension_preserves_ancestry_history_and_exact_state(tmp_path):
    from embryo.shape import extend
    sim=small();checkpoint=tmp_path/'start.npz';sim.checkpoint(checkpoint)
    source=tmp_path/'first';output=tmp_path/'extension'
    first=run(source,checkpoint=checkpoint,duration=.06,late_duration=.03,sample_interval=.015)
    original=(source/'analysis.json').read_bytes()
    result=extend(source,output,until=.12,late_duration=.03,sample_interval=.015)
    assert result['provenance']==first['provenance']
    assert len(result['history'])==9
    assert result['continuations'][0]['start_time']==.06
    assert result['summary']['late_window']==[.09,.12]
    assert (source/'analysis.json').read_bytes()==original
    payload=json.loads((output/'trajectory.json').read_text())
    assert len(payload['frames'])==9
    assert [f['metrics']['time'] for f in payload['frames']]==[r['time'] for r in result['history']]
    for _ in range(8):sim.step()
    np.testing.assert_array_equal(Simulation.restore(output/'final_state.npz').phi,sim.phi)
    with pytest.raises(ValueError,match='consistent completed source'):
        extend(source,tmp_path/'invalid',until=.03)
