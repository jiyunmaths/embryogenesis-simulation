"""Forced-response controls: external chemistry is not spontaneous patterning."""
from copy import deepcopy
from dataclasses import replace
import json

import numpy as np
import pytest

from embryo.model import Config, Simulation
from embryo.signal_patch import prescribed_patch, run, compare, CASES


def mature():
    sim=Simulation(Config(grid=16,interface_width=.16,max_cells=2,competence_cells=1,division_interval=.15))
    for _ in range(200):sim.step()
    assert len(sim.phi)==2 and not sim.divisions
    return sim


def test_patch_is_positive_mean_matched_and_rotation_translation_covariant():
    from scipy.spatial.transform import Rotation
    centers=np.array([[0,0,-1],[1,0,.2],[-1,0,.9]])
    weights=np.array([1.,2.,3.]);axis=np.array([.1,.3,1.])
    a,h=prescribed_patch(centers,weights,axis)
    assert np.all(a>0)
    assert np.average(a,weights=weights)==pytest.approx(1.)
    assert np.max(abs(a-1))==pytest.approx(.8)
    np.testing.assert_array_equal(h,np.ones(3))
    rotation=Rotation.from_rotvec([.2,.1,.4]).as_matrix()
    rotated,_=prescribed_patch(centers@rotation.T+5,weights,rotation@axis)
    np.testing.assert_allclose(a,rotated,atol=1e-14)


@pytest.mark.parametrize('kwargs',[{'axis':[0,0,0]},{'amplitude':1},{'width':0}])
def test_bad_patch_is_rejected(kwargs):
    with pytest.raises(ValueError):prescribed_patch([[0,0,-1],[0,0,1]],[1,1],**kwargs)


def test_prescribed_signals_replace_gm_but_drive_fate_and_release():
    sim=mature();reference=deepcopy(sim)
    a=np.array([.2,1.8]);h=np.ones(2)
    sim.step(prescribed_signals=(a,h));reference.step()
    np.testing.assert_array_equal(sim.activator,a)
    np.testing.assert_array_equal(sim.inhibitor,h)
    assert not np.array_equal(sim.fate,reference.fate)
    # Caller arrays are not owned by the simulation.
    sim.activator[0]=.3
    assert a[0]==.2
    sim.step()
    assert not np.array_equal(sim.inhibitor,h)


def test_invalid_clamp_does_not_modify_state():
    sim=mature();before=deepcopy(sim)
    for pair in [(np.array([-1,1]),np.ones(2)),(np.ones(3),np.ones(2)),(np.array([np.nan,1]),np.ones(2))]:
        with pytest.raises(ValueError,match='positive finite arrays'):sim.step(prescribed_signals=pair)
        np.testing.assert_array_equal(sim.phi,before.phi)
        np.testing.assert_array_equal(sim.activator,before.activator)
        assert sim.step_number==before.step_number


def test_feedback_off_mechanics_is_independent_of_forced_chemistry():
    sim=mature();sim.config=replace(sim.config,feedback=False)
    reference=deepcopy(sim)
    for _ in range(10):
        sim.step(prescribed_signals=(np.array([.2,1.8]),np.ones(2)))
        reference.step(prescribed_signals=(np.ones(2),np.ones(2)))
        np.testing.assert_array_equal(sim.phi,reference.phi)
    assert not np.array_equal(sim.fate,reference.fate)


def test_full_patch_protocol_comparison_exports_and_no_overwrite(tmp_path):
    sim=mature();checkpoint=tmp_path/'source.npz';sim.checkpoint(checkpoint)
    initial_phi=sim.phi.copy();directories=[]
    for case in CASES:
        path=tmp_path/case;directories.append(path)
        report=run(path,checkpoint,case,clamp_duration=.06,release_duration=.06,sample_interval=.015)
        assert len(report['history'])==9
        if case=='patch_no_feedback':
            assert all(row['mean_directional_tension_amplitude']==0. for row in report['history'])
        assert report['history'][4]['phase']=='clamp'
        assert report['history'][5]['phase']=='released'
        payload=json.loads((path/'trajectory.json').read_text())
        forced=np.asarray(report['protocol']['activator_by_id'])
        for frame in payload['frames'][:5]:np.testing.assert_array_equal(frame['graph']['activator'],forced)
        if case.startswith('uniform'):
            # Release includes concentration dilution under changing cell volumes.
            assert np.all(np.asarray(payload['frames'][-1]['graph']['activator']) > 0)
        else:
            assert not np.array_equal(payload['frames'][-1]['graph']['activator'],forced)
        assert 'not spontaneous symmetry breaking' in (path/'viewer.html').read_text()
        np.testing.assert_array_equal(Simulation.restore(checkpoint).phi,initial_phi)
        restored=Simulation.restore(path/'final_state.npz')
        assert restored.config.fate_tension==(.75 if case.endswith('strong') else .25)
        assert restored.config.feedback==(case!='patch_no_feedback')
    result=compare(directories,tmp_path/'comparison')
    assert len(result['comparisons'])==3
    assert (tmp_path/'comparison'/'response.png').stat().st_size>0
    assert (tmp_path/'comparison'/'final_shapes.png').stat().st_size>0
    with pytest.raises(FileExistsError):run(directories[0],checkpoint)
    mismatch=json.loads((directories[1]/'analysis.json').read_text())
    mismatch['parameters']['fate_tension']=.5
    (directories[1]/'analysis.json').write_text(json.dumps(mismatch))
    with pytest.raises(ValueError,match='match fate-tension strength'):
        compare(directories,tmp_path/'bad-strength')
    mismatch['parameters']['fate_tension']=.25
    (directories[1]/'analysis.json').write_text(json.dumps(mismatch))
    broken=json.loads((directories[1]/'analysis.json').read_text())
    broken['protocol']['axis']=[1,0,0]
    (directories[1]/'analysis.json').write_text(json.dumps(broken))
    with pytest.raises(ValueError,match='matched checkpoint'):
        compare(directories,tmp_path/'bad-comparison')
