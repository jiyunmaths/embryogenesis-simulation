"""Compare grouped and original layouts with unchanged chemical stepping."""
from pathlib import Path
import numpy as np
import pytest
from embryo import phase_carry_recipient_state_response_grouped as study
from test_phase_carry_recipient_state_response import fixture


def test_grouping_preserves_all_initial_values_and_uncoupled_transport(tmp_path):
    p,_,_,_=fixture(tmp_path)
    u,d,m,n=study.batch_layout(p)
    assert u.shape==(2,18,15) and d.shape==(18,18) and n==2
    assert [x['trial'] for x in m]==[t['key'] for t in p['trials']]
    assert len(m)==93 and u.shape[1]*u.shape[2]//n-len(m)==42
    np.testing.assert_array_equal(study.frozen.unpack(u,m,n),np.stack([study.frozen.load_context(p,t)[3] for t in p['trials']]))
    for item,trial in zip(m,p['trials']):
        _,delta,_,_,_=study.frozen.load_context(p,trial);start=item['context_index']*n
        np.testing.assert_array_equal(d[start:start+n,start:start+n],delta)
        mask=np.ones(d.shape[1],bool);mask[start:start+n]=False
        assert not np.any(d[start:start+n,mask])


def test_all_actual_toy_paths_match_original_layout_at_three_timesteps(tmp_path):
    p,_,_,_=fixture(tmp_path)
    old=tmp_path/'separate-layout';old.mkdir()
    study.write_json(old/'protocol.json',p)
    for level in p['levels']:
        study.frozen.production(old,p,level)
        study.production(tmp_path,p,level)
        with np.load(old/'production'/level['level']/'paths.npz') as a,np.load(tmp_path/'production'/level['level']/'paths.npz') as b:
            np.testing.assert_array_equal(a['trial_keys'],b['trial_keys'])
            np.testing.assert_allclose(a['trajectory'],b['trajectory'],rtol=1e-13,atol=1e-15)
        assert study.read(tmp_path/'production'/level['level']/'result.json')['padding_paths']==42


def test_grouped_checkpoint_resume_retains_exact_observation_clock(tmp_path,monkeypatch):
    p,_,_,_=fixture(tmp_path);level=p['levels'][0]
    original=study.frozen.checked_torch_step;calls=[0]
    def interrupted(*args,**kwargs):
        calls[0]+=1
        if calls[0]==90:raise RuntimeError('simulated interruption')
        return original(*args,**kwargs)
    with monkeypatch.context() as context:
        context.setattr(study.frozen,'checked_torch_step',interrupted)
        with pytest.raises(RuntimeError,match='simulated'):study.production(tmp_path,p,level)
    study.production(tmp_path,p,level)
    whole=tmp_path/'whole';whole.mkdir();study.write_json(whole/'protocol.json',p)
    study.production(whole,p,level)
    with np.load(tmp_path/'production'/'coarse'/'paths.npz') as a,np.load(whole/'production'/'coarse'/'paths.npz') as b:
        np.testing.assert_array_equal(a['times'],b['times'])
        np.testing.assert_array_equal(a['trajectory'],b['trajectory'])
