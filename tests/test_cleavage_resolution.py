import json
import numpy as np
import pytest

from embryo.model import Config
from embryo.cleavage_resolution import initial, prepare, run


def test_common_continuous_target_and_prescribed_direction():
    a=initial(Config(grid=32,extent=2.24,max_cells=2),[1,2,3])
    b=initial(Config(grid=40,extent=2.24,max_cells=2),[1,2,3])
    np.testing.assert_array_equal(a.target,b.target)
    np.testing.assert_allclose(a.divisions[0]['axis'],np.array([1,2,3])/np.sqrt(14))
    assert a.divisions[0]['start']==b.divisions[0]['start']==0
    assert len(a.phi)==len(b.phi)==1


def test_audit_isolates_abscission_from_reaction_step():
    sim=initial(Config(grid=32,extent=2.24,max_cells=2,signal_partition_noise=0.),[1,2,3])
    for _ in range(600):
        sim.step()
        if sim.abscission_audit:break
    assert sim.abscission_audit is not None
    assert max(sim.abscission_audit['relative_jumps'])<1e-6
    assert sim.abscission_audit['aggregate_occupancy_relative_jump']<1e-6
    assert len(sim.phi)==2


def test_short_unfinished_screen_records_failure_not_false_success(tmp_path):
    out=tmp_path/'study'
    p=prepare(out,grids=(32,36,40),until=.15)
    assert len(p['cases'])==8
    before=(out/'protocol.json').read_bytes()
    r=run(out)
    assert not r['passed']
    assert not r['checks']['all_divisions_complete']
    assert (out/'protocol.json').read_bytes()==before
    assert json.loads((out/'status.json').read_text())['state']=='completed'
    assert all((out/n/'final_state.npz').exists() for n in p['cases'])
    with pytest.raises(ValueError):run(out)
    with pytest.raises(FileExistsError):prepare(out)
