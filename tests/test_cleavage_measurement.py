import json
import numpy as np
import pytest
from embryo.model import Config, Simulation, occupancy
from embryo.projection_audit import project, relative_error
from embryo.cleavage_measurement import analytic, reconstruct, prepare, run


@pytest.mark.parametrize('order',[3,5])
def test_native_projection_and_constants_are_preserved_without_mutation(order):
    field=np.full((2,16,16,16),.37,dtype=np.float32)
    saved=field.copy()
    np.testing.assert_array_equal(reconstruct(field,2.24,16,order),field)
    np.testing.assert_allclose(reconstruct(field,2.24,24,order),.37,atol=1e-7)
    np.testing.assert_array_equal(field,saved)


def test_known_failed_calibration_passes_when_native_phase_is_reconstructed():
    args=dict(extent=2.24,width=.085,radius=.4,center=[0,0,0],stretch=[1,1,1],angle=0.)
    field=analytic(56,**args); exact=occupancy(analytic(88,**args))
    sim=Simulation(Config(grid=56,extent=2.24));sim.phi=field
    assert relative_error(project(sim,88,3),exact)>.0025
    assert relative_error(occupancy(reconstruct(field,2.24,88,3)),exact)<.0025


def test_protocol_refuses_changed_source_and_preserves_original(tmp_path):
    source=tmp_path/'source';source.mkdir()
    cases={}
    for axis in ('axial','oblique'):
        for grid in (16,20,24):
            name=f'{axis}-{grid}'
            sim=Simulation(Config(grid=grid,extent=1.6,interface_width=.16))
            (source/name).mkdir();sim.checkpoint(source/name/'final_state.npz')
            cases[name]={'config':vars(sim.config)}
    (source/'protocol.json').write_text(json.dumps({'cases':cases,'grids':[16,20,24]}))
    (source/'comparison.json').write_text(json.dumps({'passed':False}))
    output=tmp_path/'validation';prepare(source,output)
    with pytest.raises(FileExistsError):prepare(source,output)
    original=(source/'comparison.json').read_bytes()
    (source/'comparison.json').write_text('{}')
    with pytest.raises(ValueError,match='source changed'):run(output)
    (source/'comparison.json').write_bytes(original)
    assert json.loads((source/'comparison.json').read_text())['passed'] is False
