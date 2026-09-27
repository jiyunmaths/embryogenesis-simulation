import hashlib
import json
import numpy as np
import pytest

from embryo import Config, Simulation
from embryo.projection_audit import project, run as audit
from embryo.resolution import prepare, run as evolve


@pytest.mark.parametrize('order',[1,3,5])
def test_projection_preserves_constants_and_does_not_modify_fields(order):
    sim=Simulation(Config(grid=16,interface_width=.16))
    sim.phi[:]=.37
    original=sim.phi.copy()
    value=.37**2*(3-2*.37)
    projected=project(sim,24,order)
    np.testing.assert_allclose(projected,value,rtol=1e-6)
    np.testing.assert_array_equal(sim.phi,original)


def test_audit_calibrates_initial_fields_and_preserves_original_failure(tmp_path):
    source, output=tmp_path/'source',tmp_path/'audit'
    prepare(source,grids=(32,36,40),duration=.015)
    evolve(source)
    comparison=source/'comparison.json'
    before=hashlib.sha256(comparison.read_bytes()).hexdigest()
    original=json.loads(comparison.read_text())
    result=audit(source,output,probes=(32,40),orders=(1,3))
    assert len(result['pairs'])==8
    assert len(result['analytic_calibration'])==12
    assert result['summary']['original_screen_passed']==original['passed']
    assert hashlib.sha256(comparison.read_bytes()).hexdigest()==before
    native=[r for r in result['analytic_calibration'] if r['case']=='space-40' and r['probe_grid']==40]
    assert all(r['initial_relative_error_to_analytic']==0 for r in native)
    assert json.loads((output/'status.json').read_text())['original_comparison_unchanged']
    assert (output/'projection.png').is_file()
    with pytest.raises(FileExistsError):
        audit(source,output)
