import numpy as np
import pytest
import torch
from embryo.compensated_moving import validate_correction,MovingSimulation
from embryo.development_gpu import DevelopmentGpu,partition_chemistry,project_phase
from embryo.development_reference import partition_numpy,project_effective,ReferenceDevelopment
from embryo.model import Config,occupancy


@pytest.mark.parametrize('draws',[(0.,0.),(.04,-.03),(.2,-.2)])
def test_corrected_partition_conserves_amount_for_unequal_daughters(draws):
    u=np.array([.23,1.73]);e=.25*np.spacing(u);v=np.array([.17,.83])
    a,b=partition_numpy(u,e,1.00000004,v,draws)
    x,y=partition_chemistry(torch.from_numpy(u),torch.from_numpy(e),torch.tensor(1.00000004,dtype=torch.float64),torch.from_numpy(v),draws)
    np.testing.assert_array_equal(x.numpy(),a);np.testing.assert_array_equal(y.numpy(),b)
    np.testing.assert_allclose(a@v,u*1.00000004,rtol=2e-14,atol=0);validate_correction(a,b)
    assert np.any(b!=0)


def test_pressure_carries_transformed_sub_float32_state():
    s=ReferenceDevelopment(Config(grid=20,interface_width=.14,max_cells=2))
    p=s.phi[0];carry=np.full_like(p,1e-10,dtype=np.float64);target=float(occupancy(p.astype(float)+carry).sum()*s.dx**3)*1.001
    a,b,shift=project_effective(p,carry,target,s.dx)
    x,y,t=project_phase(torch.from_numpy(p),torch.from_numpy(carry),target,s.dx)
    np.testing.assert_allclose(x.numpy(),a,rtol=0,atol=2e-7);np.testing.assert_allclose(y.numpy(),b,rtol=0,atol=1e-7)
    assert np.any(b!=0) and shift>0
    np.testing.assert_allclose(occupancy(a.astype(float)+b).sum()*s.dx**3,target,rtol=1e-7)


def test_new_checkpoint_schema_cannot_be_read_as_mature(tmp_path):
    p=tmp_path/'bad.npz';np.savez_compressed(p,chemical_method=np.array('compensated resident progressive development v1'))
    with pytest.raises(ValueError,match='mandatory'):MovingSimulation.restore(p)
    with pytest.raises(ValueError,match='Mandatory'):DevelopmentGpu.restore(p)


def test_evolved_host_cannot_silently_reset_corrections():
    s=ReferenceDevelopment(Config(grid=20,interface_width=.14,max_cells=2));s.time=1.
    with pytest.raises(ValueError,match='mandatory'):DevelopmentGpu(s)


def test_corrected_reference_really_divides_and_preserves_lineage():
    s=ReferenceDevelopment(Config(grid=24,interface_width=.12,dt=.0075,max_cells=2,signal_partition_noise=.001))
    s.activator[:]=1.03;s.inhibitor[:]=.98;s.divide(0,[1,2,3])
    for _ in range(800):
        s.step()
        if len(s.ids)==2:break
    assert len(s.ids)==2 and not s.divisions and len(s.event_audits)==1
    a=s.event_audits[0];assert a['amount_relative_jump']<=2e-14 and a['target_relative_jump']<=2e-14
    assert a['volume_relative_jump']<=1e-6 and a['effective_occupancy_absolute_jump']<1e-14
    assert s.parents.tolist()==[0,0] and s.fate.tolist()==[0.,0.]
    assert np.any(s.phase_carry!=0) and np.any(s.chemical_correction!=0)
