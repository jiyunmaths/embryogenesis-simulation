"""Independent exact-solution and physical-geometry checks of material transport."""

import json

import numpy as np
import pytest
from scipy.integrate import quad, solve_ivp
from scipy.sparse.linalg import expm_multiply

from embryo.moving import MOTIONS, MovingTransport, motion, trajectory, time_refinement, run_benchmark
from embryo.transport import masked_cartesian_transport


@pytest.mark.parametrize('rates', [(.1,.1,.1), (.15,.05,.1), (.12,-.12,0.), (0.,0.,0.)])
def test_kinematics_match_volume_scaling_and_independent_clock_integrals(rates):
    scales, jacobian, clocks = motion(rates, 2.)
    np.testing.assert_allclose(scales, np.exp(2*np.array(rates)))
    assert jacobian == pytest.approx(np.exp(2*sum(rates)), abs=1e-15)
    for g, clock in zip(rates, clocks):
        assert clock == pytest.approx(quad(lambda t: np.exp(-2*g*t), 0., 2.)[0], rel=1e-13)


@pytest.mark.parametrize('rates', MOTIONS.values())
def test_amount_rhs_matches_actual_moving_faces_and_volumes(rates):
    model = MovingTransport(4, rates, length=1.7)
    time = 1.3
    scales = np.exp(np.array(rates)*time)
    current = masked_cartesian_transport(edges=tuple(e*s for e,s in zip(model.edges,scales)), mask=model.mask)
    concentration = np.random.default_rng(7).uniform(.5,1.5,len(current.volumes))
    amounts = current.volumes*concentration
    # Independently sum physical face fluxes G_ij (c_j-c_i), without using
    # the moving solver's directional matrices or amount-space generator.
    degree = np.asarray(current.conductance.sum(axis=1)).ravel()
    expected = .02*(current.conductance @ concentration - degree*concentration)
    np.testing.assert_allclose(model.rhs(time,amounts), expected, atol=2e-16, rtol=1e-12)
    assert abs(model.rhs(time,amounts).sum()) < 1e-15
    np.testing.assert_allclose(model.volumes(time), current.volumes, rtol=1e-14)


def test_uniform_dilution_includes_material_advection_and_geometric_divergence():
    model=MovingTransport(4, (.15,.05,.1))
    initial=model.reference.volumes.copy()
    report,fields=trajectory(model,duration=.2,dt=.002,snapshots=3)
    assert report['maximum_uniform_relative_error'] < 1e-13
    assert report['maximum_relative_amount_drift'] < 1e-13
    # Uniform quantities have no diffusive flux. Their concentration derivative
    # is -div(v)c even though the amount derivative is zero.
    t=.4;epsilon=1e-5
    c=initial/model.volumes(t)
    derivative=(initial/model.volumes(t+epsilon)-initial/model.volumes(t-epsilon))/(2*epsilon)
    np.testing.assert_allclose(derivative,-.3*c,rtol=1e-9)
    np.testing.assert_allclose(model.rhs(t,initial),0,atol=1e-16)


def test_no_diffusion_preserves_each_amount_during_anisotropic_deformation():
    model=MovingTransport(4,(.2,-.07,.03),diffusivity=0.)
    q=np.random.default_rng(12).uniform(.1,1.,len(model.reference.volumes))
    original=q.copy()
    for i in range(10):
        q=model.step(q,.1*i,.1)
    np.testing.assert_array_equal(q,original)
    c=q/model.volumes(1.)
    np.testing.assert_allclose(c/(original/model.reference.volumes),np.exp(-.16),rtol=1e-14)


def test_stationary_limit_matches_existing_semidiscrete_diffusion():
    model=MovingTransport(4,(0.,0.,0.))
    report,fields=trajectory(model,duration=.2,dt=.001,snapshots=2)
    exact=expm_multiply(.2*.02*model.reference.delta,model.exact(0.))
    np.testing.assert_allclose(fields['concentration'][-1],exact,atol=2e-8)
    assert report['maximum_uniform_relative_error'] < 1e-13


def test_exact_moving_cell_averages_against_gauss_quadrature_and_integrated_clock():
    rates=(.15,.05,-.08);model=MovingTransport(4,rates,length=1.3)
    t=.7;scales=np.exp(np.array(rates)*t);jacobian=np.prod(scales)
    clocks=np.array([quad(lambda s: np.exp(-2*g*s),0,t)[0] for g in rates])
    nodes,weights=np.polynomial.legendre.leggauss(8)
    values=[]
    for ijk in np.argwhere(model.mask):
        physical=[]
        for axis,k in enumerate(ijk):
            lo,hi=model.edges[axis][k:k+2]*scales[axis]
            physical.append((lo+hi)/2+nodes*(hi-lo)/2)
        x=np.meshgrid(*physical,indexing='ij')
        pointwise=np.ones(x[0].shape)
        for amplitude,mode in ((.1,(2,0,0)),(.07,(0,2,1)),(.04,(2,2,2))):
            factor=np.exp(-.02*np.pi**2*np.dot(np.square(mode),clocks)/1.3**2)
            pattern=np.ones_like(pointwise)
            for axis in range(3):
                pattern*=np.cos(mode[axis]*np.pi*x[axis]/(1.3*scales[axis]))
            pointwise+=amplitude*factor*pattern
        values.append(np.einsum('i,j,k,ijk',weights/2,weights/2,weights/2,pointwise)/jacobian)
    np.testing.assert_allclose(model.exact(t),values,atol=2e-14)


def test_heun_matches_independent_ode_with_rebuilt_geometry_at_each_rhs():
    model=MovingTransport(4,(.15,.05,.1))
    q0=model.reference.volumes*model.exact(0.)
    def independent_rhs(t,q):
        scales=np.exp(np.array([.15,.05,.1])*t)
        transport=masked_cartesian_transport(edges=tuple(e*s for e,s in zip(model.edges,scales)),mask=model.mask)
        c=q/transport.volumes
        return .02*(transport.conductance @ c - np.asarray(transport.conductance.sum(axis=1)).ravel()*c)
    exact=solve_ivp(independent_rhs,(0.,.4),q0,method='DOP853',rtol=1e-12,atol=1e-14)
    assert exact.success
    errors=[]
    for dt in (.02,.01,.005):
        q=q0.copy()
        for k in range(round(.4/dt)):
            q=model.step(q,k*dt,dt)
        errors.append(np.linalg.norm(q-exact.y[:,-1]))
    assert all(3.8<a/b<4.4 for a,b in zip(errors,errors[1:]))


def test_unsafe_step_and_invalid_amounts_rejected_without_mutation():
    model=MovingTransport(8,(.1,.1,.1))
    q=model.reference.volumes.copy();original=q.copy()
    with pytest.raises(ValueError,match='positivity bound'):
        model.step(q,0.,10.)
    np.testing.assert_array_equal(q,original)
    for bad in (-q,np.full_like(q,np.nan),q[:-1],q.astype(complex)):
        with pytest.raises(ValueError):
            model.step(bad,0.,.001)


def test_complete_report_retains_failed_short_expansion_control_and_snapshots(tmp_path):
    output=tmp_path/'moving'
    report=run_benchmark(output,resolutions=(4,8,16),duration=.2,dt=.002)
    assert json.loads((output/'analysis.json').read_text())==report
    assert report['checks']['amount_conserved_all_steps']
    assert report['checks']['uniform_dilution_matches_inverse_volume']
    assert not report['checks']['omitting_dilution_detectably_violates_balance']
    for name in ('fields.npz','convergence.png','balance.png','moving_fields.png'):
        assert (output/name).stat().st_size>0
    with np.load(output/'fields.npz',allow_pickle=False) as archive:
        assert archive['anisotropic_n16_concentration'].shape==(21,3072)
        assert all(np.isfinite(archive[key]).all() for key in archive.files)
    original=(output/'analysis.json').read_bytes()
    with pytest.raises(FileExistsError):
        run_benchmark(output)
    assert (output/'analysis.json').read_bytes()==original


@pytest.mark.parametrize('kwargs', [
    {'resolutions':(4,8)}, {'resolutions':(4,8,12)}, {'resolutions':(4,8,64)},
    {'duration':0}, {'duration':float('nan')}, {'dt':True}, {'dt':.03},
    {'dt':1.}, {'diffusivity':-1}, {'grading':1.}, {'length':0.},
])
def test_invalid_benchmarks_fail_before_creating_output(tmp_path,kwargs):
    output=tmp_path/'invalid'
    with pytest.raises(ValueError):
        run_benchmark(output,**kwargs)
    assert not output.exists()
