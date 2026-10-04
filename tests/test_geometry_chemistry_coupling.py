import numpy as np
import pytest
import torch
from scipy.integrate import solve_ivp

from embryo.geometry_chemistry_coupling import (Schedule, amount_rhs, amount_jacobian,
    production_step, reference_worker, worker, assess, CRITERIA)
from embryo.feedback_long import digest
from embryo.resolution import write_json


def schedule():
    times = np.array([0., .3, .6]); masses = np.array([[.8, 1.2], [1., 1.05], [1.2, .9]])
    g = np.array([[[0., w], [w, 0.]] for w in (.2, .3, .4)])
    return Schedule(times, masses, g)


def test_interpolation_preserves_conservation_and_positivity():
    s = schedule()
    for t in np.linspace(0., .6, 17):
        volume, delta = s.delta(t)
        np.testing.assert_allclose(volume@delta, 0, atol=1e-15)
        np.testing.assert_allclose(delta.sum(1), 0, atol=1e-15)
        assert np.all(volume > 0)
    with pytest.raises(ValueError, match='outside'): s.at(.61)
    broken = s.conductance.copy(); broken[0, 0, 1] += .1
    with pytest.raises(ValueError): Schedule(s.times, s.masses, broken)


def test_amount_jacobian_matches_finite_difference_with_changing_volumes():
    s = schedule(); t = .22; amount = np.array([.7, 1.4, .9, 1.2])
    fun = amount_rhs(s, 2., .02, .55); jac = amount_jacobian(s, 2., .02, .55)
    epsilon = 1e-6; directions = np.eye(4)*epsilon
    finite = np.column_stack([(fun(t, amount+d)-fun(t, amount-d))/(2*epsilon) for d in directions])
    np.testing.assert_allclose(jac(t, amount), finite, rtol=1e-8, atol=1e-8)


@pytest.mark.parametrize('method', ['beginning','midpoint'])
def test_volume_conversion_preserves_amount_exactly_without_reaction_or_flux(method):
    s = schedule(); a = torch.tensor([.7, 1.3], dtype=torch.float64); b = torch.tensor([1., .8], dtype=torch.float64)
    identity = lambda a,b,*args:(a,b)
    x,y,error = production_step(a,b,s,0.,.3,method,stepper=identity)
    before, after = s.at(0.)[0], s.at(.3)[0]
    np.testing.assert_allclose(after*np.stack((x.numpy(),y.numpy())), before*np.stack((a.numpy(),b.numpy())), rtol=1e-15)
    assert error < 2e-14


def test_beginning_is_first_order_and_midpoint_second_order_for_same_moving_input():
    s = schedule(); initial = np.array([[.7, 1.3], [1., .8]])
    ref = solve_ivp(amount_rhs(s,2.,.02,.55), (0.,.6), (initial*s.masses[0]).ravel(),
        method='DOP853', rtol=1e-12, atol=1e-14).y[:,-1].reshape(2,-1)/s.masses[-1]
    orders = {}
    for method in ('beginning','midpoint'):
        errors = []
        for dt in (.00375,.001875,.0009375):
            a,b = [torch.from_numpy(x.copy()) for x in initial]
            for step in range(round(.6/dt)):
                a,b,_ = production_step(a,b,s,step*dt,dt,method)
            errors.append(abs(np.log(np.stack((a.numpy(),b.numpy()))/ref)).max())
        orders[method] = np.log2(np.array(errors[:-1])/errors[1:])
    np.testing.assert_allclose(orders['beginning'], [1.,1.], atol=.02)
    np.testing.assert_allclose(orders['midpoint'], [2.,2.], atol=.02)


def test_end_to_end_references_resume_and_semantic_clock_validation(tmp_path):
    s = schedule(); file = tmp_path/'source.npz'
    # Sampling at .15 matches the experiment; the underlying test path is linear.
    t = np.arange(5)*.15; masses = np.array([s.at(x)[0] for x in t])
    g = np.array([[[0., .2+x/3], [.2+x/3, 0.]] for x in t]); initial=np.array([[.7,1.3],[1.,.8]])
    np.savez(file, times=t, masses=masses, conductance=g, initial=initial, ids=[3,8])
    references = [dict(key=f'test_{spacing}', context='test', spacing=spacing) for spacing in ('fine','coarse')]
    jobs = [dict(key=f'{spacing}_{method}_{dt}',context='test',spacing=spacing,method=method,dt=dt)
        for spacing in ('fine','coarse') for method in ('beginning','midpoint') for dt in (.00375,.001875,.0009375)]
    p=dict(contexts=[dict(key='test', source=str(file), duration=.6)], references=references, jobs=jobs,
        dts=[.00375,.001875,.0009375],interval=.15,beta=2.,da=.02,db=.55,criteria=CRITERIA,
        source_sha256={},input_sha256={str(file):digest(file)},interpretation='Manufactured numerical test')
    write_json(tmp_path/'protocol.json',p)
    for job in references: assert reference_worker((tmp_path,job))['passed']
    for job in jobs: assert worker((tmp_path,job))['accuracy_pass']
    report = assess(tmp_path)
    assert report['time_stepping_pass'] and report['spacing_screen_pass']
    job = jobs[0]; folder=tmp_path/job['key']; expected=digest(folder/'paths.npz')
    # Resume from the saved amount/concentration checkpoint, not the source.
    (folder/'result.json').unlink(); worker((tmp_path,job))
    assert digest(folder/'paths.npz') == expected
    with np.load(folder/'paths.npz') as z: payload={k:z[k].copy() for k in z.files}
    payload['times'][1]+=.01; np.savez_compressed(folder/'paths.npz',**payload)
    r=__import__('json').loads((folder/'result.json').read_text()); r['paths_sha256']=digest(folder/'paths.npz')
    write_json(folder/'result.json',r)
    with pytest.raises(ValueError,match='clock'): assess(tmp_path)
