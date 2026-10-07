import os
from pathlib import Path

import numpy as np
import pytest
import torch

from embryo import moving_initiation_controls as study
from embryo.gpu_initiation_controls import ARMS, conservative_delta, checkpoint_fields, InitiationControlSimulation
from embryo.gpu_backend import gm_step
from embryo.gpu_precision_control import PrecisionSimulation
from embryo.resolution import write_json
from test_gpu_backend import source
from test_polarity_robustness import protocol, trace


def design():
    return study.job_design({s:dict(coarse=f'{s}-c',fine=f'{s}-f',chemical=f'{s}-chem',transport=f'{s}-g') for s in (7,8)})


def test_design_balances_two_histories_four_arms_and_two_timesteps():
    jobs=design()
    assert len(jobs)==16 and len({j['key'] for j in jobs})==16
    assert all(j['point']['chi']==0 and j['point']['ratio']==27.5 and j['family']=='uniform' for j in jobs)
    for seed in (7,8):
        assert {(j['arm'],j['dt']) for j in jobs if j['seed']==seed}=={(a,dt) for a in ARMS for dt in (.00375,.001875)}


def test_fixed_conductances_conserve_amount_with_changing_unequal_volumes():
    g=torch.tensor([[0.,.4,.2],[.4,0.,.7],[.2,.7,0.]],dtype=torch.float64)
    old=torch.tensor([.4,1.3,.8],dtype=torch.float64)
    new=torch.tensor([.6,.9,1.7],dtype=torch.float64)
    initial=conservative_delta(g,old); moving=conservative_delta(g,new)
    np.testing.assert_allclose(new.numpy()@moving.numpy(),0,atol=1e-15)
    np.testing.assert_allclose(moving.numpy().sum(1),0,atol=1e-15)
    assert not torch.equal(initial,moving)
    assert np.linalg.norm(new.numpy()@initial.numpy())>.1 # freezing Delta would be wrong
    x=torch.tensor([.6,1.4,2.],dtype=torch.float64)
    assert abs(float(new@(moving@x)))<1e-14


def test_independent_chemical_reference_matches_substepped_torch_on_irregular_capacities():
    g=torch.tensor([[0.,4.,2.],[4.,0.,7.],[2.,7.,0.]],dtype=torch.float64)
    delta=conservative_delta(g,torch.tensor([.4,1.3,.8],dtype=torch.float64))
    state=np.array([[.7,1.5,.9],[.8,1.3,1.1]])
    expected=study.cpu_chemical_step(state,delta.numpy(),.3,2.,.02,.55)
    actual=gm_step(*[torch.from_numpy(x.copy()) for x in state],delta,.3,2.,.02,.55)
    np.testing.assert_allclose(np.array(actual),expected,rtol=1e-14,atol=1e-14)
    assert np.min(expected)>0


@pytest.mark.parametrize('values,qualified,expected',[
    ((False,False,False,False),True,'dilution_and_contact_ablation_insufficient'),
    ((False,True,False,True),True,'removing_dilution_restores_formation'),
    ((False,False,True,True),True,'fixing_contact_conductances_restores_formation'),
    ((False,True,True,True),True,'each_single_intervention_restores_formation'),
    ((False,False,False,True),True,'joint_removal_required_in_tested_context'),
    ((False,True,False,False),True,'single_arm_rescue_but_combined_nonrescue_investigate_nonmonotonicity'),
    ((True,True,True,True),True,'baseline_forms_not_a_failure_rescue'),
    ((False,True,True,True),False,'unresolved_numerical_or_physical_checks')])
def test_rescue_classification_preserves_negative_and_nonmonotonic_results(values,qualified,expected):
    assert study.rescue_interpretation(dict(zip(ARMS,values)),qualified)==expected


def test_failed_pilot_blocks_only_its_own_new_pair(tmp_path,monkeypatch):
    jobs=design();p=dict(jobs=jobs,histories=[7,8],arms=list(ARMS),reused={j['key']:{} for j in jobs if j['arm']=='baseline'},pilot_duration=60.,duration=240.)
    calls=[]
    def advance(root,job,p,horizon):
        calls.append((job['seed'],job['arm'],horizon));return [dict(elapsed=horizon)]
    monkeypatch.setattr(study,'advance',advance)
    monkeypatch.setattr(study,'pair_report',lambda root,p,seed,arm,horizon,histories:dict(passed=not(seed==7 and arm=='fixed-conductances')))
    class Executor:
        def submit(self,*args):return object()
    done,pending,failed=study.run_pairs(tmp_path,p,Executor())
    assert done==10 and len(pending)==10 and failed==[dict(seed=7,arm='fixed-conductances')]
    assert all(h==60 for seed,arm,h in calls if seed==7 and arm=='fixed-conductances')
    assert any(h==240 for seed,arm,h in calls if seed==8 and arm=='fixed-conductances')


def test_pair_records_are_immutable_and_check_geometry_as_well_as_chemical_transport(tmp_path):
    p=protocol();write_json(tmp_path/'protocol.json',p)
    histories={level:[dict(r,geometric_delta=r['delta']) for r in trace()] for level in ('coarse','fine')}
    saved=study.pair_report(tmp_path,p,7,'combined',.3,histories)
    assert saved['passed'] and 'history 7' in saved['interpretation']
    assert study.pair_report(tmp_path,p,7,'combined',.3,histories)==saved
    write_json(Path(next(iter(saved['history_sha256']))),[])
    with pytest.raises(ValueError,match='immutable'):study.pair_report(tmp_path,p,7,'combined',.3,histories)


def payload():
    shape=(2,4,4,4)
    return dict(precision_arm=np.array('phase_carry'),initiation_arm=np.array('combined'),
        phase_carry=np.zeros(shape,dtype=np.float64),rounding=np.zeros(2,dtype=np.int32),
        initial_conductance=np.array([[0.,1.],[1.,0.]]),initial_volume=np.ones(2),
        amount_source=np.zeros((2,2)),last_amount_source=np.zeros((2,2)),
        last_relative_volume_change=np.zeros(2),conversion_error=np.array(0.))


@pytest.mark.parametrize('key,bad',[
    ('initial_conductance',np.array([[0.,1.],[2.,0.]])),
    ('initial_volume',np.array([1.,0.])),
    ('phase_carry',np.zeros((2,4,4,4),dtype=np.float32)),
    ('amount_source',np.array([[0.,np.nan],[0.,0.]]))])
def test_invalid_restart_capacity_carry_and_accounting_are_rejected_before_upload(key,bad):
    p=payload();checkpoint_fields(p,(2,4,4,4));p[key]=bad
    with pytest.raises(ValueError):checkpoint_fields(p,(2,4,4,4))


gpu=pytest.mark.skipif(os.environ.get('EMBRYO_CUDA_TESTS')!='1',reason='GPU opt-in required')


@gpu
def test_new_baseline_matches_unchanged_carry_multistep():
    new=InitiationControlSimulation(source());old=PrecisionSimulation(source(),'phase_carry')
    for _ in range(16):
        new.step();old.step();study.exact_state(new,old)
    assert new.audit()['volume_conversion_amount_error']<2e-14


@gpu
@pytest.mark.parametrize('arm',ARMS)
def test_changed_chemistry_amount_source_and_exact_restart_on_nonuniform_state(tmp_path,arm):
    sim=InitiationControlSimulation(source(),arm)
    for _ in range(8):
        old=sim.geometry[:,0].cpu().numpy().copy()
        _,d,_,g=sim.matrices();effective,_=sim.chemical_transport(d,g)
        state=np.array([sim.activator.cpu().numpy(),sim.inhibitor.cpu().numpy()])
        expected=study.cpu_chemical_step(state,effective.cpu().numpy(),sim.config.dt,
                                       sim.config.signal_beta,sim.config.signal_da,sim.config.signal_dh)
        before=expected.copy();sim.step();new=sim.geometry[:,0].cpu().numpy()
        if sim.dilution:expected*=old/new
        np.testing.assert_allclose(np.array([sim.activator.cpu().numpy(),sim.inhibitor.cpu().numpy()]),expected,rtol=1e-13,atol=1e-14)
        expected_source=np.zeros_like(expected) if sim.dilution else before*(new-old)[None,:]
        np.testing.assert_allclose(sim.last_amount_source.cpu().numpy(),expected_source,rtol=1e-11,atol=1e-14)
        assert sim.audit()['volume_conversion_amount_error']<2e-14
    initial=sim.initial_conductance.clone();path=tmp_path/'state.npz';sim.checkpoint(path)
    restored=InitiationControlSimulation.restore(path)
    for _ in range(3):
        sim.step();restored.step();study.exact_state(sim,restored,True)
    assert torch.equal(initial,restored.initial_conductance)
    row=restored.observe(.1)
    if sim.fixed_contacts:
        np.testing.assert_array_equal(row['chemical_conductance'],initial.cpu().numpy())
        assert not np.array_equal(row['chemical_conductance'],row['geometric_conductance'])
