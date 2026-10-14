import json
import numpy as np
import pytest
from embryo import moving_recipient_amount_response as study
from embryo.compensated_moving import retime
from test_phase_carry_exchange_response import source_fixture


def test_initial_physical_source_retiming_preserves_states_and_random_streams(tmp_path):
    host,source=source_fixture(tmp_path,start=510.)
    before={k:getattr(host,k).copy() for k in ('phi','polarity','target','ids','parents','due','activator','inhibitor')}
    time=host.time;old_step=host.step_number;old_save=host.config.save_every
    rng=host.rng.bit_generator.state.copy();retime(host,host.config.dt/2)
    for k,x in before.items():np.testing.assert_array_equal(getattr(host,k),x)
    assert host.time==time and host.step_number==2*old_step and host.config.save_every==2*old_save and host.rng.bit_generator.state==rng


def response_fixture(tmp_path):
    times=np.arange(401)*.15;jobs=[]
    for C,tau in [('S',1.),('R',4.)]:
        for dose in (None,'fractional','amount'):
            factor=1. if dose is None else .9 if dose=='fractional' or C=='S' else .95
            source=tmp_path/(C+'_'+str(dose)+'.npz');study.save_npz(source,initial=np.ones((2,2)))
            study.write_json(source.with_suffix('.json'),dict(record=dict(factor=factor)))
            for level,_ in study.LEVELS:
                key=C+'_'+str(dose)+'_'+level;job=dict(key=key,context=C,stage='response',dose=dose,level=level,source=str(source));jobs.append(job)
                base=np.exp(.01*times[:,None,None])*np.ones((401,2,2));path=base.copy()
                if dose:path[:,0,0]*=np.exp(np.log(factor)*np.exp(-times/tau))
                rows=[dict(elapsed=float(t),chemistry=u.tolist(),ids=[19,20],volumes=[1.,2.]) for t,u in zip(times,path)]
                study.write_json(tmp_path/key/'history.json',rows)
    study.write_json(tmp_path/'protocol.json',{})
    return dict(jobs=jobs,pilot_duration=6.),dict(seed=8,cell=19,contexts=dict(S='S',R='R'))


def test_signed_response_uses_own_moving_control_and_actual_dose(tmp_path):
    p,pair=response_fixture(tmp_path);report=study.response_contrast(tmp_path,p,pair,60.)
    assert report['passed'] and all(x['fine_above_screen'] for x in report['rows'])
    a,b=report['rows'];assert a['effects']['fine']==pytest.approx(b['effects']['fine'],abs=1e-14)
    assert all(x['direct_signed_waveform_error']==0 for x in report['rows'])


def test_history_rejects_an_altered_initial_intervention(tmp_path):
    source=tmp_path/'initial.npz';study.save_npz(source,initial=np.ones((2,2)),ids=np.array([19,20]))
    job=dict(source=str(source),start=510.)
    rows=[dict(elapsed=t,time=510+t,chemistry=np.ones((2,2)).tolist(),ids=[19,20],volumes=[1.,2.],delta=[[-1.,1.],[.5,-.5]]) for t in [0.,.15,.3]]
    study.check_history(rows,job,dict(interval=.15),.3)
    rows[0]['chemistry'][0][0]=.8
    with pytest.raises(ValueError,match='intervention'):study.check_history(rows,job,dict(interval=.15),.3)
