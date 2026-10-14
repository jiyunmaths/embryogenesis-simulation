"""Checkpoint integration tests for the isolated chemistry branch."""
from pathlib import Path
import json
import numpy as np
import pytest
import torch
from embryo import compensated_recipient_response as study


def setup(root):
    root.mkdir();u=np.array([[.8,1.2],[1.1,.9]]);d=np.array([[-.2,.2],[.1,-.1]]);v=np.array([1.,2.]);ids=np.array([11,17])
    c=dict(key='context',source=str(root/'context.npz'))
    study.save_npz(c['source'],delta=d,volumes=v,initial=u,ids=ids)
    trials=[]
    for target in (None,11):
        key='control' if target is None else 'pulse';initial=u.copy()
        if target is not None:initial[0,0]*=.9
        t=dict(key=key,context=c['key'],target=target,source=str(root/(key+'.npz')));study.save_npz(t['source'],initial=initial);trials.append(t)
    p=dict(contexts=[c],trials=trials,beta=2.,da=.02,db=.55,duration=.6,interval=.15,checkpoint_interval=.3,
        runtime=dict(torch=str(torch.__version__),torch_cuda=torch.version.cuda),device='cpu',
        input_sha256=study.frozen.pin([c['source'],*[t['source'] for t in trials]]))
    study.write_json(root/'protocol.json',p);return p


def test_interrupted_production_resumes_with_identical_state_correction_and_history(tmp_path):
    a,b=tmp_path/'a',tmp_path/'b';pa,pb=setup(a),setup(b);l=dict(level='fine',dt=.001875)
    study.production(a,pa,l);r=study.production(b,pb,l,stop_after=160);assert r['interrupted']
    with np.load(b/'production/fine/checkpoint.npz') as z:assert np.any(z['correction']!=0)
    study.production(b,pb,l)
    with np.load(a/'production/fine/checkpoint.npz') as x,np.load(b/'production/fine/checkpoint.npz') as y:
        for k in ('state','correction','history'):np.testing.assert_array_equal(x[k],y[k])
    with np.load(a/'production/fine/paths.npz') as x,np.load(b/'production/fine/paths.npz') as y:np.testing.assert_array_equal(x['trajectory'],y['trajectory'])


def test_resume_rejects_missing_correction_even_with_recommitted_hash(tmp_path):
    root=tmp_path/'run';p=setup(root);l=dict(level='fine',dt=.001875);study.production(root,p,l,stop_after=160)
    file=root/'production/fine/checkpoint.npz'
    with np.load(file) as z:data={k:z[k].copy() for k in z.files if k!='correction'}
    study.save_npz(file,**data);commit=study.read(root/'production/fine/checkpoint-commit.json');commit['checkpoint_sha256']=study.digest(file);study.write_json(root/'production/fine/checkpoint-commit.json',commit)
    with pytest.raises(ValueError,match='mandatory checkpoint correction'):study.production(root,p,l)
