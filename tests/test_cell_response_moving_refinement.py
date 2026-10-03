import json
import numpy as np
import pytest
from embryo.cell_response_moving_refinement import response_error, assess, CRITERIA
from embryo.feedback_long import digest


def test_response_error_cancels_common_drift_but_detects_changed_response():
    base=np.ones((3,2,2));pulse=base.copy();pulse[0,0,0]=1.1
    drift=np.exp(np.arange(3)[:,None,None]*.1)
    assert response_error(pulse,base,pulse*drift,base*drift,1.1)<1e-14
    changed=pulse*drift;changed[1,0,0]*=1.01
    assert response_error(pulse,base,changed,base*drift,1.1)>.1
    with pytest.raises(ValueError):response_error(pulse,base,changed[:2],base,1.1)


def test_assessment_checks_recovery_classification_and_response(tmp_path):
    base=tmp_path/'base';fine=tmp_path/'fine';base.mkdir();fine.mkdir()
    jobs=[dict(family='pattern',target=None,factor=1.),dict(family='pattern',target=1,factor=1.1)]
    from embryo.cell_response_moving import name
    for root in (base,fine):
        p=dict(source_sha256={},input_sha256={},baseline=str(base),refinement_criteria=CRITERIA,jobs=jobs,duration=60.,interval=30.)
        (root/'protocol.json').write_text(json.dumps(p))
        np.savez(root/'initial_states.npz',ids=[1,2])
        for job in jobs:
            folder=root/name(job);folder.mkdir()
            (folder/'result.json').write_text(json.dumps(dict(quality_pass=True,protocol_sha256=digest(root/'protocol.json'))))
            h=[dict(elapsed=t,ids=[1,2],chemistry=[[1.,1.],[1.,1.]]) for t in (0.,30.,60.)]
            if job['target']:h[0]['chemistry'][0][0]=1.1
            (folder/'history.json').write_text(json.dumps(h))
        report=dict(trials=[dict(job=jobs[1],moving=dict(target_activator_log_auc_per_log_pulse=2.,target_recovery_time=3.,network_recovery_time=4.))])
        (root/'comparison.json').write_text(json.dumps(report))
    assert assess(fine)['passed']
    report['trials'][0]['moving']['target_recovery_time']=None
    (fine/'comparison.json').write_text(json.dumps(report))
    assert not assess(fine)['passed']
