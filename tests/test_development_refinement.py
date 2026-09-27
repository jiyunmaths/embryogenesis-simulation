from copy import deepcopy
import json
import numpy as np
import pytest

from embryo.model import Config
from embryo.development_refinement import initial, prepare, run, record, pair_metrics, DevelopmentSimulation


def test_zygotes_share_target_and_random_state_without_field_resampling():
    a=initial(Config(grid=32,extent=2.24))
    b=initial(Config(grid=40,extent=2.24))
    np.testing.assert_array_equal(a.target,b.target)
    np.testing.assert_array_equal(a.due,b.due)
    assert a.rng.bit_generator.state==b.rng.bit_generator.state
    assert a.signal_rng.bit_generator.state==b.signal_rng.bit_generator.state
    assert a.phi.shape==(1,32,32,32) and b.phi.shape==(1,40,40,40)
    assert not a.divisions and not b.divisions


def test_instrumented_division_records_conservation_without_changing_dynamics():
    c=Config(grid=32,extent=2.24,max_cells=2)
    a=initial(c,DevelopmentSimulation);a.event_audits=[]
    b=initial(Config(**vars(c)))
    for s in (a,b):s.divide(0,[1,2,3])
    for _ in range(600):
        a.step();b.step()
        if a.event_audits:break
    assert len(a.event_audits)==1
    assert max(a.event_audits[0]['relative_jumps'])<1e-6
    for key in ('phi','activator','fate','polarity'):
        np.testing.assert_array_equal(getattr(a,key),getattr(b,key))


def test_distribution_comparison_does_not_assume_cell_id_correspondence():
    s=initial(Config(grid=32,extent=2.24));r=record(s)
    r['activator']=[.9,1.1];r['inhibitor']=[1.,1.2];r['volumes']=[1.,2.]
    r['metrics'].update(cells=2,fate_a=1,fate_b=1,uncommitted=0)
    a={'history':[r],'config':{'max_cells':2},'lineage':[{'division':1.}]}
    b=deepcopy(a)
    for key in ('activator','inhibitor','volumes'):b['history'][0][key].reverse()
    b['history'][0]['ids']=[500,99]
    p=pair_metrics(a,b)
    assert p['activator_max_wasserstein']==0
    assert p['inhibitor_max_wasserstein']==0


def test_short_pipeline_preserves_protocol_and_reports_unfinished_development(tmp_path):
    output=tmp_path/'screen'
    p=prepare(output,grids=(32,36,40),until=.15,interval=.15)
    assert len(p['cases'])==5
    assert {p['cases'][n]['dt'] for n in p['space_cases']}=={.0075}
    assert {p['cases'][n]['grid'] for n in p['time_cases']}=={36}
    frozen=(output/'protocol.json').read_bytes()
    r=run(output)
    assert not r['passed']
    assert not r['comparison_checks']['space_maturity_timing_agrees']
    assert (output/'protocol.json').read_bytes()==frozen
    assert len(json.loads((output/'status.json').read_text())['completed_cases'])==5
    assert all((output/n/'viewer.html').is_file() for n in p['cases'])
    with pytest.raises(ValueError):run(output)
    with pytest.raises(FileExistsError):prepare(output)
