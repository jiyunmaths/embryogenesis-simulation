import json
import numpy as np
import pytest

from embryo.model import Config, Simulation
from embryo.resolution import manufactured
from embryo.cutoff_dynamics import prepare, run


def test_explicit_baseline_polarity_cutoff_preserves_trajectory_and_checkpoint(tmp_path):
    inherited=manufactured(Config(grid=32,extent=2.24))
    explicit=manufactured(Config(grid=32,extent=2.24,polarity_contact_cutoff=.02))
    for _ in range(3):
        inherited.step(); explicit.step()
    for key in ('phi','activator','inhibitor','polarity','fate'):
        np.testing.assert_array_equal(getattr(inherited,key),getattr(explicit,key))
    explicit.checkpoint(tmp_path/'state.npz')
    restored=Simulation.restore(tmp_path/'state.npz')
    assert restored.config.polarity_contact_cutoff==.02
    explicit.step(); restored.step()
    np.testing.assert_array_equal(explicit.phi,restored.phi)


def test_signaling_cutoff_does_not_change_fixed_polarity_graph(monkeypatch):
    import embryo.model as model
    graphs=[]
    original=model.evolve
    def record(p,c,a,delta,*args):
        graphs.append(delta.copy())
        return original(p,c,a,delta,*args)
    monkeypatch.setattr(model,'evolve',record)
    for cutoff in (.01,.04):
        sim=manufactured(Config(grid=32,extent=2.24,graph_contact_cutoff=cutoff,polarity_contact_cutoff=.02))
        sim.step()
    np.testing.assert_array_equal(graphs[0],graphs[1])


@pytest.mark.parametrize('value',[-.5,1.,float('nan'),float('inf')])
def test_invalid_polarity_cutoff(value):
    with pytest.raises(ValueError,match='polarity_contact_cutoff'):
        Config(polarity_contact_cutoff=value).validate()


def test_matched_study_completes_and_preserves_source(tmp_path):
    checkpoint=tmp_path/'source.npz'
    manufactured(Config(grid=32,extent=2.24)).checkpoint(checkpoint)
    before=checkpoint.read_bytes()
    output=tmp_path/'study'
    prepare(checkpoint,output,until=.03,interval=.015)
    protocol=(output/'protocol.json').read_bytes()
    result=run(output)
    assert result['source_unchanged']
    assert checkpoint.read_bytes()==before
    assert (output/'protocol.json').read_bytes()==protocol
    assert len(result['pairs'])==2
    for cutoff in (.02,.01,.04):
        report=json.loads((output/f'cutoff-{cutoff:g}'/'analysis.json').read_text())
        assert len(report['history'])==3
        assert report['config']['polarity_contact_cutoff']==.02
        assert report['config']['graph_contact_cutoff']==cutoff
        assert (output/f'cutoff-{cutoff:g}'/'viewer.html').exists()
    assert json.loads((output/'status.json').read_text())['state']=='completed'
    with pytest.raises(ValueError):run(output)
    with pytest.raises(FileExistsError):prepare(checkpoint,output)
