from dataclasses import asdict
import numpy as np
from embryo.attribute_development import AttributeSimulation
from embryo.feedback_survival import initialize,similarity
from embryo.model import Config


def test_switch_changes_only_feedback_flag_and_keeps_off_exactly_continuable(tmp_path):
    source=AttributeSimulation(Config(grid=24,max_cells=1,interface_width=.12),'no_feedback')
    source.activator[:]=1.3;source.inhibitor[:]=.8;source.polarity[:]=[.1,.2,0]
    path=tmp_path/'source.npz';source.checkpoint(path)
    off=initialize(path,'keep_off');on=initialize(path,'switch_on')
    before=asdict(source.config);after=asdict(on.config)
    assert [k for k in before if before[k]!=after[k]]==['feedback']
    for sim in [off,on]:
        for key in ['phi','activator','inhibitor','polarity','fate','ids','parents','target','due']:
            np.testing.assert_array_equal(getattr(sim,key),getattr(source,key))
        assert sim.time==source.time and sim.step_number==source.step_number
        assert sim.signal_rng.bit_generator.state==source.signal_rng.bit_generator.state
    for _ in range(2):source.step();off.step()
    for key in ['phi','activator','inhibitor','polarity']:
        np.testing.assert_array_equal(getattr(off,key),getattr(source,key))


def test_similarity_distinguishes_rearrangement_and_uniformity():
    x=np.array([[.2,1.,2.],[.5,1.,2.]])
    v=np.ones(3)
    same=similarity(x,x,v)
    assert abs(same['initial_log_activator_correlation']-1)<1e-14
    assert same['initial_chemical_log_rms']==0
    assert similarity(x[:,::-1],x,v)['initial_log_activator_correlation']<0
    assert similarity(np.ones_like(x),x,v)['initial_log_activator_correlation'] is None
