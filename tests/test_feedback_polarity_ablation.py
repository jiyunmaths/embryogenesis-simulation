from dataclasses import asdict
import numpy as np
from embryo.attribute_development import AttributeSimulation
from embryo.model import Config
from embryo import feedback_long as original
from embryo.feedback_polarity_ablation import initialize, worker, assess, ARM
from embryo.resolution import write_json


def test_ablation_changes_only_polar_tension_and_worker_uses_same_state(tmp_path):
    sim=AttributeSimulation(Config(grid=24,max_cells=1,interface_width=.12),'no_feedback')
    sim.polarity[:]=[.1,.2,0]
    source=tmp_path/'source.npz';sim.checkpoint(source)
    signals=np.array([[1.3],[.9]])
    full=original.initialize(source,signals,'full');ablated=initialize(source,signals)
    a,b=asdict(full.config),asdict(ablated.config)
    assert [k for k in a if a[k]!=b[k]]==['polarity_tension']
    assert b['polarity_enabled'] and b['polarity_tension']==0
    for name in ['phi','polarity','activator','inhibitor','target','ids','parents','due']:
        np.testing.assert_array_equal(getattr(full,name),getattr(ablated,name))
    assert full.rng.bit_generator.state==ablated.rng.bit_generator.state
    assert full.signal_rng.bit_generator.state==ablated.signal_rng.bit_generator.state
    # Exercise real reused worker and restart with a short actual integration.
    root=tmp_path/'run';root.mkdir()
    np.savez_compressed(root/'initial_states.npz',formation=signals)
    write_json(root/'protocol.json',dict(checkpoint=str(source),start=0.,duration=.03,dt=.015,interval=.015,late_duration=.015))
    job=dict(state='formation',arm=ARM)
    result=worker((root,job));assert result['quality_pass'];assert ARM not in original.ARMS
    restored=AttributeSimulation.restore(root/f'formation_{ARM}'/'latest_state.npz')
    ablated.step();ablated.step()
    for name in ['phi','polarity','activator','inhibitor']:
        np.testing.assert_array_equal(getattr(ablated,name),getattr(restored,name))
    assert worker((root,job))==result


def test_assessment_does_not_assume_rescue():
    passed=dict(quality_pass=True,persistent=True);failed=dict(quality_pass=True,persistent=False)
    assert assess(passed,failed,passed)=='removing_polar_mechanics_restores_formation'
    assert assess(failed,failed,passed)=='suppression_persists_without_polar_mechanics'
    assert assess(passed,passed,passed)=='inconclusive_reference_outcomes'
    assert assess(dict(quality_pass=False,persistent=True),failed,passed)=='inconclusive_numerical_quality'
