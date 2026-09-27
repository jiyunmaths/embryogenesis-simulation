"""The optimized kernel must retain complete trajectories, including cleavage."""
from copy import deepcopy
from dataclasses import asdict

import numpy as np
import pytest

from embryo.model import Config, Simulation
from reference_mechanics import ReferenceSimulation


@pytest.mark.parametrize('options', [
    {}, {'feedback':False}, {'polarity_enabled':False}, {'polarity_tension':0.},
    {'fate_tension':.75},
])
def test_optimized_trajectory_is_exact_through_cytokinesis(options):
    sim=Simulation(Config(grid=16,interface_width=.16,max_cells=2,
                          competence_cells=1,division_interval=.15,**options))
    sim.fate[:]=.4
    sim.activator[:]=1.3
    sim.polarity[:]=[.05,-.1,.15]
    reference=deepcopy(sim);reference.__class__=ReferenceSimulation
    for _ in range(200):
        sim.step();reference.step()
        for key in ('phi','fate','activator','inhibitor','polarity','target','ids'):
            np.testing.assert_array_equal(getattr(sim,key),getattr(reference,key))
        assert sim.clipped_fraction==reference.clipped_fraction
        assert sim.volume_projection_max==reference.volume_projection_max
    assert len(sim.phi)==2
    assert sim.lineage==reference.lineage
    assert sim.divisions==reference.divisions
    assert asdict(sim.config)==asdict(reference.config)
    for name in ('rng','fate_rng','signal_rng'):
        assert getattr(sim,name).bit_generator.state==getattr(reference,name).bit_generator.state
