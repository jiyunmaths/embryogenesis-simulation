import numpy as np
from embryo.fast_mechanics import polar_fields,FastAttributeSimulation
from embryo.polarity import tension_field,flux_divergence
from embryo.attribute_development import AttributeSimulation
from embryo.model import Config


def test_compiled_polar_flux_matches_reference_and_conserves():
    rng=np.random.default_rng(3);n=12;extent=1.4;dx=2*extent/n
    axis=(np.arange(n)+.5)*dx-extent;xyz=np.array(np.meshgrid(axis,axis,axis,indexing='ij'))
    field=rng.random((n,n,n),dtype=np.float32);center=np.array([.1,-.2,.07]);p=np.array([.2,-.1,.3])
    g,f=polar_fields(field,dx,extent,center,p,1.2,.35,.085)
    reference=tension_field(xyz-center[:,None,None,None],p,1.2,.35,.085)
    np.testing.assert_allclose(g,reference,rtol=2e-15,atol=2e-15)
    np.testing.assert_allclose(f,flux_divergence(field,reference,dx),rtol=2e-13,atol=1e-12)
    assert abs(f.sum())<1e-10


def test_cached_values_refresh_and_do_not_escape_step():
    sim=FastAttributeSimulation(Config(grid=24,interface_width=.12,max_cells=1,dt=.0075),'direct')
    sim._fast_cache={};old=sim.volumes();sim.phi=sim.phi*.99;new=sim.volumes();assert np.all(new<old)
    sim._fast_cache=None;sim.step();assert sim._fast_cache is None
    old=sim.volumes();sim.phi*=.99;assert np.all(sim.volumes()<old)


def test_steps_and_checkpoint_match_reference(tmp_path):
    c=Config(grid=24,interface_width=.12,max_cells=2,dt=.0075)
    a=AttributeSimulation(c,'direct');a.polarity[:]=[.2,.1,-.1];path=tmp_path/'start.npz';a.checkpoint(path)
    b=FastAttributeSimulation.restore(path)
    for _ in range(20):a.step();b.step()
    np.testing.assert_allclose(a.phi,b.phi,atol=2e-7,rtol=2e-6)
    for key in ('activator','inhibitor','polarity'):np.testing.assert_allclose(getattr(a,key),getattr(b,key),atol=1e-7,rtol=1e-6)
    b.checkpoint(tmp_path/'fast.npz');d=FastAttributeSimulation.restore(tmp_path/'fast.npz');b.step();d.step()
    np.testing.assert_array_equal(b.phi,d.phi)


def test_division_fallback_matches_original(tmp_path):
    a=AttributeSimulation(Config(grid=24,interface_width=.12,max_cells=2,dt=.0075),'direct')
    a.polarity[:]=[.1,0,0];a.divide(0,direction=[1,0,0]);path=tmp_path/'source.npz';a.checkpoint(path);b=FastAttributeSimulation.restore(path)
    for _ in range(12):a.step();b.step()
    np.testing.assert_array_equal(a.phi,b.phi)
    np.testing.assert_array_equal(a.activator,b.activator)


def test_cache_invalidation_through_abscission(tmp_path):
    a=AttributeSimulation(Config(grid=24,interface_width=.12,max_cells=2,dt=.0075),'direct')
    a.divide(0,direction=[1,0,0]);path=tmp_path/'division.npz';a.checkpoint(path);b=FastAttributeSimulation.restore(path)
    for _ in range(240):
        a.step();b.step()
        if len(a.ids)==2 and not a.divisions:break
    assert len(a.ids)==len(b.ids)==2 and not b.divisions
    np.testing.assert_array_equal(a.ids,b.ids)
    np.testing.assert_allclose(a.phi,b.phi,atol=2e-7,rtol=2e-6)
    np.testing.assert_allclose(a.activator,b.activator,atol=1e-7,rtol=1e-6)
    assert a.lineage==b.lineage
