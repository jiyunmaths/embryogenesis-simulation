from copy import deepcopy
import numpy as np
from embryo.model import Config,Simulation
from embryo.surface import cell_mesh,viewer_template


def test_sphere_mesh_is_closed_outward_and_at_correct_physical_location():
    sim=Simulation(Config(grid=32,max_cells=1))
    m=cell_mesh(sim.phi[0],sim.dx,sim.config.extent)
    v=np.asarray(m['vertices']);f=np.asarray(m['faces']);n=np.asarray(m['normals'])
    assert m['closed'] and m['boundary_edges']==0
    edges=np.sort(np.concatenate([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]]),axis=1)
    assert np.all(np.unique(edges,axis=0,return_counts=True)[1]==2)
    assert len(v)-len(edges)//2+len(f)==2
    np.testing.assert_allclose(v.mean(axis=0),0,atol=1e-6)
    np.testing.assert_allclose(np.linalg.norm(v,axis=1),.8,atol=.015)
    assert np.all(np.einsum('ij,ij->i',v,n)>0)
    volume=np.einsum('ij,ij->i',v[f[:,0]],np.cross(v[f[:,1]],v[f[:,2]])).sum()/6
    assert abs(volume/(4*np.pi*.8**3/3)-1)<.03


def test_nonconvex_neck_is_preserved():
    sim=Simulation(Config(grid=40,max_cells=1))
    radius1=np.sqrt((sim.xyz[0]-.4)**2+sim.xyz[1]**2+sim.xyz[2]**2)
    radius2=np.sqrt((sim.xyz[0]+.4)**2+sim.xyz[1]**2+sim.xyz[2]**2)
    field=.5*(1-np.tanh((np.minimum(radius1,radius2)-.55)/.1))
    m=cell_mesh(field,sim.dx,sim.config.extent);v=np.asarray(m['vertices'])
    assert m['closed']
    neck=v[np.abs(v[:,0])<sim.dx]
    assert len(neck)>0 and np.max(np.linalg.norm(neck[:,1:],axis=1))<.48
    assert np.max(np.linalg.norm(v[:,1:],axis=1))>.53


def test_boundary_crossing_is_not_artificially_closed():
    field=np.zeros((12,12,12));field[6:]=1.
    m=cell_mesh(field,.1,.6)
    assert len(m['faces'])>0 and not m['closed'] and m['boundary_edges']>0
    assert cell_mesh(np.zeros_like(field),.1,.6)['faces']==[]


def test_mesh_export_preserves_state_rng_and_future_steps():
    sim=Simulation(Config(grid=16,interface_width=.16,max_cells=1));before=deepcopy(sim)
    cells=sim.surfaces();assert cells[0]['mesh']['closed']
    assert cells[0]['points']==sim.surfaces(include_mesh=False)[0]['points']
    assert 'mesh' not in sim.surfaces(include_mesh=False)[0]
    sim.step();before.step()
    for key in ('phi','fate','activator','inhibitor','polarity'):
        np.testing.assert_array_equal(getattr(sim,key),getattr(before,key))
    for name in ('rng','fate_rng','signal_rng'):
        assert getattr(sim,name).bit_generator.state==getattr(before,name).bit_generator.state


def test_saved_viewer_inlines_renderer_without_cdn():
    template=viewer_template()
    assert '__CELL_SURFACE_RENDERER__' not in template
    assert 'root.CellSurface=' in template and 'CellSurface.draw' in template
    assert 'src="http' not in template
