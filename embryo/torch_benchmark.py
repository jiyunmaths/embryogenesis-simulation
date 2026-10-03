"""Opt-in PyTorch/CUDA snapshot benchmarks; production backend is untouched.

Inputs stay on GPU during timing. CUDA graphs replay fixed inputs, not evolving
embryos. No torch.compile/Triton dependency (Pascal is unsupported there).
"""
import argparse
import gc
import json
from pathlib import Path
import time
import numpy as np
import torch
from .native_mechanics import NativeSimulation, kernel
from .cuda_benchmark import inputs, arguments, update
from .cuda_spatial_bench import library, snapshot
from .polarity import exposure_cue
from .feedback_long import digest
from .resolution import write_json


def arrays(phi):
    h=phi*phi*(3-2*phi)
    shell=phi*(1-phi)
    return h,shell,(shell*shell).double(),h.sum(dim=0)


def contacts(shell,dx):
    flat=shell.flatten(1)
    result=(flat@flat.T).double()*(dx**3)
    result.fill_diagonal_(0)
    return result


def attraction(shell2,adhesion):
    return (adhesion@shell2.flatten(1)).reshape(shell2.shape)


def geometry(phi,h,shell,occupied,xyz,dx):
    mass=h.sum(dim=(1,2,3),dtype=torch.float64)
    centers=(h.flatten(1).double()@xyz.flatten(1).T)/mass[:,None]
    gradients=torch.stack(torch.gradient(phi,spacing=(dx,dx,dx),dim=(1,2,3),edge_order=1),dim=1)
    length=(gradients*gradients).sum(dim=1).sqrt().clamp_min(1e-10)
    normal=-gradients/length[:,None]
    free=1-(2*(occupied[None]-h)).clamp(0,1)
    cue=(normal*(shell*free)[:,None]).sum(dim=(2,3,4),dtype=torch.float64)
    cue=cue/shell.sum(dim=(1,2,3),dtype=torch.float64).clamp_min(1e-12)[:,None]
    cue=torch.where(torch.linalg.vector_norm(cue,dim=1)[:,None]<1e-6,torch.zeros_like(cue),cue)
    return mass*(dx**3),centers,cue


def mechanics(phi,excluded,attract,centers,polarity,tensions,vf,xyz,dx,width,contrast,repulsion,dt,precision=64):
    dtype=torch.float64 if precision==64 else torch.float32
    # Keep explicit operation order and reference float32 products.
    relative=xyz.to(dtype)[None]-centers.to(dtype)[:,:,None,None,None]
    x,y,z=relative[:,0],relative[:,1],relative[:,2]
    distance=(x*x+y*y+z*z+width*width).sqrt()
    p=polarity.to(dtype)
    dot=p[:,0,None,None,None]*(x/distance)+p[:,1,None,None,None]*(y/distance)+p[:,2,None,None,None]*(z/distance)
    gamma=tensions.to(dtype)[:,None,None,None]*(1-contrast*dot)
    flux=torch.zeros_like(gamma)
    for axis in (1,2,3):
        left=[slice(None)]*4;right=left.copy();left[axis]=slice(None,-1);right[axis]=slice(1,None)
        left,right=tuple(left),tuple(right)
        face=.5*(gamma[left]+gamma[right])*(phi[right]-phi[left])/(dx*dx)
        flux[left]+=face;flux[right]-=face
    derivative=2*phi*(1-phi)*(1-2*phi)
    force=width*width*flux-gamma*derivative
    force+=vf.to(dtype)[:,None,None,None]*6*phi*(1-phi)
    force-=repulsion*phi*excluded
    force+=derivative*attract.to(dtype)
    return (phi+dt*force).float()


def measure(fn,repeats):
    for _ in range(3):result=fn()
    torch.cuda.synchronize()
    values=[];wall=[]
    for _ in range(3):
        start,end=torch.cuda.Event(enable_timing=True),torch.cuda.Event(enable_timing=True)
        before=time.perf_counter();start.record()
        for _ in range(repeats):result=fn()
        end.record();end.synchronize()
        values.append(start.elapsed_time(end)/repeats);wall.append((time.perf_counter()-before)*1000/repeats)
    return result,dict(cuda_event_ms=float(np.median(values)),synchronized_wall_ms=float(np.median(wall)),trials_ms=values)


def graph_measure(fn,repeats):
    # Capture excludes allocation/capture setup from replay timing.
    stream=torch.cuda.Stream();stream.wait_stream(torch.cuda.current_stream())
    with torch.cuda.stream(stream):
        for _ in range(3):fn()
    torch.cuda.current_stream().wait_stream(stream);torch.cuda.synchronize()
    graph=torch.cuda.CUDAGraph()
    with torch.cuda.graph(graph):result=fn()
    _,timing=measure(graph.replay,repeats)
    cpu=tuple(x.cpu().numpy().copy() for x in result)
    del graph,result;gc.collect()
    return cpu,timing


def compare(reference,actual):
    return float(np.max(np.abs(np.asarray(reference)-np.asarray(actual))))


@torch.inference_mode()
def benchmark(sim,repeats=20):
    tensor=lambda a:torch.from_numpy(np.ascontiguousarray(a)).to('cuda:0')
    phi=tensor(sim.phi);xyz=tensor(sim.xyz);tension,adhesion=sim.material_coefficients();a=tensor(adhesion)
    prepared=inputs(sim);gpu=[tensor(x) for x in prepared]
    c=sim.config;common=(*gpu,xyz,sim.dx,c.interface_width,c.polarity_tension,c.repulsion,c.dt)
    # Match the mechanics function argument ordering to the original kernel.
    def mechanical(precision):return mechanics(*common,precision=precision)
    data,t0=measure(lambda:arrays(phi),repeats);h,shell,shell2,occupied=data
    contact,t1=measure(lambda:contacts(shell,sim.dx),repeats)
    attract,t2=measure(lambda:attraction(shell2,a),repeats)
    geom,t3=measure(lambda:geometry(phi,h,shell,occupied,xyz,sim.dx),repeats)
    cpu_contact,_=sim.contacts();cpu_cue=exposure_cue(sim.phi,sim.dx)
    cpu_attract=prepared[2]
    actual_contact=contact.cpu().numpy();vol,cent,cue=[x.cpu().numpy() for x in geom]
    errors=dict(volume_relative=float(abs(vol/sim.volumes()-1).max()),center_abs=compare(sim.centers(),cent),
        cue_abs=compare(cpu_cue,cue),contact_relative_frobenius=float(np.linalg.norm(actual_contact-cpu_contact)/max(np.linalg.norm(cpu_contact),1e-30)),adhesion_abs=compare(cpu_attract,attract.cpu().numpy()))
    cpu=np.empty_like(sim.phi);kernel()(*arguments(sim,prepared,cpu),1)
    forces=[]
    for precision in (64,32):
        result,timing=measure(lambda:mechanical(precision),repeats)
        other,seconds=update(sim,prepared,precision,50)
        forces.append(dict(precision=precision,pytorch=timing,custom_cuda_ms=seconds*1000,
            pytorch_phi_abs=compare(cpu,result.cpu().numpy()),custom_phi_abs=compare(cpu,other)))
    def pipeline():
        hh,ss,ss2,oo=arrays(phi)
        cc=contacts(ss,sim.dx);aa=attraction(ss2,a)
        vv,xx,pp=geometry(phi,hh,ss,oo,xyz,sim.dx)
        return vv,xx,pp,cc,aa
    eager,eager_time=measure(pipeline,repeats)
    graph_result,graph_time=graph_measure(pipeline,repeats)
    graph_error=max(compare(x.cpu().numpy(),y) for x,y in zip(eager,graph_result))
    custom=snapshot(library(),sim,50)
    return dict(cells=len(sim.ids),grid=c.grid,spatial_errors=errors,
        pytorch_component_timing=dict(zip(('occupancy_shell','contacts','adhesion','geometry_cues'),(t0,t1,t2,t3))),
        custom_cuda_spatial=custom,mechanics=forces,spatial_pipeline_eager=eager_time,
        spatial_pipeline_cuda_graph=graph_time,graph_vs_eager_max_abs=graph_error)


def run(output,repeats):
    torch.set_num_threads(1);torch.backends.cuda.matmul.allow_tf32=False
    if not torch.cuda.is_available():raise RuntimeError('CUDA PyTorch required')
    root=Path(output);root.mkdir(parents=True,exist_ok=False);rows=[]
    paths=[Path('outputs/fertilization-cue/screen-dt-0.00375/latest_state.npz'),Path('outputs/cell-exchange-response-moving/unexchanged/unexchanged_control/latest_state.npz')]
    for path in paths:
        sim=NativeSimulation.restore(path)
        row=benchmark(sim,repeats);row.update(checkpoint=str(path.resolve()),checkpoint_sha256=digest(path),model_time=sim.time)
        rows.append(row);write_json(root/'cases.json',rows);print(json.dumps(row),flush=True)
        gc.collect();torch.cuda.empty_cache()
    report=dict(device=torch.cuda.get_device_name(0),pytorch=torch.__version__,cuda_build=torch.version.cuda,repeats=repeats,
        source_sha256={str(p.resolve()):digest(p) for p in [Path(__file__),Path(__file__).with_name('native_mechanics.py'),Path(__file__).with_name('cuda_mechanics.cu'),Path(__file__).with_name('cuda_spatial_bench.cu')]},
        cases=rows,scope='Eager PyTorch and CUDA graph replay on fixed, GPU-resident inputs. CPU copies, warmup and graph capture excluded. Three repeated timing windows under ongoing CPU jobs. Custom CUDA comparison rerun in same session. No torch.compile, gradients, coupled evolution, division, or complete-step speedup claim. Spatial pipeline excludes mechanics, chemical integration, and scientific diagnostics.')
    write_json(root/'comparison.json',report)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True);parser.add_argument('--repeats',type=int,default=20);a=parser.parse_args()
    if a.repeats<1:parser.error('repeats must be positive')
    run(a.output,a.repeats)
