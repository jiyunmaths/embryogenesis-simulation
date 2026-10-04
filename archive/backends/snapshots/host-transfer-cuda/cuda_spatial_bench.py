"""Snapshot and transfer microbenchmarks, without changing live simulations."""
import ctypes as C
import hashlib
import json
from pathlib import Path
import subprocess
import time
import numpy as np
from .native_mechanics import NativeSimulation
from .polarity import exposure_cue
from .cuda_benchmark import inputs,arguments,update,device
from .native_mechanics import kernel
from .feedback_long import digest
from .resolution import write_json


def library():
    source=Path(__file__).with_suffix('.cu')
    flags=['-O3','-std=c++17','-arch=sm_61','--fmad=false','--ftz=false','--prec-div=true','--prec-sqrt=true','-ccbin','/usr/bin/g++-12','-shared','-Xcompiler','-fPIC','-lcublas']
    key=hashlib.sha256(source.read_bytes()+source.with_name('cuda_mechanics.cu').read_bytes()+repr(flags).encode()+subprocess.check_output(['nvcc','--version'])).hexdigest()[:20]
    path=source.parent.parent/'outputs/cuda-kernel-cache'/f'spatial-{key}.so'
    path.parent.mkdir(exist_ok=True,parents=True)
    if not path.exists():subprocess.run(['nvcc',*flags,str(source),'-o',str(path)],check=True)
    lib=C.CDLL(str(path));f=np.ctypeslib.ndpointer(dtype=np.float32,flags='C_CONTIGUOUS');d=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
    lib.spatial_benchmark.argtypes=[f,d,C.c_int,C.c_int,C.c_double,C.c_double,C.c_int,d,f,d,f]
    lib.spatial_benchmark.restype=C.c_int;lib.transfer_benchmark.argtypes=[C.c_int,C.c_int,C.c_int,f]
    lib.transfer_benchmark.restype=C.c_int;lib.cuda_mechanics_error.restype=C.c_char_p
    return lib


def snapshot(lib,sim,repeats=20):
    cells=len(sim.phi);geom=np.empty((cells,7));contact=np.empty((cells,cells),np.float32);attraction=np.empty(sim.phi.shape);times=np.empty(4,np.float32)
    tension,adhesion=sim.material_coefficients()
    if lib.spatial_benchmark(sim.phi,adhesion,cells,sim.config.grid,sim.dx,sim.config.extent,repeats,geom,contact,attraction,times):raise RuntimeError(lib.cuda_mechanics_error().decode())
    reference_contact,_=sim.contacts();gpu_contact=contact.astype(float)*sim.dx**3;np.fill_diagonal(gpu_contact,0)
    refcue=exposure_cue(sim.phi,sim.dx);shell2=(sim.phi*(1-sim.phi))**2
    refatt=(adhesion@shell2.reshape(cells,-1)).reshape(sim.phi.shape)
    errors=dict(volume_relative=float(abs(geom[:,0]/sim.volumes()-1).max()),center_abs=float(abs(geom[:,1:4]-sim.centers()).max()),cue_abs=float(abs(geom[:,4:]-refcue).max()),contact_abs=float(abs(gpu_contact-reference_contact).max()),contact_relative_frobenius=float(np.linalg.norm(gpu_contact-reference_contact)/max(np.linalg.norm(reference_contact),1e-30)),adhesion_abs=float(abs(attraction-refatt).max()))
    return dict(cells=cells,grid=sim.config.grid,component_ms=dict(zip(('occupancy_shell','contact_sgemm','adhesion_dgemm','geometry_and_cue_reductions'),map(float,times))),errors=errors)


def run(output):
    root=Path(output);root.mkdir(parents=True,exist_ok=False);lib=library();print(device(),flush=True)
    paths=[Path('outputs/fertilization-cue/screen-dt-0.00375/latest_state.npz'),Path('outputs/fertilization-cue/baseline/latest_state.npz'),Path('outputs/cell-exchange-response-moving/unexchanged/source.npz'),Path('outputs/cell-exchange-response-moving/unexchanged/unexchanged_control/latest_state.npz')]
    rows=[]
    for path in paths:
        sim=NativeSimulation.restore(path);row=snapshot(lib,sim);row.update(checkpoint=str(path.resolve()),sha256=digest(path),time=sim.time)
        args=inputs(sim);out=np.empty_like(sim.phi);kernel()(*arguments(sim,args,out),1)
        row['mechanical_precision']=[]
        for precision in (64,32):
            gpu,seconds=update(sim,args,precision,30)
            row['mechanical_precision'].append(dict(precision=precision,resident_ms=1000*seconds,phi_abs=float(abs(gpu-out).max())))
        rows.append(row);write_json(root/'snapshots.json',rows);print(json.dumps(row),flush=True)
    transfers=[]
    for size in (24,96):
        for pinned in (0,1):
            times=np.empty(2,np.float32);start=time.perf_counter()
            if lib.transfer_benchmark(size,20,pinned,times):raise RuntimeError(lib.cuda_mechanics_error().decode())
            row=dict(mebibytes=size,pinned=bool(pinned),h2d_ms=float(times[0]),d2h_ms=float(times[1]),total_host_seconds=time.perf_counter()-start)
            transfers.append(row);print(row,flush=True)
    write_json(root/'transfers.json',transfers)
    # Performance-only replicated/resampled fields; these are not embryo predictions.
    from scipy.ndimage import zoom
    from .model import Config
    base=NativeSimulation.restore(paths[-1]);scaling=[]
    for cells,n in ((16,48),(16,96),(32,72)):
        sim=NativeSimulation(Config(grid=n,extent=base.config.extent,interface_width=base.config.interface_width), 'direct')
        fields=base.phi if n==72 else zoom(base.phi,(1,n/72,n/72,n/72),order=1)
        sim.phi=np.ascontiguousarray(np.tile(fields,(cells//16,1,1,1)));sim.ids=np.arange(cells);sim.activator=np.tile(base.activator,cells//16)
        row=snapshot(lib,sim,20);scaling.append(row);print('scaling',json.dumps(row),flush=True)
    write_json(root/'scaling.json',scaling)
    write_json(root/'metadata.json',dict(device=device(),scope='Frozen spatial components and synchronous transfer microbenchmarks under active CPU load. GPU reductions/cuBLAS change summation order. Synthetic scaling fields are not physical refinement studies. No complete GPU timestep, clipping, division, or long-term validation.',source_sha256={str(p.resolve()):digest(p) for p in [Path(__file__),Path(__file__).with_suffix('.cu'),Path(__file__).with_name('cuda_mechanics.cu')]}))


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True);a=p.parse_args();run(a.output)
