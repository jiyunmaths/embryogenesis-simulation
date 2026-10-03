"""Isolated CUDA precision pilot. Does not modify production simulation classes.

Set CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1 to use the 1080 Ti.
"""
import argparse
import ctypes
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time
import numpy as np
from .native_mechanics import NativeSimulation, kernel as cpu_kernel
from .feedback_long import digest
from .resolution import write_json

_LIBRARY = None


def library():
    global _LIBRARY
    if _LIBRARY is None:
        source = Path(__file__).with_name('cuda_mechanics.cu')
        flags = ['-O3', '-std=c++17', '-arch=sm_61', '--fmad=false', '--ftz=false',
                 '--prec-div=true', '--prec-sqrt=true', '-ccbin', '/usr/bin/g++-12',
                 '-shared', '-Xcompiler', '-fPIC']
        version = subprocess.check_output(['nvcc', '--version'])
        key = hashlib.sha256(source.read_bytes()+repr(flags).encode()+version).hexdigest()[:20]
        folder = source.parent.parent/'outputs/cuda-kernel-cache'
        folder.mkdir(parents=True, exist_ok=True)
        binary = folder/f'mechanics-{key}.so'
        if not binary.exists():
            temporary = binary.with_suffix(f'.{os.getpid()}.tmp.so')
            subprocess.run(['nvcc', *flags, str(source), '-o', str(temporary)], check=True)
            temporary.replace(binary)
        lib = ctypes.CDLL(str(binary))
        f = np.ctypeslib.ndpointer(dtype=np.float32, flags='C_CONTIGUOUS')
        d = np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS')
        lib.cuda_mechanics_update.argtypes = [f,f,d,d,d,d,d,f,ctypes.c_int,ctypes.c_int,
            *([ctypes.c_double]*6),ctypes.c_int,ctypes.c_int,ctypes.POINTER(ctypes.c_float)]
        lib.cuda_mechanics_update.restype = ctypes.c_int
        lib.cuda_mechanics_error.restype = ctypes.c_char_p
        lib.cuda_mechanics_device.argtypes = [ctypes.c_char_p,ctypes.c_int]
        lib.cuda_mechanics_device.restype = ctypes.c_int
        _LIBRARY = lib
    return _LIBRARY


def device():
    lib=library();name=ctypes.create_string_buffer(256)
    if lib.cuda_mechanics_device(name,256):
        raise RuntimeError(lib.cuda_mechanics_error().decode())
    return name.value.decode()


def inputs(sim):
    c=sim.config;phi=sim.phi
    shell2=(phi*(1-phi))**2;tension,adhesion=sim.material_coefficients()
    attract=np.ascontiguousarray((adhesion@shell2.reshape(len(phi),-1)).reshape(phi.shape))
    excluded=np.sum(phi**2,axis=0)[None]-phi**2
    vf=c.volume_stiffness*(sim.target-sim.volumes())/sim.target
    return (phi,excluded,attract,np.ascontiguousarray(sim.centers()),
            np.ascontiguousarray(sim.polarity),tension,np.ascontiguousarray(vf))


def arguments(sim, arrays, out):
    c=sim.config
    return (*arrays,out,len(sim.ids),c.grid,sim.dx,c.extent,c.interface_width,
            c.polarity_tension,c.repulsion,c.dt)


def update(sim, arrays, precision=64, repeats=1):
    if precision not in (32,64) or repeats<1:raise ValueError('Invalid precision/repeats')
    out=np.empty_like(sim.phi);elapsed=ctypes.c_float()
    lib=library()
    code=lib.cuda_mechanics_update(*arguments(sim,arrays,out),precision,repeats,ctypes.byref(elapsed))
    if code:raise RuntimeError(lib.cuda_mechanics_error().decode())
    return out,elapsed.value/1000


class CudaMechanicalPilot(NativeSimulation):
    precision=64

    def mechanical_step(self):
        c=self.config
        if self.divisions or not c.feedback or not c.polarity_enabled or np.any(np.linalg.norm(self.polarity,axis=1)<=1e-10):
            return super().mechanical_step()
        out,_=update(self,inputs(self),self.precision)
        self.volume_projection_max=0.
        self.clipped_fraction=float(np.mean((out<0)|(out>1)))
        self.phi=np.clip(out,0,1)


def discrepancies(a,b):
    return dict(phi_max_abs=float(abs(a.phi-b.phi).max()),
        chemical_max_log=float(abs(np.log(np.array([a.activator,a.inhibitor])/np.array([b.activator,b.inhibitor]))).max()),
        polarity_max_abs=float(abs(a.polarity-b.polarity).max()),
        volume_max_relative=float(abs(a.volumes()/b.volumes()-1).max()))


def run(checkpoint,output,steps=48,repeats=100):
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    name=device();print('GPU:',name,flush=True)
    reference=NativeSimulation.restore(checkpoint);reference.native_threads=1
    if reference.divisions:raise ValueError('Requires mature checkpoint')
    arrays=inputs(reference);cpu=np.empty_like(reference.phi);cpu_times=[]
    for _ in range(4):
        start=time.perf_counter();cpu_kernel()(*arguments(reference,arrays,cpu),1)
        cpu_times.append(time.perf_counter()-start)
    rows=[]
    for precision in (64,32):
        out,resident=update(reference,arrays,precision,repeats)
        calls=[]
        for _ in range(5):
            start=time.perf_counter();update(reference,arrays,precision,1);calls.append(time.perf_counter()-start)
        diff=out.astype(float)-cpu
        row=dict(precision=precision,resident_kernel_seconds=resident,
                 host_call_seconds=float(np.median(calls)),phi_max_abs=float(abs(diff).max()),
                 phi_rms=float(np.sqrt(np.mean(diff**2))),different_voxels=int(np.count_nonzero(diff)),
                 finite=bool(np.isfinite(out).all()))
        rows.append(row);print(json.dumps(row),flush=True)
    sims={p:CudaMechanicalPilot.restore(checkpoint) for p in (64,32)}
    for p,sim in sims.items():sim.precision=p;sim.native_threads=1
    maxima={p:dict.fromkeys(discrepancies(reference,sim),0.) for p,sim in sims.items()}
    timings={str(p):[] for p in sims};timings['cpu']=[]
    for step in range(steps):
        for label,sim in [('cpu',reference),*[(str(p),sim) for p,sim in sims.items()]]:
            start=time.perf_counter();sim.step();timings[label].append(time.perf_counter()-start)
        for p,sim in sims.items():
            for k,v in discrepancies(reference,sim).items():maxima[p][k]=max(maxima[p][k],v)
        if (step+1)%8==0:print(json.dumps(dict(step=step+1,errors=maxima)),flush=True)
    report=dict(device=name,checkpoint=str(Path(checkpoint).resolve()),checkpoint_sha256=digest(checkpoint),
        grid=reference.config.grid,cells=len(reference.ids),dt=reference.config.dt,steps=steps,
        kernel_benchmark=rows,cpu_single_thread_kernel_seconds=float(np.median(cpu_times[1:])),
        trajectory_maximum_errors=maxima,
        hybrid_step_seconds={k:float(np.median(v)) for k,v in timings.items()},
        scope='GPU mechanics only; polarity, contacts, chemistry and reductions stay on CPU. CUDA event timings exclude allocation/transfers; host call includes allocation, transfers, warmup and one measured kernel pair. Repeated kernel inputs are frozen; separate coupled trajectories advance normally. All host timings collected under three active production workers; not comparable to idle baseline as end-to-end speedup.',
        compiler=subprocess.check_output(['nvcc','--version'],text=True),
        sources={str(p.resolve()):digest(p) for p in [Path(__file__),Path(__file__).with_name('cuda_mechanics.cu'),Path(__file__).with_name('native_mechanics.py'),Path(__file__).with_name('native_mechanics.cpp')]})
    write_json(output/'comparison.json',report);print(json.dumps(report,indent=2),flush=True)
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--checkpoint',type=Path,default=Path('outputs/cell-exchange-response-moving/unexchanged/unexchanged_control/latest_state.npz'))
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--steps',type=int,default=48)
    parser.add_argument('--repeats',type=int,default=100)
    a=parser.parse_args()
    if min(a.steps,a.repeats)<1:parser.error('Positive steps/repeats required')
    run(a.checkpoint,a.output,a.steps,a.repeats)
