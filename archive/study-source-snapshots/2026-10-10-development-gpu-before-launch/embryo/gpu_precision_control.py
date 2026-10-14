"""Opt-in numerical precision controls, separate from accepted GPU physics.

Carry sub-float32 phase updates in float64; independently accumulate contacts
in float64. This is not an all-float64 backend or scientific acceptance.
"""
import ctypes
import hashlib
import os
from pathlib import Path
import subprocess

import numpy as np
import torch

from .gpu_backend import GpuSimulation, call, gm_step, transport_matrices

ARMS = ('baseline', 'phase_carry', 'contact64', 'both')
_LIBRARY = None


def library():
    global _LIBRARY
    if _LIBRARY is None:
        source = Path(__file__).with_name('gpu_precision_kernels.cu')
        files = [source, *[source.with_name(f) for f in ('gpu_kernels.cu', 'cuda_spatial_bench.cu', 'cuda_mechanics.cu')]]
        flags = ['-O3', '-std=c++17', '-arch=sm_61', '--fmad=false', '--ftz=false',
            '--prec-div=true', '--prec-sqrt=true', '-ccbin', '/usr/bin/g++-12', '-shared', '-Xcompiler', '-fPIC']
        key = hashlib.sha256(b''.join(f.read_bytes() for f in files)+repr(flags).encode()+
            subprocess.check_output(['nvcc', '--version'])).hexdigest()[:20]
        folder = source.parent.parent/'outputs/cuda-kernel-cache'; folder.mkdir(parents=True, exist_ok=True)
        binary = folder/f'precision-control-{key}.so'
        if not binary.exists():
            temporary = binary.with_suffix(f'.{os.getpid()}.tmp.so')
            subprocess.run(['nvcc', *flags, str(source), '-lcublas', '-o', str(temporary)], check=True)
            temporary.replace(binary)
        lib = ctypes.CDLL(str(binary)); ptr, integer, real = ctypes.c_void_p, ctypes.c_int, ctypes.c_double
        lib.precision_force.argtypes = [ptr]*12+[integer]*3+[real]*6+[ptr, integer]
        lib.precision_force.restype = integer; lib.cuda_mechanics_error.restype = ctypes.c_char_p
        _LIBRARY = lib
    return _LIBRARY


def force(tensors, scalars):
    if len(tensors) != 12 or not all(t.is_cuda and t.is_contiguous() for t in tensors):
        raise ValueError('Twelve contiguous CUDA buffers required')
    device = tensors[0].device
    if any(t.device != device for t in tensors): raise ValueError('One CUDA device required')
    f, d, i = torch.float32, torch.float64, torch.int32
    if any(t.dtype != kind for t, kind in zip(tensors, (f, f, d, d, d, d, d, d, f, i, d, i))):
        raise ValueError('Precision-control buffer dtype mismatch')
    cells, n, carry = scalars[:3]; field = (cells, n, n, n)
    shapes = (field, field, field, (cells, 3), (cells, 3), (cells,), (cells,), field, field, (2,), field, (2,))
    if carry not in (0, 1) or any(tuple(t.shape) != shape for t, shape in zip(tensors, shapes)):
        raise ValueError('Precision-control buffer shape or carry flag mismatch')
    lib = library()
    code = lib.precision_force(*[t.data_ptr() for t in tensors], *scalars,
        torch.cuda.current_stream(device).cuda_stream, device.index)
    if code: raise RuntimeError(lib.cuda_mechanics_error().decode())


class PrecisionSimulation(GpuSimulation):
    def __init__(self, host, arm='baseline', device='cuda:0'):
        if arm not in ARMS: raise ValueError('Unknown numerical control')
        self.precision_arm = arm
        super().__init__(host, device)
        self.phase_carry = torch.zeros_like(self.phi, dtype=torch.float64)
        self.rounding = torch.zeros(2, dtype=torch.int32, device=self.device)

    def checkpoint(self, path):
        """Save the carry along with the visible state; never reset it on resume."""
        path = Path(path); temporary = path.with_name(path.stem+'.precision.tmp.npz')
        super().checkpoint(temporary)
        with np.load(temporary, allow_pickle=False) as z:
            payload = {key: z[key].copy() for key in z.files}
        payload.update(precision_arm=np.array(self.precision_arm),
            phase_carry=self.phase_carry.cpu().numpy(), rounding=self.rounding.cpu().numpy())
        np.savez_compressed(temporary, **payload); temporary.replace(path)

    @classmethod
    def restore(cls, path, device='cuda:0'):
        from .attribute_development import AttributeSimulation
        with np.load(path, allow_pickle=False) as z:
            arm = str(z['precision_arm']); carry = z['phase_carry'].copy(); counts = z['rounding'].copy()
        host = AttributeSimulation.restore(path)
        if (carry.dtype != np.float64 or carry.shape != host.phi.shape or not np.isfinite(carry).all()
                or counts.dtype != np.int32 or counts.shape != (2,) or np.any(counts < 0)
                or (arm in ('baseline', 'contact64') and np.any(carry != 0))):
            raise ValueError('Invalid precision checkpoint carry or counters')
        sim = cls(host, arm, device)
        sim.phase_carry.copy_(torch.from_numpy(carry).to(sim.device))
        sim.rounding.copy_(torch.from_numpy(counts).to(sim.device))
        return sim

    def matrices(self):
        if self.precision_arm not in ('contact64', 'both'):
            return super().matrices()
        c = self.config; flat = self.shell.double().flatten(1)
        contact = (flat@flat.T)*self.dx**3; contact.fill_diagonal_(0)
        pc = c.graph_contact_cutoff if c.polarity_contact_cutoff == -1 else c.polarity_contact_cutoff
        delta, normalized, conductance = transport_matrices(contact, self.geometry[:, 0].contiguous(),
            self.geometry[:, 1:4].contiguous(), c.interface_width, c.graph_contact_cutoff, pc)
        return contact, delta, normalized, conductance

    @torch.inference_mode()
    def step(self):
        # Match the accepted splitting/order, changing only two named rounding paths.
        c = self.config; _, delta, normalized, _ = self.matrices()
        self.activator, self.inhibitor = gm_step(self.activator, self.inhibitor, delta, c.dt,
            c.signal_beta, c.signal_da, c.signal_dh)
        aligned = (normalized@self.polarity).contiguous(); p = torch.empty_like(self.polarity)
        call('gpu_polarity', [self.polarity, self.geometry[:, 4:].contiguous(), self.activator, aligned, p],
            (len(self.ids), c.dt, c.polarity_rate, c.polarity_alignment, c.polarity_decay))
        self.polarity = p; response = torch.tanh(self.activator-1)
        tensions = c.surface_tension*(1+c.fate_tension*response)
        adhesion = c.adhesion*(1+c.fate_adhesion*response[:, None]*response[None, :]); adhesion.fill_diagonal_(0)
        attraction = (adhesion@self.shell2.flatten(1)).reshape(self.phi.shape)
        old_volume = self.geometry[:, 0].clone(); centers = self.geometry[:, 1:4].contiguous()
        vf = c.volume_stiffness*(self.target-old_volume)/self.target
        amounts = torch.stack([old_volume@self.activator, old_volume@self.inhibitor])
        force([self.phi, self.excluded, attraction, centers, self.polarity, tensions, vf,
            self.gamma, self.updated, self.quality, self.phase_carry, self.rounding],
            (len(self.ids), c.grid, int(self.precision_arm in ('phase_carry', 'both')), self.dx,
             c.extent, c.interface_width, c.polarity_tension, c.repulsion, c.dt))
        self.phi, self.updated = self.updated, self.phi; self.refresh()
        new_volume = self.geometry[:, 0]
        self.activator *= old_volume/new_volume; self.inhibitor *= old_volume/new_volume
        new_amounts = torch.stack([new_volume@self.activator, new_volume@self.inhibitor])
        self._dilution_error = ((new_amounts/amounts)-1).abs().max()
        self.step_number += 1; self.time = self.step_number*c.dt

    def precision_diagnostics(self):
        counts = self.rounding.cpu().numpy(); active, stuck = map(int, counts)
        return dict(active_force_voxels=active, uncarried_stuck_voxels=stuck,
            uncarried_stuck_fraction=stuck/max(active, 1),
            phase_carry_abs_max=float(self.phase_carry.abs().max()),
            phase_carry_l2=float(torch.linalg.vector_norm(self.phase_carry)))
