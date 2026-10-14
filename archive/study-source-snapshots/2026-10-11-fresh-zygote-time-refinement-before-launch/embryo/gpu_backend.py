"""Opt-in resident PyTorch/custom-CUDA backend for mature attribute embryos.

No cleavage, fate switch, or silent CPU fallback. Spatial fields stay on GPU;
only small diagnostics and explicit standard checkpoints cross to the host.
"""
import ctypes
import hashlib
import math
import os
from pathlib import Path
import subprocess

import numpy as np
import torch

from .attribute_development import AttributeSimulation

_LIBRARY = None


def library():
    global _LIBRARY
    if _LIBRARY is None:
        source = Path(__file__).with_name('gpu_kernels.cu')
        dependencies = [source, source.with_name('cuda_spatial_bench.cu'),
                        source.with_name('cuda_mechanics.cu')]
        flags = ['-O3', '-std=c++17', '-arch=sm_61', '--fmad=false', '--ftz=false',
                 '--prec-div=true', '--prec-sqrt=true', '-ccbin', '/usr/bin/g++-12',
                 '-shared', '-Xcompiler', '-fPIC']
        version = subprocess.check_output(['nvcc', '--version'])
        key = hashlib.sha256(b''.join(f.read_bytes() for f in dependencies)
                             + repr(flags).encode() + version).hexdigest()[:20]
        folder = source.parent.parent / 'outputs/cuda-kernel-cache'
        folder.mkdir(parents=True, exist_ok=True)
        binary = folder / f'resident-{key}.so'
        if not binary.exists():
            temporary = binary.with_suffix(f'.{os.getpid()}.tmp.so')
            subprocess.run(['nvcc', *flags, str(source), '-lcublas', '-o', str(temporary)], check=True)
            temporary.replace(binary)
        lib = ctypes.CDLL(str(binary))
        ptr, integer, real = ctypes.c_void_p, ctypes.c_int, ctypes.c_double
        signatures = {
            'gpu_arrays': [ptr]*6 + [integer]*2 + [ptr, integer],
            'gpu_geometry': [ptr]*6 + [integer]*2 + [real]*2 + [ptr, integer],
            'gpu_force': [ptr]*10 + [integer]*2 + [real]*6 + [ptr, integer],
            'gpu_polarity': [ptr]*5 + [integer] + [real]*4 + [ptr, integer],
            'gpu_union': [ptr]*2 + [integer] + [real]*2 + [ptr, integer],
        }
        for name, signature in signatures.items():
            function = getattr(lib, name)
            function.argtypes, function.restype = signature, integer
        lib.cuda_mechanics_error.restype = ctypes.c_char_p
        _LIBRARY = lib
    return _LIBRARY


def call(name, tensors, scalars=()):
    """Launch on PyTorch's current stream; tensors own every device allocation."""
    if not tensors or not all(t.is_cuda and t.is_contiguous() for t in tensors):
        raise ValueError('Custom kernels require contiguous CUDA tensors')
    device = tensors[0].device
    if any(t.device != device for t in tensors):
        raise ValueError('All buffers must be on the same CUDA device')
    f, d, i = torch.float32, torch.float64, torch.int32
    specifications = {
        'gpu_arrays': (f, f, f, d, f, f),
        'gpu_geometry': (f, f, f, f, d, d),
        'gpu_force': (f, f, d, d, d, d, d, d, f, i),
        'gpu_polarity': (d, d, d, d, d),
        'gpu_union': (f, d),
    }
    expected = specifications[name]
    if len(tensors) != len(expected) or any(t.dtype != dtype for t, dtype in zip(tensors, expected)):
        raise ValueError('Kernel buffer dtype/count mismatch')
    if name in ('gpu_arrays', 'gpu_geometry', 'gpu_force'):
        cells, n = scalars[:2]
        field, grid = (cells, n, n, n), (n, n, n)
        chunks = (n**3+255)//256
        shapes = {
            'gpu_arrays': (field, field, field, field, grid, field),
            'gpu_geometry': (field, field, field, grid, (cells, 8, chunks), (cells, 7)),
            'gpu_force': (field, field, field, (cells, 3), (cells, 3), (cells,), (cells,), field, field, (2,)),
        }[name]
    elif name == 'gpu_polarity':
        cells = scalars[0]
        shapes = ((cells, 3), (cells, 3), (cells,), (cells, 3), (cells, 3))
    else:
        n = scalars[0]
        shapes = ((n, n, n), (10, (n**3+255)//256))
    if any(tuple(t.shape) != shape for t, shape in zip(tensors, shapes)):
        raise ValueError('Kernel buffer shape mismatch')
    lib = library()
    code = getattr(lib, name)(*[t.data_ptr() for t in tensors], *scalars,
                             torch.cuda.current_stream(device).cuda_stream, device.index)
    if code:
        raise RuntimeError(lib.cuda_mechanics_error().decode())


def transport_matrices(contacts, volumes, centers, width, cutoff, polarity_cutoff):
    """Same geometric conductance closure and separate normalized polarity graph."""
    w = .5*(contacts+contacts.T)
    w = w.clone()
    w.fill_diagonal_(0)
    w = torch.where(w < cutoff*w.max(), 0., w)
    distance = torch.linalg.vector_norm(centers[:, None]-centers[None], dim=-1)
    active = w > 0
    # Invalid touching centers remain nonfinite and are caught by the step audit.
    safe_distance = torch.where(active, distance, torch.ones_like(distance))
    conductance = 6*math.sqrt(2)*w/width/safe_distance
    degree = conductance.sum(dim=1)
    delta = (conductance-torch.diag(degree))/volumes[:, None]
    p = contacts.clone()
    p.fill_diagonal_(0)
    threshold = torch.maximum(polarity_cutoff*p.max(), p.new_tensor(1e-12))
    p = torch.where(p < threshold, 0., p)
    pd = p.sum(dim=1)
    normalized = p/pd.clamp_min(1e-300)[:, None]-torch.diag((pd > 0).double())
    return delta, normalized, conductance


def gm_step(a, h, delta, dt, beta, da, dh):
    """Reference SSP-RK2; positivity restriction uses the measured exit rate."""
    exit_rate = float((-delta.diagonal()).max())
    substeps = max(1, math.ceil(dt*max(1+da*exit_rate, beta+dh*exit_rate)/.2))
    step = dt/substeps

    def rhs(x, y):
        return x*x/y-x+da*(delta@x), beta*(x*x-y)+dh*(delta@y)

    for _ in range(substeps):
        fa, fh = rhs(a, h)
        first_a, first_h = a+step*fa, h+step*fh
        fa, fh = rhs(first_a, first_h)
        a = .5*a+.5*(first_a+step*fa)
        h = .5*h+.5*(first_h+step*fh)
    return a, h


class GpuSimulation:
    def __init__(self, host, device='cuda:0'):
        c = host.config
        if (host.divisions or len(host.ids) != c.max_cells or
                host.attribute_mode != 'direct' or not c.feedback or not c.polarity_enabled or
                not c.signaling or c.signal_transport != 'conservative' or c.differentiation or
                np.any(host.fate != 0) or np.any(np.linalg.norm(host.polarity, axis=1) <= 1e-10)):
            raise ValueError('GPU backend requires mature, polar, direct attribute state with conservative signaling')
        if host.phi.dtype != np.float32 or host.phi.shape != (len(host.ids), c.grid, c.grid, c.grid):
            raise ValueError('Expected cubic float32 phase fields')
        if not all(np.isfinite(x).all() for x in (host.phi, host.activator, host.inhibitor, host.polarity)):
            raise ValueError('Invalid starting state')
        if np.any(host.phi < 0) or np.any(host.phi > 1):
            raise ValueError('Phase fields must lie in [0,1]')
        if np.any(host.activator <= 0) or np.any(host.inhibitor <= 0):
            raise ValueError('Positive chemical concentrations required')
        self.device = torch.device(device)
        if self.device.index is None:
            self.device = torch.device('cuda', torch.cuda.current_device())
        if self.device.type != 'cuda':
            raise ValueError('CUDA device required')
        torch.backends.cuda.matmul.allow_tf32 = False
        self.host, self.config = host, c
        self.ids = host.ids.copy()
        self.time, self.step_number, self.dx = host.time, host.step_number, host.dx
        upload = lambda x: torch.from_numpy(np.ascontiguousarray(x)).to(self.device)
        self.phi, self.activator, self.inhibitor, self.polarity, self.target = [
            upload(x) for x in (host.phi, host.activator, host.inhibitor, host.polarity, host.target)]
        if any(t.dtype != torch.float64 for t in (self.activator, self.inhibitor, self.polarity, self.target)):
            raise ValueError('Chemical, polarity and target arrays must be float64')
        self.h = torch.empty_like(self.phi)
        self.shell = torch.empty_like(self.phi)
        self.shell2 = torch.empty_like(self.phi, dtype=torch.float64)
        self.occupied = torch.empty_like(self.phi[0])
        self.excluded = torch.empty_like(self.phi)
        self.gamma = torch.empty_like(self.shell2)
        self.updated = torch.empty_like(self.phi)
        chunks = (c.grid**3+255)//256
        self.partial = torch.empty((len(self.ids), 8, chunks), dtype=torch.float64, device=self.device)
        self.geometry = torch.empty((len(self.ids), 7), dtype=torch.float64, device=self.device)
        self.union_partial = torch.empty((10, chunks), dtype=torch.float64, device=self.device)
        self.quality = torch.zeros(2, dtype=torch.int32, device=self.device)
        self.clipped_fraction = 0.
        self.volume_projection_max = 0.
        self.dilution_amount_error = 0.
        self.refresh()

    @classmethod
    def restore(cls, path, device='cuda:0'):
        return cls(AttributeSimulation.restore(path), device)

    def refresh(self):
        call('gpu_arrays', [self.phi, self.h, self.shell, self.shell2, self.occupied, self.excluded],
             (len(self.ids), self.config.grid))
        call('gpu_geometry', [self.phi, self.h, self.shell, self.occupied, self.partial, self.geometry],
             (len(self.ids), self.config.grid, self.dx, self.config.extent))

    def matrices(self):
        c = self.config
        flat = self.shell.flatten(1)
        contact = (flat@flat.T).double()*self.dx**3
        contact.fill_diagonal_(0)
        volume = self.geometry[:, 0].contiguous()
        centers = self.geometry[:, 1:4].contiguous()
        pc = c.graph_contact_cutoff if c.polarity_contact_cutoff == -1 else c.polarity_contact_cutoff
        delta, normalized, conductance = transport_matrices(contact, volume, centers,
                                                           c.interface_width, c.graph_contact_cutoff, pc)
        return contact, delta, normalized, conductance

    @torch.inference_mode()
    def step(self):
        c = self.config
        _, delta, normalized, _ = self.matrices()
        self.activator, self.inhibitor = gm_step(self.activator, self.inhibitor, delta, c.dt,
                                               c.signal_beta, c.signal_da, c.signal_dh)
        aligned = (normalized@self.polarity).contiguous()
        p = torch.empty_like(self.polarity)
        call('gpu_polarity', [self.polarity, self.geometry[:, 4:].contiguous(), self.activator, aligned, p],
             (len(self.ids), c.dt, c.polarity_rate, c.polarity_alignment, c.polarity_decay))
        self.polarity = p
        response = torch.tanh(self.activator-1)
        tensions = c.surface_tension*(1+c.fate_tension*response)
        adhesion = c.adhesion*(1+c.fate_adhesion*response[:, None]*response[None, :])
        adhesion.fill_diagonal_(0)
        attraction = (adhesion@self.shell2.flatten(1)).reshape(self.phi.shape)
        old_volume = self.geometry[:, 0].clone()
        centers = self.geometry[:, 1:4].contiguous()
        vf = c.volume_stiffness*(self.target-old_volume)/self.target
        amounts = torch.stack([old_volume@self.activator, old_volume@self.inhibitor])
        call('gpu_force', [self.phi, self.excluded, attraction, centers, self.polarity, tensions, vf,
                           self.gamma, self.updated, self.quality],
             (len(self.ids), c.grid, self.dx, c.extent, c.interface_width,
              c.polarity_tension, c.repulsion, c.dt))
        self.phi, self.updated = self.updated, self.phi
        self.refresh()
        new_volume = self.geometry[:, 0]
        dilution = old_volume/new_volume
        self.activator *= dilution
        self.inhibitor *= dilution
        new_amounts = torch.stack([new_volume@self.activator, new_volume@self.inhibitor])
        self._dilution_error = ((new_amounts/amounts)-1).abs().max()
        self.step_number += 1
        self.time = self.step_number*c.dt

    def audit(self):
        """Small host transfer; every step checks quality and chemical validity."""
        volume = self.geometry[:, 0]
        finite = (torch.isfinite(self.geometry).all() & torch.isfinite(self.polarity).all()
                  & torch.isfinite(self.activator).all() & torch.isfinite(self.inhibitor).all())
        valid = finite & (self.activator > 0).all() & (self.inhibitor > 0).all() & (volume > 0).all()
        values = torch.stack([(volume/self.target-1).abs().max(),
                              ((3*volume/(4*math.pi))**(1/3)/self.dx).min(),
                              self.quality[0].double()/self.phi.numel(), self.quality[1].double(),
                              valid.double(), getattr(self, '_dilution_error', volume.new_tensor(0.))]).cpu().numpy()
        self.clipped_fraction = float(values[2])
        self.dilution_amount_error = float(values[5])
        if values[3] or not values[4]:
            raise FloatingPointError('Nonfinite or nonpositive GPU state')
        if values[0] >= .05 or values[1] < 4 or values[2] > 0:
            raise FloatingPointError('GPU numerical quality screen failed: '+str(values))
        return dict(max_volume_error=float(values[0]), min_radius=float(values[1]),
                    max_clipping=self.clipped_fraction, dilution_amount_error=self.dilution_amount_error)

    def observe(self, elapsed):
        c = self.config
        call('gpu_union', [self.occupied, self.union_partial], (c.grid, self.dx, c.extent))
        sums = self.union_partial.sum(dim=1).cpu().numpy()
        center = sums[1:4]/sums[0]
        covariance = np.array([[sums[4], sums[7], sums[8]], [sums[7], sums[5], sums[9]],
                               [sums[8], sums[9], sums[6]]])/sums[0]-np.outer(center, center)
        eigen = np.linalg.eigvalsh(covariance)
        edge = torch.stack([self.occupied[0].max(), self.occupied[-1].max(),
                            self.occupied[:, 0].max(), self.occupied[:, -1].max(),
                            self.occupied[:, :, 0].max(), self.occupied[:, :, -1].max()]).max().clamp_max(1)
        delta = self.matrices()[1].cpu().numpy()
        geometry = self.geometry.cpu().numpy()
        chemistry = torch.stack([self.activator, self.inhibitor]).cpu().numpy()
        return dict(elapsed=elapsed, time=self.time, ids=self.ids.tolist(), chemistry=chemistry.tolist(),
                    volumes=geometry[:, 0].tolist(), centers=geometry[:, 1:4].tolist(),
                    polarity=self.polarity.cpu().numpy().tolist(), delta=delta.tolist(),
                    axis_ratio=float(np.sqrt(eigen[-1]/max(eigen[0], 1e-12))),
                    boundary_occupancy=float(edge), log_activator_sd=float(np.std(np.log(chemistry[0]))))

    def to_cpu(self):
        host = self.host
        for name in ('phi', 'activator', 'inhibitor', 'polarity'):
            setattr(host, name, getattr(self, name).cpu().numpy().copy())
        host.time, host.step_number = self.time, self.step_number
        host.clipped_fraction, host.volume_projection_max = self.clipped_fraction, 0.
        return host

    def checkpoint(self, path):
        self.to_cpu().checkpoint(path)
