"""Optional second backend: fused surface polarity cue plus fast mechanics.

Original fast backend is unchanged. Float64 cue accumulation changes reduction
order, so this backend must pass separate numerical agreement gates.
"""
import ctypes
import hashlib
from pathlib import Path
import platform
import subprocess
import numpy as np
from .fast_mechanics import FastAttributeSimulation
from .model import occupancy
from .signaling import normalized_graph,integrate
from .polarity import evolve

_LIBRARY=None


def compiled_exposure(phi,dx,h=None):
    global _LIBRARY
    if phi.dtype!=np.float32 or phi.ndim!=4 or len(set(phi.shape[1:]))!=1 or phi.shape[1]<2 or not np.isfinite(dx) or dx<=0:
        raise ValueError('Requires cubic float32 cell fields and positive spacing')
    if _LIBRARY is None:
        source=Path(__file__).with_suffix('.c');key=hashlib.sha256(source.read_bytes()+platform.machine().encode()).hexdigest()[:20]
        folder=Path(__file__).resolve().parent.parent/'outputs/mechanics-kernel-cache';folder.mkdir(parents=True,exist_ok=True)
        binary=folder/f'cue-{key}.so'
        if not binary.exists():
            import os
            tmp=binary.with_suffix(f'.{os.getpid()}.tmp.so')
            subprocess.run(['cc','-O3','-std=c99','-shared','-fPIC','-ffp-contract=off',str(source),'-lm','-o',str(tmp)],check=True,capture_output=True);tmp.replace(binary)
        lib=ctypes.CDLL(str(binary));f=lib.surface_cue
        array=np.ctypeslib.ndpointer(dtype=np.float32,flags='C_CONTIGUOUS')
        f.argtypes=[array,array,array,array,np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS'),ctypes.c_int,ctypes.c_int,ctypes.c_float];f.restype=None;_LIBRARY=lib
    phi=np.ascontiguousarray(phi);h=occupancy(phi) if h is None else np.ascontiguousarray(h)
    if h.shape!=phi.shape or h.dtype!=np.float32:raise ValueError('Occupancy mismatch')
    occupied=h.sum(axis=0);shell_sums=(phi*(1-phi)).sum(axis=(1,2,3))
    result=np.empty((len(phi),3),dtype=np.float64)
    _LIBRARY.surface_cue(phi,h,occupied,shell_sums,result,len(phi),phi.shape[1],dx)
    return result


class FastPolaritySimulation(FastAttributeSimulation):
    def step(self,**kwargs):
        self._fast_cache={}
        try:return self._advance(**kwargs)
        finally:self._fast_cache=None

    def _advance(self, *, prescribed_signals=None):
        """Advance once, optionally clamping regulators for a forced-response test.

        A prescribed (activator, inhibitor) pair replaces graph reaction and
        transport for this step only. Fate, polarity, and mechanics still evolve.
        This externally maintained chemical pattern is not spontaneous signaling.
        """
        if prescribed_signals is not None:
            if not self.config.signaling:
                raise ValueError("prescribed signals require signaling-enabled downstream coupling")
            if len(prescribed_signals) != 2:
                raise ValueError("provide activator and inhibitor arrays")
            prescribed = tuple(np.asarray(x, dtype=float) for x in prescribed_signals)
            if any(x.shape != self.activator.shape or not np.isfinite(x).all() or np.any(x <= 0)
                   for x in prescribed):
                raise ValueError("prescribed signals must be positive finite arrays matching cells")
        contact, exposure = self.contacts()
        c = self.config
        graph = self.signaling_graph(contact)
        if prescribed_signals is not None:
            self.activator, self.inhibitor = (x.copy() for x in prescribed)
        elif c.signaling:
            self.activator, self.inhibitor = integrate(self.activator, self.inhibitor, graph, c.dt,
                                                       c.signal_beta, c.signal_da, c.signal_dh)
        if c.polarity_enabled:
            cue = compiled_exposure(self.phi, self.dx, self._geometry()['h'])
            self.polarity = evolve(self.polarity, cue, self.activator, normalized_graph(contact, c.graph_contact_cutoff if c.polarity_contact_cutoff == -1 else c.polarity_contact_cutoff).delta, c.dt,
                                   c.polarity_rate, c.polarity_alignment, c.polarity_decay)
        self.update_fate(contact, exposure)
        old_volumes = self.volumes()
        self.mechanical_step()
        if c.signaling and c.signal_transport == "conservative" and prescribed_signals is None:
            # Lagrangian compartments: mechanics redistributes concentration,
            # never creates molecular amount. Reactions were advanced above.
            dilution = old_volumes / self.volumes()
            self.activator *= dilution
            self.inhibitor *= dilution
        self.step_number += 1
        self.time = self.step_number * self.config.dt
        self._finish_divisions()
        while len(self.phi) + len(self.divisions) < self.config.max_cells and np.any(self.due <= self.time):
            self.divide(int(np.argmin(self.due)))
        if not np.all(np.isfinite(self.phi)) or not np.all(np.isfinite(self.fate)):
            raise FloatingPointError("non-finite state; reduce time step")
