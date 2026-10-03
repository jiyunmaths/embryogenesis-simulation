"""Opt-in fused C++/OpenMP mechanics; existing reference backends unchanged."""
import ctypes
import hashlib
import os
from pathlib import Path
import platform
import subprocess
import numpy as np
from .fast_mechanics import FastAttributeSimulation
from .signaling import normalized_graph,integrate
from .polarity import evolve

_LIB=None


def kernel():
    global _LIB
    if _LIB is None:
        source=Path(__file__).with_suffix('.cpp');key=hashlib.sha256(source.read_bytes()+platform.machine().encode()+b'O3-openmp-contract-off-no-errno').hexdigest()[:20]
        folder=source.parent.parent/'outputs/mechanics-kernel-cache';folder.mkdir(parents=True,exist_ok=True);binary=folder/f'mechanics-{key}.so'
        if not binary.exists():
            tmp=binary.with_suffix(f'.{os.getpid()}.tmp.so');subprocess.run(['c++','-O3','-std=c++17','-shared','-fPIC','-fopenmp','-ffp-contract=off','-fno-math-errno',str(source),'-o',str(tmp)],check=True,capture_output=True);tmp.replace(binary)
        lib=ctypes.CDLL(str(binary));f=lib.mechanical_update
        a=np.ctypeslib.ndpointer(dtype=np.float32,flags='C_CONTIGUOUS');d=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
        f.argtypes=[a,a,d,d,d,d,d,a,ctypes.c_int,ctypes.c_int,*([ctypes.c_double]*6),ctypes.c_int];f.restype=None
        cue=lib.surface_weights;cue.argtypes=[a,a,a,a,ctypes.c_int,ctypes.c_int,ctypes.c_float,ctypes.c_int];cue.restype=None
        _LIB=lib
    return _LIB.mechanical_update


class NativeSimulation(FastAttributeSimulation):
    """Opt-in backend. Set native_threads for each worker (not checkpointed).

    Limit BLAS threads to one when using multiple OpenMP threads or workers.
    """
    native_threads=1

    def step(self,**kwargs):
        self._validate_native()
        self._fast_cache={}
        try:return self._advance(**kwargs)
        finally:self._fast_cache=None

    def contacts(self):
        cache=self._geometry()
        if cache is None:return super().contacts()
        shell=self.phi*(1-self.phi);cache['shell']=shell
        flat=shell.reshape(len(shell),-1);contacts=(flat@flat.T).astype(float)*self.dx**3
        np.fill_diagonal(contacts,0)
        h=cache['h'];occupied=h.sum(axis=0);cache['occupied']=occupied
        blocked=np.clip(2*(occupied[None]-h),0,1)
        sums=shell.sum(axis=(1,2,3));cache['shell_sums']=sums
        exposure=1-(shell*blocked).sum(axis=(1,2,3))/np.maximum(sums,1e-12)
        return contacts,exposure

    def native_exposure(self):
        self._validate_native()
        cache=self._geometry()
        if cache is None:
            self._fast_cache={}
            try:return self.native_exposure()
            finally:self._fast_cache=None
        if 'shell' not in cache:self.contacts()
        field=self.phi
        out=np.empty((len(field),3,*field.shape[1:]),dtype=np.float32)
        kernel();_LIB.surface_weights(field,cache['h'],cache['occupied'],out,len(field),self.config.grid,self.dx,self.native_threads)
        # Retain NumPy float32 reduction order and normalization of the reference.
        cue=out.sum(axis=(2,3,4))/np.maximum(cache['shell_sums'],1e-12)[:,None]
        cue[np.linalg.norm(cue,axis=1)<1e-6]=0.
        return cue

    def _validate_native(self):
        if not isinstance(self.native_threads,int) or not 1<=self.native_threads<=16:
            raise ValueError('Use 1–16 native threads')
        if self.phi.dtype!=np.float32 or not self.phi.flags.c_contiguous:
            raise ValueError('Expected contiguous float32 phase fields')
        if self.phi.shape != (len(self.ids),self.config.grid,self.config.grid,self.config.grid) or self.config.grid<2:
            raise ValueError('Expected cubic phase fields matching configured grid')

    def mechanical_step(self):
        c=self.config
        # Preserve exact original stencils for nonpolar cells and cytokinesis.
        if self.divisions or not c.feedback or not c.polarity_enabled or np.any(np.linalg.norm(self.polarity,axis=1)<=1e-10):return super().mechanical_step()
        phi=self.phi
        if phi.dtype!=np.float32 or not phi.flags.c_contiguous:raise ValueError('Expected contiguous float32 phase fields')
        if not 1<=self.native_threads<=16:raise ValueError('Use 1–16 native threads')
        cache=self._geometry();shell=cache.get('shell') if cache is not None else None
        if shell is None:shell=phi*(1-phi)
        shell2=shell**2;tension,adhesion=self.material_coefficients()
        attract=np.ascontiguousarray((adhesion@shell2.reshape(len(phi),-1)).reshape(phi.shape),dtype=np.float64)
        exclude=np.sum(phi**2,axis=0)[None]-phi**2
        vf=c.volume_stiffness*(self.target-self.volumes())/self.target
        updated=np.empty_like(phi);self.volume_projection_max=0.
        kernel()(phi,exclude,attract,np.ascontiguousarray(self.centers()),np.ascontiguousarray(self.polarity),tension,np.ascontiguousarray(vf),updated,len(phi),c.grid,self.dx,c.extent,c.interface_width,c.polarity_tension,c.repulsion,c.dt,self.native_threads)
        self.clipped_fraction=float(np.mean((updated<0)|(updated>1)));self.phi=np.clip(updated,0,1)

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
            cue = self.native_exposure()
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
