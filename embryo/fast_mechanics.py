"""Opt-in compiled polar flux and step-local geometry reuse.

Default production classes remain unchanged. Requires a C99 compiler on first
use; built library is cached by source/platform hash under outputs/.
"""
import ctypes
import hashlib
from pathlib import Path
import platform
import subprocess
import numpy as np
from .attribute_development import AttributeSimulation
from .model import occupancy
from scipy.ndimage import laplace

_LIBRARY=None


def kernel():
    global _LIBRARY
    if _LIBRARY is None:
        source=Path(__file__).with_suffix('.c')
        key=hashlib.sha256(source.read_bytes()+platform.machine().encode()+b'O3-no-fp-contract-v1').hexdigest()[:20]
        folder=Path(__file__).resolve().parent.parent/'outputs'/'mechanics-kernel-cache';folder.mkdir(parents=True,exist_ok=True)
        binary=folder/f'polar-{key}.so'
        if not binary.exists():
            import os
            temp=binary.with_suffix(f'.{os.getpid()}.tmp.so')
            subprocess.run(['cc','-O3','-std=c99','-shared','-fPIC','-ffp-contract=off',str(source),'-lm','-o',str(temp)],check=True,capture_output=True)
            temp.replace(binary)
        library=ctypes.CDLL(str(binary));f=library.polar_flux
        f.argtypes=[np.ctypeslib.ndpointer(dtype=np.float32,flags='C_CONTIGUOUS'),np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS'),np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS'),ctypes.c_int,ctypes.c_double,ctypes.c_double,np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS'),np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS'),ctypes.c_double,ctypes.c_double,ctypes.c_double]
        f.restype=None;_LIBRARY=library
    return _LIBRARY.polar_flux


def polar_fields(field,dx,extent,center,p,base,contrast,width):
    field=np.ascontiguousarray(field,dtype=np.float32)
    if field.ndim!=3 or len(set(field.shape))!=1:raise ValueError('Requires cubic cell field')
    gamma=np.empty(field.shape,dtype=np.float64);flux=np.empty_like(gamma)
    kernel()(field,gamma,flux,field.shape[0],dx,extent,np.ascontiguousarray(center,dtype=np.float64),np.ascontiguousarray(p,dtype=np.float64),base,contrast,width)
    return gamma,flux


class FastAttributeSimulation(AttributeSimulation):
    """Same integrator with opt-in polar kernel; original cleavage mechanics."""
    def step(self,**kwargs):
        self._fast_cache={}
        try:return super().step(**kwargs)
        finally:self._fast_cache=None

    def _geometry(self):
        cache=getattr(self,'_fast_cache',None)
        if cache is None:return None
        if cache.get('phi') is not self.phi:
            cache.clear();cache['phi']=self.phi
        if 'h' not in cache:cache['h']=occupancy(self.phi)
        return cache

    def volumes(self):
        cache=self._geometry()
        if cache is None:return super().volumes()
        if 'volumes' not in cache:cache['volumes']=cache['h'].sum(axis=(1,2,3),dtype=np.float64)*self.dx**3
        return cache['volumes']

    def centers(self):
        cache=self._geometry()
        if cache is None:return super().centers()
        if 'centers' not in cache:cache['centers']=np.einsum('nijk,dijk->nd',cache['h'],self.xyz)*self.dx**3/self.volumes()[:,None]
        return cache['centers']

    def mechanical_step(self):
        # Do not change cytokinesis, ring forces or projection arithmetic.
        if self.divisions:
            super().mechanical_step()
            if getattr(self,'_fast_cache',None) is not None:self._fast_cache.clear()
            return
        c=self.config;phi=self.phi
        shell2=(phi*(1-phi))**2;derivative=2*phi*(1-phi)*(1-2*phi)
        tension,adhesion=self.material_coefficients()
        attract=(adhesion@shell2.reshape(len(phi),-1)).reshape(phi.shape)
        exclude=np.sum(phi**2,axis=0)[None]-phi**2
        volumes=self.volumes();centers=self.centers()
        volume_force=c.volume_stiffness*(self.target-volumes)/self.target
        updated=np.empty_like(phi);self.volume_projection_max=0.
        for i in range(len(phi)):
            if c.feedback and c.polarity_enabled and np.linalg.norm(self.polarity[i])>1e-10:
                gamma,flux=polar_fields(phi[i],self.dx,c.extent,centers[i],self.polarity[i],tension[i],c.polarity_tension,c.interface_width)
                force=c.interface_width**2*flux-gamma*derivative[i]
            else:
                diffusion=laplace(phi[i],mode='nearest')/self.dx**2
                force=tension[i]*(c.interface_width**2*diffusion-derivative[i])
            force+=volume_force[i]*6*phi[i]*(1-phi[i])
            force-=c.repulsion*phi[i]*exclude[i]
            force+=derivative[i]*attract[i]
            updated[i]=phi[i]+c.dt*force
        self.clipped_fraction=float(np.mean((updated<0)|(updated>1)));self.phi=np.clip(updated,0,1)
