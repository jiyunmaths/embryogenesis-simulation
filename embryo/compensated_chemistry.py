"""Isolated float64 SSP-RK2 with retained accumulation rounding error.

No historical kernel imports or monkeypatches. Reactions and the conservative
matrix are unchanged; only addition of the RK2 increment is compensated.
The stored correction is the error to SUBTRACT from the next increment.
"""
import math
import numpy as np
import torch

METHOD = 'float64 SSP-RK2 increment with Kahan accumulation v1'


def step_count(delta, dt, beta=2., da=.02, db=.55, reaction=True):
    if (isinstance(dt, (bool, np.bool_)) or not np.isfinite([dt,beta,da,db]).all()
            or dt<=0 or beta<=0 or min(da,db)<0):
        raise ValueError('Finite positive timestep/beta and nonnegative diffusion required')
    if isinstance(delta, torch.Tensor):
        if delta.dtype != torch.float64 or delta.ndim!=2 or delta.shape[0]!=delta.shape[1]:
            raise ValueError('Square float64 transport matrix required')
        exit_rate=float((-delta.diagonal()).max().item())
    else:
        if delta.dtype != np.float64 or delta.ndim!=2 or delta.shape[0]!=delta.shape[1]:
            raise ValueError('Square float64 transport matrix required')
        exit_rate=float(np.max(-np.diag(delta)))
    if not np.isfinite(exit_rate) or exit_rate<0:raise ValueError('Invalid transport exit rate')
    loss=max((1. if reaction else 0.)+da*exit_rate,(beta if reaction else 0.)+db*exit_rate)
    required=dt*loss/.2
    if not np.isfinite(required):raise ValueError('Nonfinite positivity substeps')
    return max(1,math.ceil(required))


def numpy_add(values, correction, increment):
    """Return accumulated high part and its signed rounding correction."""
    adjusted=increment-correction
    updated=values+adjusted
    return updated,(updated-values)-adjusted


def torch_add(values, correction, increment):
    adjusted=increment-correction
    updated=values+adjusted
    return updated,(updated-values)-adjusted


def _numpy_positive(value, label):
    if not np.isfinite(value).all() or np.any(value<=0):
        raise FloatingPointError(label+' lost positivity/finiteness; no clipping')


def _torch_positive(value, label):
    if not bool(torch.isfinite(value).all() & (value>0).all()):
        raise FloatingPointError(label+' lost positivity/finiteness; no clipping')


class NumpyStepper:
    """Fixed-geometry NumPy implementation, with optional pure-exchange assay."""
    def __init__(self,delta,dt,beta=2.,da=.02,db=.55,reaction=True):
        self.count=step_count(delta,dt,beta,da,db,reaction)
        if not np.isfinite(delta).all():raise ValueError('Nonfinite transport matrix')
        self.delta=delta.copy();self.dt=float(dt)/self.count
        self.beta,self.da,self.db=map(float,(beta,da,db));self.reaction=bool(reaction)

    def rhs(self,values):
        a,b=values
        fa=self.da*(self.delta@a);fb=self.db*(self.delta@b)
        if self.reaction:
            fa=a*a/b-a+fa;fb=self.beta*(a*a-b)+fb
        return np.stack((fa,fb))

    def __call__(self,values,correction):
        if (values.dtype!=np.float64 or correction.dtype!=np.float64
                or values.shape!=correction.shape or values.ndim not in (2,3)
                or values.shape[:2]!=(2,self.delta.shape[0])):
            raise ValueError('Matching float64 species/state/correction shapes required')
        _numpy_positive(values,'Initial state')
        if not np.isfinite(correction).all():raise FloatingPointError('Nonfinite correction')
        with np.errstate(over='raise',divide='raise',invalid='raise'):
            for _ in range(self.count):
                first_rhs=self.rhs(values);first=values+self.dt*first_rhs
                _numpy_positive(first,'First Euler stage')
                second_rhs=self.rhs(first);second=first+self.dt*second_rhs
                _numpy_positive(second,'Second Euler stage')
                increment=(.5*self.dt)*(first_rhs+second_rhs)
                values,correction=numpy_add(values,correction,increment)
                _numpy_positive(values,'Compensated final stage')
                if not np.isfinite(correction).all():raise FloatingPointError('Nonfinite correction')
        return values,correction


class TorchStepper:
    """Resident float64 PyTorch arrays/matmul; no CPU scientific fallback."""
    def __init__(self,delta,dt,beta=2.,da=.02,db=.55,reaction=True):
        if not isinstance(delta,torch.Tensor):raise ValueError('PyTorch matrix required')
        self.count=step_count(delta,dt,beta,da,db,reaction)
        if not bool(torch.isfinite(delta).all()):raise ValueError('Nonfinite transport matrix')
        self.delta=delta.detach().clone();self.dt=float(dt)/self.count
        self.beta,self.da,self.db=map(float,(beta,da,db));self.reaction=bool(reaction)

    def rhs(self,values):
        a,b=values
        fa=self.da*(self.delta@a);fb=self.db*(self.delta@b)
        if self.reaction:
            fa=a*a/b-a+fa;fb=self.beta*(a*a-b)+fb
        return torch.stack((fa,fb))

    def __call__(self,values,correction):
        if (values.dtype!=torch.float64 or correction.dtype!=torch.float64
                or values.shape!=correction.shape or values.ndim not in (2,3)
                or values.shape[:2]!=(2,self.delta.shape[0])
                or values.device!=self.delta.device or correction.device!=values.device):
            raise ValueError('Matching float64 species/state/correction shapes and devices required')
        _torch_positive(values,'Initial state')
        if not bool(torch.isfinite(correction).all()):raise FloatingPointError('Nonfinite correction')
        for _ in range(self.count):
            first_rhs=self.rhs(values);first=values+self.dt*first_rhs
            _torch_positive(first,'First Euler stage')
            second_rhs=self.rhs(first);second=first+self.dt*second_rhs
            _torch_positive(second,'Second Euler stage')
            increment=(.5*self.dt)*(first_rhs+second_rhs)
            values,correction=torch_add(values,correction,increment)
            _torch_positive(values,'Compensated final stage')
            if not bool(torch.isfinite(correction).all()):raise FloatingPointError('Nonfinite correction')
        return values,correction
