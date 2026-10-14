"""Isolated compensated chemical arithmetic with the existing CUDA mechanics."""
import json,math
from pathlib import Path
import numpy as np
import torch
from . import compensated_chemistry as chemical
from .gpu_backend import call
from .gpu_precision_control import PrecisionSimulation,force
from .attribute_development import AttributeSimulation
from .phase_carry_exchange_response import payload
from .resolution import _steps

METHOD='compensated SSP-RK2 and error-free dilution, phase-carry CUDA v1'
SPLITTER=float(2**27+1)


def numpy_scale(u,e,r):
    """Scale the represented concentration u-e; retain product rounding."""
    with np.errstate(over='raise',invalid='raise',divide='raise'):
        product=u*r;cu=SPLITTER*u;cr=SPLITTER*r
        uh=cu-(cu-u);ul=u-uh;rh=cr-(cr-r);rl=r-rh
        error=((uh*rh-product)+uh*rl+ul*rh)+ul*rl
        return chemical.numpy_add(product,np.zeros_like(u),error-e*r)


def torch_scale(u,e,r):
    product=u*r;cu=SPLITTER*u;cr=SPLITTER*r
    uh=cu-(cu-u);ul=u-uh;rh=cr-(cr-r);rl=r-rh
    error=((uh*rh-product)+uh*rl+ul*rh)+ul*rl
    return chemical.torch_add(product,torch.zeros_like(u),error-e*r)


def validate_correction(u,e):
    if (u.dtype!=np.float64 or e.dtype!=np.float64 or u.shape!=e.shape or
            not np.isfinite(u).all() or np.any(u<=0) or not np.isfinite(e).all() or
            np.any(abs(e)>2*np.spacing(u))):raise ValueError('Invalid chemical concentration/correction pair')


def retime(host,dt):
    old=host.config.dt
    host.step_number=_steps(host.time,dt)
    host.config.steps=_steps(host.config.steps*old,dt)
    host.config.save_every=_steps(host.config.save_every*old,dt)
    host.config.dt=dt;host.config.validate()


class MovingSimulation(PrecisionSimulation):
    def __init__(self,host,device='cuda:0'):
        super().__init__(host,'phase_carry',device)
        self.chemical_correction=torch.zeros((2,len(self.ids)),dtype=torch.float64,device=self.device)

    @classmethod
    def from_physical(cls,path,dt,initial,device='cuda:0'):
        """Explicit first intervention; historical physical residuals retained."""
        z=payload(path);host=AttributeSimulation.restore(path);retime(host,dt)
        if initial.dtype!=np.float64 or initial.shape!=(2,len(host.ids)):raise ValueError('Float64 initial chemistry required')
        validate_correction(initial,np.zeros_like(initial));host.activator,host.inhibitor=initial.copy()
        sim=cls(host,device);sim.phase_carry.copy_(torch.from_numpy(z['phase_carry']).to(sim.device));sim.rounding.copy_(torch.from_numpy(z['rounding']).to(sim.device))
        return sim

    def checkpoint(self,path):
        path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);tmp=path.with_name(path.stem+'.chemical.tmp.npz')
        super().checkpoint(tmp)
        with np.load(tmp) as z:out={k:z[k].copy() for k in z.files}
        u=np.stack([out['activator'],out['inhibitor']]);e=self.chemical_correction.cpu().numpy();validate_correction(u,e)
        out.update(chemical_correction=e,chemical_method=np.array(METHOD))
        np.savez_compressed(tmp,**out);tmp.replace(path)

    @classmethod
    def restore(cls,path,device='cuda:0'):
        # Read mandatory fields before allocating GPU memory. Missing correction
        # cannot masquerade as a historical first-start initialization.
        with np.load(path) as z:
            if 'chemical_correction' not in z or str(z.get('chemical_method',''))!=METHOD:
                raise ValueError('Missing/wrong mandatory chemical correction or method')
            e=z['chemical_correction'].copy();u=np.array([z['activator'],z['inhibitor']]);validate_correction(u,e)
        raw=payload(path);host=AttributeSimulation.restore(path);sim=cls(host,device)
        sim.phase_carry.copy_(torch.from_numpy(raw['phase_carry']).to(sim.device));sim.rounding.copy_(torch.from_numpy(raw['rounding']).to(sim.device))
        sim.chemical_correction.copy_(torch.from_numpy(e).to(sim.device));return sim

    @torch.inference_mode()
    def step(self):
        # Same beginning-geometry splitting and CUDA polarity/mechanics as the
        # accepted phase-carry implementation; no midpoint/order upgrade here.
        c=self.config;_,delta,normalized,_=self.matrices()
        step=chemical.TorchStepper(delta,c.dt,c.signal_beta,c.signal_da,c.signal_dh)
        u,e=step(torch.stack([self.activator,self.inhibitor]),self.chemical_correction)
        self.activator,self.inhibitor=u;self.chemical_correction=e
        aligned=(normalized@self.polarity).contiguous();p=torch.empty_like(self.polarity)
        call('gpu_polarity',[self.polarity,self.geometry[:,4:].contiguous(),self.activator,aligned,p],
            (len(self.ids),c.dt,c.polarity_rate,c.polarity_alignment,c.polarity_decay))
        self.polarity=p;response=torch.tanh(self.activator-1)
        tensions=c.surface_tension*(1+c.fate_tension*response)
        adhesion=c.adhesion*(1+c.fate_adhesion*response[:,None]*response[None,:]);adhesion.fill_diagonal_(0)
        attraction=(adhesion@self.shell2.flatten(1)).reshape(self.phi.shape)
        old=self.geometry[:,0].clone();centers=self.geometry[:,1:4].contiguous();vf=c.volume_stiffness*(self.target-old)/self.target
        amounts=torch.stack([old@self.activator,old@self.inhibitor])
        force([self.phi,self.excluded,attraction,centers,self.polarity,tensions,vf,self.gamma,self.updated,self.quality,self.phase_carry,self.rounding],
            (len(self.ids),c.grid,1,self.dx,c.extent,c.interface_width,c.polarity_tension,c.repulsion,c.dt))
        self.phi,self.updated=self.updated,self.phi;self.refresh();new=self.geometry[:,0]
        u,e=torch_scale(torch.stack([self.activator,self.inhibitor]),self.chemical_correction,old/new)
        if not bool(torch.isfinite(u).all() & (u>0).all() & torch.isfinite(e).all()):raise FloatingPointError('Invalid compensated dilution')
        self.activator,self.inhibitor=u;self.chemical_correction=e
        updated=torch.stack([new@self.activator,new@self.inhibitor]);self._dilution_error=(updated/amounts-1).abs().max()
        self.step_number+=1;self.time=self.step_number*c.dt

    def audit(self):
        result=super().audit()
        if not bool(torch.isfinite(self.chemical_correction).all()):raise FloatingPointError('Nonfinite chemical correction')
        return result

    def precision_diagnostics(self):
        result=super().precision_diagnostics();result['chemical_correction_abs_max']=float(self.chemical_correction.abs().max());result['chemical_method']=METHOD
        return result
