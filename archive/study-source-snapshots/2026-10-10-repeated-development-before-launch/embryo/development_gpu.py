"""Isolated resident development backend; qualification is required before science.

The mature backend and its restrictions remain intact. Event metadata/RNG live
on the host; spatial fields, ring mechanics, projection and partition stay on
the GPU. Visible fields define measured capacities, as in mature experiments.
"""
import ctypes
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess

import numpy as np
from scipy.ndimage import label
import torch
from . import compensated_chemistry as chemical
from .compensated_moving import MovingSimulation, torch_scale, validate_correction
from .gpu_backend import GpuSimulation, call
from .gpu_precision_control import force
from .attribute_development import AttributeSimulation

METHOD = 'compensated resident progressive development v1'
SCHEMA = 1
_LIBRARY = None


def library():
    global _LIBRARY
    if _LIBRARY is None:
        source=Path(__file__).with_name('gpu_development_kernels.cu')
        names=('gpu_precision_kernels.cu','gpu_kernels.cu','cuda_spatial_bench.cu','cuda_mechanics.cu')
        files=[source,*[source.with_name(n) for n in names]]
        flags=['-O3','-std=c++17','-arch=sm_61','--fmad=false','--ftz=false',
               '--prec-div=true','--prec-sqrt=true','-ccbin','/usr/bin/g++-12','-shared','-Xcompiler','-fPIC']
        key=hashlib.sha256(b''.join(f.read_bytes() for f in files)+repr(flags).encode()+subprocess.check_output(['nvcc','--version'])).hexdigest()[:20]
        folder=source.parent.parent/'output/kernel-cache/development';folder.mkdir(parents=True,exist_ok=True)
        binary=folder/f'development-{key}.so'
        if not binary.exists():
            temporary=binary.with_suffix(f'.{os.getpid()}.tmp.so')
            subprocess.run(['nvcc',*flags,str(source),'-lcublas','-o',str(temporary)],check=True)
            temporary.replace(binary)
        lib=ctypes.CDLL(str(binary));ptr,integer,real=ctypes.c_void_p,ctypes.c_int,ctypes.c_double
        lib.development_force.argtypes=[ptr]*13+[integer]*2+[real]*7+[ptr,integer]
        lib.development_force.restype=integer;lib.cuda_mechanics_error.restype=ctypes.c_char_p
        lib.binary_path=str(binary);_LIBRARY=lib
    return _LIBRARY


def ring_force(tensors,scalars):
    cells,n=scalars[:2];field=(cells,n,n,n)
    shapes=(field,field,field,(cells,3),(cells,3),(cells,),(cells,),field,field,(2,),field,(2,),(cells,7))
    kinds=(torch.float32,torch.float32,torch.float64,torch.float64,torch.float64,torch.float64,torch.float64,
           torch.float64,torch.float32,torch.int32,torch.float64,torch.int32,torch.float64)
    if len(tensors)!=13 or len(scalars)!=9:raise ValueError('Wrong developmental kernel argument count')
    device=tensors[0].device
    if any(not t.is_cuda or t.device!=device or not t.is_contiguous() or t.dtype!=k or tuple(t.shape)!=s
           for t,k,s in zip(tensors,kinds,shapes)):raise ValueError('Invalid developmental kernel buffer')
    lib=library();code=lib.development_force(*[t.data_ptr() for t in tensors],*scalars,
        torch.cuda.current_stream(device).cuda_stream,device.index)
    if code:raise RuntimeError(lib.cuda_mechanics_error().decode())


def occupancy64(p):return p*p*(3-2*p)


def store_phase(effective):
    """Retain the complete transformed phase pair, including pressure/splitting."""
    visible=effective.float();return visible,effective-visible.double()


def project_phase(field,carry,volume,dx):
    effective=field.double()+carry;weight=6*effective*(1-effective);shift=0.
    for _ in range(20):
        candidate=(effective+shift*weight).clamp(0,1)
        error=volume-float(occupancy64(candidate).sum()*dx**3)
        if abs(error)<=1e-7*volume:return (*store_phase(candidate),abs(shift))
        slope=float((6*candidate*(1-candidate)*weight).sum()*dx**3)
        if slope<=1e-12:break
        shift+=error/slope
    raise FloatingPointError('Developmental pressure projection did not converge')


def partition_chemistry(u,e,old_volume,new_volumes,draws):
    """Daughter concentration/correction pairs preserve mother species amounts.

    The signed correction is subtracted from the concentration. Anticorrelated
    partition fluctuations have zero volume-weighted mean for each species.
    """
    frac=new_volumes[0]/new_volumes.sum()
    base,err=torch_scale(u,e,old_volume/new_volumes.sum())
    draws=torch.as_tensor(draws,dtype=torch.float64,device=u.device).clamp(-.1,.1)
    factors=torch.stack((1+2*draws*(1-frac),1-2*draws*frac),dim=1)
    return torch_scale(base[:,None].expand(-1,2),err[:,None].expand(-1,2),factors)


class DevelopmentGpu(MovingSimulation):
    def __init__(self,host,device='cuda:0',*,_restore=False):
        # Do not weaken GpuSimulation's accepted mature-state constructor.
        if host.time!=0 and not _restore:raise ValueError('Evolved state requires mandatory DevelopmentGpu.restore corrections')
        c=host.config;c.validate();count=len(host.ids)
        if (host.attribute_mode not in ('direct','no_feedback') or not c.signaling or
                c.signal_transport!='conservative' or c.differentiation or np.any(host.fate!=0) or
                not 1<=count<=c.max_cells or count+len(host.divisions)>c.max_cells):
            raise ValueError('Development requires unfated conservative attribute cells')
        if host.phi.dtype!=np.float32 or host.phi.shape!=(count,c.grid,c.grid,c.grid):raise ValueError('Wrong phase dtype/shape')
        if any(not np.isfinite(x).all() for x in (host.phi,host.activator,host.inhibitor,host.polarity,host.target)):
            raise ValueError('Nonfinite developmental start')
        if np.any(host.phi<0) or np.any(host.phi>1) or np.any(host.activator<=0) or np.any(host.inhibitor<=0) or np.any(host.target<=0):
            raise ValueError('Invalid developmental start')
        self.device=torch.device(device)
        if self.device.type!='cuda':raise ValueError('Resident development requires CUDA')
        if self.device.index is None:self.device=torch.device('cuda',torch.cuda.current_device())
        torch.backends.cuda.matmul.allow_tf32=False
        self.host,self.config=host,c;self.ids=host.ids.copy()
        self.time,self.step_number,self.dx=host.time,host.step_number,host.dx
        upload=lambda x:torch.from_numpy(np.ascontiguousarray(x)).to(self.device)
        self.phi,self.activator,self.inhibitor,self.polarity,self.target=[upload(x) for x in (host.phi,host.activator,host.inhibitor,host.polarity,host.target)]
        if any(t.dtype!=torch.float64 for t in (self.activator,self.inhibitor,self.polarity,self.target)):raise ValueError('Float64 attributes/targets required')
        self.phase_carry=torch.zeros_like(self.phi,dtype=torch.float64)
        self.chemical_correction=torch.zeros((2,count),dtype=torch.float64,device=self.device)
        self.rounding=torch.zeros(2,dtype=torch.int32,device=self.device);self.precision_arm='phase_carry'
        self.clipped_fraction=self.volume_projection_max=self.dilution_amount_error=0.
        self.event_audits=[];self._allocate()
        axis=(torch.arange(c.grid,dtype=torch.float64,device=self.device)+.5)*self.dx-c.extent
        self.xyz=torch.stack(torch.meshgrid(axis,axis,axis,indexing='ij'));self.refresh()

    def _allocate(self):
        n=self.config.grid;count=len(self.ids);chunks=(n**3+255)//256
        self.h=torch.empty_like(self.phi);self.shell=torch.empty_like(self.phi)
        self.shell2=torch.empty_like(self.phi,dtype=torch.float64);self.occupied=torch.empty_like(self.phi[0])
        self.excluded=torch.empty_like(self.phi);self.gamma=torch.empty_like(self.shell2);self.updated=torch.empty_like(self.phi)
        self.partial=torch.empty((count,8,chunks),dtype=torch.float64,device=self.device)
        self.geometry=torch.empty((count,7),dtype=torch.float64,device=self.device)
        self.union_partial=torch.empty((10,chunks),dtype=torch.float64,device=self.device)
        self.quality=torch.zeros(2,dtype=torch.int32,device=self.device)

    @torch.inference_mode()
    def divide(self,index,direction=None):
        parent=int(self.ids[index]);host=self.host;c=self.config
        if parent in host.divisions or len(self.ids)+len(host.divisions)>=c.max_cells:return False
        relative=self.xyz-self.geometry[index,1:4,None,None,None]
        if direction is None:
            sample=host.rng.normal(size=3)
            if c.division_orientation=='shape':
                offset=relative.flatten(1);h=self.h[index].double().flatten()
                covariance=((offset*h)@offset.T/h.sum()).cpu().numpy();val,vec=np.linalg.eigh(covariance)
                basis=vec[:,val>=val[-1]*(1-c.axis_degeneracy)];sample=basis@(basis.T@sample)
            direction=sample
        direction=np.asarray(direction,dtype=float)
        if direction.shape!=(3,) or not np.isfinite(direction).all() or np.linalg.norm(direction)==0:raise ValueError('Finite nonzero division axis required')
        direction=direction/np.linalg.norm(direction);axis=torch.as_tensor(direction,device=self.device)
        axial=torch.einsum('d,dijk->ijk',axis,relative);mother=occupancy64(self.phi[index].double()+self.phase_carry[index])
        lo,hi=float(axial.min()),float(axial.max());total=mother.sum()
        for _ in range(35):
            offset=(lo+hi)/2;weight=.5*(1+torch.tanh((axial-offset)/c.interface_width))
            if bool((mother*weight).sum()>.5*total):lo=offset
            else:hi=offset
        radial=(relative.square().sum(0)-axial.square()).clamp_min(0).sqrt()
        inside=self.phi[index]>=.5
        if not bool(inside.any()):raise FloatingPointError('No resolved mother interior')
        host.divisions[parent]=dict(start=self.time,axis=direction.tolist(),offset=offset,radius=float(radial[inside].max()),volume=float(self.geometry[index,0]))
        host.due[index]=np.inf;host.lineage[parent].update(division_start=self.time,division_axis=direction.tolist());return True

    def _cleavage(self,index,event):
        relative=self.xyz-self.geometry[index,1:4,None,None,None]
        axis=torch.tensor(event['axis'],dtype=torch.float64,device=self.device)
        axial=torch.einsum('d,dijk->ijk',axis,relative)-event['offset']
        effective=self.phi[index].double()+self.phase_carry[index]
        mother=occupancy64(effective);weight=.5*(1+torch.tanh(axial/self.config.interface_width))
        first=mother*weight;second=mother-first
        hs=torch.stack((first,second));daughters=.5-torch.sin(torch.asin((1-2*hs).clamp(-1,1))/3)
        visible,carry=store_phase(daughters)
        band=axial.abs()<=max(self.config.interface_width*.5,self.dx*.5)
        neck=float(self.phi[index][band].max());overlap=float((visible[0].double().square()*visible[1].double().square()).sum()*self.dx**3/event['volume'])
        fraction=float(first.sum()/mother.sum())
        jump=float((occupancy64(daughters).sum(0)-mother).abs().max())
        return visible,carry,neck,overlap,fraction,jump

    @torch.inference_mode()
    def _finish_divisions(self):
        c=self.config
        for parent,event in list(self.host.divisions.items()):
            if self.time-event['start']<c.cytokinesis_duration:continue
            index=int(np.flatnonzero(self.ids==parent)[0]);daughters,carry,neck,overlap,fraction,jump=self._cleavage(index,event)
            if neck>c.neck_threshold or overlap>c.division_overlap_tolerance or not .1<fraction<.9:continue
            self._complete_division(index,daughters,carry,neck,overlap,fraction,jump)

    def _complete_division(self,index,daughters,carry,neck,overlap,fraction,jump):
        host=self.host;c=self.config;parent=int(self.ids[index]);before_ids=self.ids.tolist()
        old=self.geometry[:,0].clone();u=torch.stack((self.activator,self.inhibitor));e=self.chemical_correction
        before=(u*old).sum(1);target_before=float(self.target.sum());volume_before=float(old.sum())
        keep=np.arange(len(self.ids))!=index;mask=torch.as_tensor(keep,device=self.device)
        # Carry is transformed through occupancy partition, never zeroed.
        self.phi=torch.cat((self.phi[mask],daughters));self.phase_carry=torch.cat((self.phase_carry[mask],carry))
        self.ids=np.r_[self.ids[keep],host.next_id,host.next_id+1];children=self.ids[-2:].tolist();host.next_id+=2
        self._allocate();self.refresh();new=self.geometry[-2:,0];actual_fraction=new[0]/new.sum()
        draw=[host.signal_rng.normal(scale=c.signal_partition_noise) for _ in range(2)]
        daughter_u,daughter_e=partition_chemistry(u[:,index],e[:,index],old[index],new,draw)
        combined=torch.cat((u[:,mask],daughter_u),dim=1);self.activator,self.inhibitor=combined
        self.chemical_correction=torch.cat((e[:,mask],daughter_e),dim=1)
        first_target=self.target[index]*actual_fraction;second_target=self.target[index]-first_target
        self.target=torch.cat((self.target[mask],torch.stack((first_target,second_target))))
        self.polarity=torch.cat((self.polarity[mask],self.polarity[index].expand(2,-1))).contiguous()
        host.parents=np.r_[host.parents[keep],parent,parent]
        host.due=np.r_[host.due[keep],self.time+host._cycle(),self.time+host._cycle()]
        host.lineage[parent].update(division=self.time,neck_at_abscission=neck,overlap_at_abscission=overlap)
        host.lineage.extend(dict(id=child,parent=parent,birth=self.time,division=None) for child in children)
        del host.divisions[parent];host.ids=self.ids.copy();host.fate=np.zeros(len(self.ids))
        after=(combined*self.geometry[:,0]).sum(1)
        audit=dict(time=self.time,parent=parent,children=children,ids_before=before_ids,ids_after=self.ids.tolist(),
            neck=neck,overlap=overlap,fraction=fraction,actual_fraction=float(actual_fraction),draws=draw,
            volumes_before=old.cpu().numpy().tolist(),volumes_after=self.geometry[:,0].cpu().numpy().tolist(),
            chemistry_before=u.cpu().numpy().tolist(),chemistry_after=combined.cpu().numpy().tolist(),
            chemical_correction_before=e.cpu().numpy().tolist(),chemical_correction_after=self.chemical_correction.cpu().numpy().tolist(),
            amount_relative_jump=float((after/before-1).abs().max()),target_relative_jump=abs(float(self.target.sum())/target_before-1),
            volume_relative_jump=abs(float(self.geometry[:,0].sum())/volume_before-1),effective_occupancy_absolute_jump=jump,
            daughter_components=[int(label(f>=.5)[1]) for f in daughters.cpu().numpy()])
        self.event_audits.append(audit);host.graph_events.append(audit.copy())

    @torch.inference_mode()
    def step(self):
        c=self.config;_,delta,normalized,_=self.matrices()
        u,e=chemical.TorchStepper(delta,c.dt,c.signal_beta,c.signal_da,c.signal_dh)(torch.stack((self.activator,self.inhibitor)),self.chemical_correction)
        self.activator,self.inhibitor=u;self.chemical_correction=e
        if c.polarity_enabled:
            p=torch.empty_like(self.polarity);aligned=(normalized@self.polarity).contiguous()
            call('gpu_polarity',[self.polarity,self.geometry[:,4:].contiguous(),self.activator,aligned,p],
                 (len(self.ids),c.dt,c.polarity_rate,c.polarity_alignment,c.polarity_decay));self.polarity=p
        response=torch.tanh(self.activator-1)
        tensions=torch.full_like(self.activator,c.surface_tension);adhesion=torch.full_like(delta,c.adhesion)
        if c.feedback:
            tensions*=1+c.fate_tension*response;adhesion*=1+c.fate_adhesion*response[:,None]*response[None,:]
        adhesion.fill_diagonal_(0);attraction=(adhesion@self.shell2.flatten(1)).reshape(self.phi.shape)
        old=self.geometry[:,0].clone();centers=self.geometry[:,1:4].contiguous();vf=c.volume_stiffness*(self.target-old)/self.target
        before=(u*old).sum(1);contrast=c.polarity_tension if c.feedback and c.polarity_enabled else 0.
        tensors=[self.phi,self.excluded,attraction,centers,self.polarity,tensions,vf,self.gamma,self.updated,self.quality,self.phase_carry,self.rounding]
        if self.host.divisions:
            events=np.zeros((len(self.ids),7))
            for i,cell in enumerate(self.ids):
                ev=self.host.divisions.get(int(cell))
                if ev is not None:
                    progress=np.clip((self.time-ev['start'])/c.cytokinesis_duration,0,1);ramp=progress**2*(3-2*progress)
                    events[i]=[1,*ev['axis'],ev['offset'],ev['radius'],ramp]
            packed=torch.as_tensor(events,device=self.device)
            ring_force(tensors+[packed],(len(self.ids),c.grid,self.dx,c.extent,c.interface_width,contrast,c.repulsion,c.dt,c.ring_strength))
        else:
            force(tensors,(len(self.ids),c.grid,1,self.dx,c.extent,c.interface_width,contrast,c.repulsion,c.dt))
        self.phi,self.updated=self.updated,self.phi
        self.clipped_fraction=float(self.quality[0])/self.phi.numel();self.volume_projection_max=0.
        if bool(self.quality[1]):raise FloatingPointError('Nonfinite developmental force')
        if self.clipped_fraction:raise FloatingPointError('Developmental force clipped; preserve failure and refine')
        for i,cell in enumerate(self.ids):
            ev=self.host.divisions.get(int(cell))
            if ev is not None:
                field,carry,shift=project_phase(self.phi[i],self.phase_carry[i],ev['volume'],self.dx)
                self.phi[i].copy_(field);self.phase_carry[i].copy_(carry);self.volume_projection_max=max(self.volume_projection_max,shift)
        self.refresh();u,e=torch_scale(u,self.chemical_correction,old/self.geometry[:,0])
        self.activator,self.inhibitor=u;self.chemical_correction=e
        self._dilution_error=(((u*self.geometry[:,0]).sum(1)/before)-1).abs().max()
        self.step_number+=1;self.time=self.step_number*c.dt;self._finish_divisions()
        self.host.time,self.host.step_number=self.time,self.step_number
        while len(self.ids)+len(self.host.divisions)<c.max_cells and np.any(self.host.due<=self.time):self.divide(int(np.argmin(self.host.due)))

    def audit(self):
        result=GpuSimulation.audit(self)
        if not bool(torch.isfinite(self.phase_carry).all() & torch.isfinite(self.chemical_correction).all()):raise FloatingPointError('Invalid development corrections')
        result['max_clipping']=self.clipped_fraction;result['volume_projection_max']=self.volume_projection_max
        return result

    def to_cpu(self):
        host=super().to_cpu();host.ids=self.ids.copy();host.target=self.target.cpu().numpy().copy()
        host.fate=np.zeros(len(self.ids));host.volume_projection_max=self.volume_projection_max;return host

    def checkpoint(self,path):
        path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);tmp=path.with_name(path.stem+'.development.tmp.npz')
        self.to_cpu().checkpoint(tmp)
        with np.load(tmp,allow_pickle=False) as z:out={k:z[k].copy() for k in z.files}
        correction=self.chemical_correction.cpu().numpy();validate_correction(np.stack((out['activator'],out['inhibitor'])),correction)
        out.update(phase_carry=self.phase_carry.cpu().numpy(),rounding=self.rounding.cpu().numpy(),precision_arm=np.array('phase_carry'),
                   chemical_correction=correction,chemical_method=np.array(METHOD),development_schema=np.array(SCHEMA),event_audits=np.array(json.dumps(self.event_audits)))
        np.savez_compressed(tmp,**out);tmp.replace(path)

    @classmethod
    def restore(cls,path,device='cuda:0'):
        with np.load(path,allow_pickle=False) as z:
            if str(z.get('chemical_method',''))!=METHOD or int(z.get('development_schema',-1))!=SCHEMA:raise ValueError('Mandatory developmental correction/method/schema required')
            carry=z['phase_carry'].copy();counts=z['rounding'].copy();e=z['chemical_correction'].copy()
            validate_correction(np.stack((z['activator'],z['inhibitor'])),e)
            if carry.dtype!=np.float64 or carry.shape!=z['phi'].shape or not np.isfinite(carry).all() or counts.dtype!=np.int32 or counts.shape!=(2,) or np.any(counts<0):raise ValueError('Invalid developmental carry/counters')
            events=json.loads(str(z['event_audits']))
        host=AttributeSimulation.restore(path);sim=cls(host,device,_restore=True)
        sim.phase_carry.copy_(torch.as_tensor(carry,device=sim.device));sim.rounding.copy_(torch.as_tensor(counts,device=sim.device))
        sim.chemical_correction.copy_(torch.as_tensor(e,device=sim.device));sim.event_audits=events;return sim
