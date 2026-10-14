"""Independent NumPy developmental arithmetic for GPU event qualification.

Uses the existing CPU geometry/transport/polarity definitions. Ring and phase
updates are computed independently; zero polarity uses the same conservative
face-flux formula as nonzero polarity. No GPU or development-kernel calls.
"""
import numpy as np
from .attribute_development import AttributeSimulation
from .model import occupancy
from .polarity import exposure_cue,evolve,tension_field,flux_divergence
from .signaling import normalized_graph
from .compensated_chemistry import NumpyStepper
from .compensated_moving import numpy_scale


def project_effective(field,carry,volume,dx):
    p=field.astype(np.float64)+carry;weight=6*p*(1-p);shift=0.
    for _ in range(20):
        candidate=np.clip(p+shift*weight,0,1)
        error=volume-float(occupancy(candidate).sum()*dx**3)
        if abs(error)<=1e-7*volume:
            visible=candidate.astype(np.float32);return visible,candidate-visible.astype(np.float64),abs(shift)
        slope=float((6*candidate*(1-candidate)*weight).sum()*dx**3)
        if slope<=1e-12:break
        shift+=error/slope
    raise FloatingPointError('Reference pressure projection did not converge')


def partition_numpy(u,e,old_volume,new_volumes,draws):
    f=new_volumes[0]/new_volumes.sum();base,err=numpy_scale(u,e,old_volume/new_volumes.sum())
    draw=np.clip(np.asarray(draws),-.1,.1)
    factors=np.stack((1+2*draw*(1-f),1-2*draw*f),axis=1)
    return numpy_scale(base[:,None],err[:,None],factors)


class ReferenceDevelopment(AttributeSimulation):
    def __init__(self,config=None):
        super().__init__(config);self.phase_carry=np.zeros_like(self.phi,dtype=np.float64)
        self.chemical_correction=np.zeros((2,len(self.ids)));self.event_audits=[];self.dilution_amount_error=0.

    def divide(self,index,direction=None):
        c=self.config;parent=int(self.ids[index])
        if parent in self.divisions or len(self.ids)+len(self.divisions)>=c.max_cells:return False
        direction=self.division_direction(index) if direction is None else np.asarray(direction,dtype=float)
        if direction.shape!=(3,) or not np.isfinite(direction).all() or np.linalg.norm(direction)==0:raise ValueError('Invalid division axis')
        direction=direction/np.linalg.norm(direction);relative=self.xyz-self.centers()[index,:,None,None,None]
        axial=np.einsum('d,dijk->ijk',direction,relative);mother=occupancy(self.phi[index].astype(float)+self.phase_carry[index])
        lo,hi=float(axial.min()),float(axial.max());total=mother.sum()
        for _ in range(35):
            offset=(lo+hi)/2;weight=.5*(1+np.tanh((axial-offset)/c.interface_width))
            if (mother*weight).sum()>.5*total:lo=offset
            else:hi=offset
        radial=np.sqrt(np.maximum(np.sum(relative**2,axis=0)-axial**2,0))
        self.divisions[parent]=dict(start=self.time,axis=direction.tolist(),offset=offset,radius=float(radial[self.phi[index]>=.5].max()),volume=float(self.volumes()[index]))
        self.due[index]=np.inf;self.lineage[parent].update(division_start=self.time,division_axis=direction.tolist());return True

    def mechanical_step(self):
        c=self.config;phi=self.phi;h=occupancy(phi);volumes=self.volumes();centers=self.centers()
        derivative=2*phi*(1-phi)*(1-2*phi);shell2=(phi*(1-phi))**2
        tension,adhesion=self.material_coefficients();attraction=(adhesion@shell2.reshape(len(phi),-1)).reshape(phi.shape)
        excluded=np.sum(phi**2,axis=0)[None]-phi**2;vf=c.volume_stiffness*(self.target-volumes)/self.target
        updated=np.empty_like(phi);carry=np.empty_like(self.phase_carry);clipped=0
        self.volume_projection_max=0.
        for i in range(len(phi)):
            if c.feedback and c.polarity_enabled:
                gamma=tension_field(self.xyz-centers[i,:,None,None,None],self.polarity[i],tension[i],c.polarity_tension,c.interface_width)
            else:gamma=np.full_like(phi[i],tension[i],dtype=float)
            f=c.interface_width**2*flux_divergence(phi[i],gamma,self.dx)-gamma*derivative[i]
            f+=vf[i]*6*phi[i]*(1-phi[i]);f-=c.repulsion*phi[i]*excluded[i];f+=derivative[i]*attraction[i]
            ev=self.divisions.get(int(self.ids[i]))
            if ev is not None:
                relative=self.xyz-centers[i,:,None,None,None];axial=np.einsum('d,dijk->ijk',ev['axis'],relative)
                radial=np.sqrt(np.maximum(np.sum(relative**2,axis=0)-axial**2,0))
                progress=np.clip((self.time-ev['start'])/c.cytokinesis_duration,0,1);ramp=progress**2*(3-2*progress)
                band=np.exp(-.5*((axial-ev['offset'])/c.interface_width)**2)
                outside=.5*(1+np.tanh((radial-ev['radius']*(1-ramp))/c.interface_width))
                f-=c.ring_strength*ramp*band*outside*6*phi[i]*(1-phi[i])
            effective=phi[i].astype(float)+(self.phase_carry[i]+c.dt*f)
            rounded=effective.astype(np.float32);outside=(rounded<0)|(rounded>1);clipped+=int(outside.sum())
            carry[i]=np.where(outside,0.,effective-rounded.astype(float));updated[i]=np.clip(rounded,0,1)
            if ev is not None:
                updated[i],carry[i],shift=project_effective(updated[i],carry[i],ev['volume'],self.dx)
                self.volume_projection_max=max(self.volume_projection_max,shift)
        self.phi,self.phase_carry=updated,carry;self.clipped_fraction=clipped/phi.size

    def _finish_divisions(self):
        c=self.config
        for parent,ev in list(self.divisions.items()):
            if self.time-ev['start']<c.cytokinesis_duration:continue
            i=int(np.flatnonzero(self.ids==parent)[0]);axial,_=self._division_coordinates(i,ev)
            mother=occupancy(self.phi[i].astype(float)+self.phase_carry[i]);w=.5*(1+np.tanh(axial/c.interface_width))
            hs=np.stack((mother*w,mother-mother*w));effective=.5-np.sin(np.arcsin(np.clip(1-2*hs,-1,1))/3)
            fields=effective.astype(np.float32);carry=effective-fields.astype(float)
            band=abs(axial)<=max(c.interface_width*.5,self.dx*.5);neck=float(self.phi[i][band].max())
            overlap=float((fields[0].astype(float)**2*fields[1].astype(float)**2).sum()*self.dx**3/ev['volume'])
            fraction=float(hs[0].sum()/mother.sum())
            if neck>c.neck_threshold or overlap>c.division_overlap_tolerance or not .1<fraction<.9:continue
            self._complete(i,fields,carry,neck,overlap,fraction,float(abs(occupancy(effective).sum(0)-mother).max()))

    def _complete(self,i,fields,carry,neck,overlap,fraction,jump):
        c=self.config;parent=int(self.ids[i]);old=self.volumes();u=np.stack((self.activator,self.inhibitor));e=self.chemical_correction
        before=(u*old).sum(1);volume_before=old.sum();target_before=self.target.sum();mask=np.arange(len(self.ids))!=i
        self.phi=np.concatenate((self.phi[mask],fields));self.phase_carry=np.concatenate((self.phase_carry[mask],carry))
        new=self.volumes();draw=[self.signal_rng.normal(scale=c.signal_partition_noise) for _ in range(2)]
        a,b=partition_numpy(u[:,i],e[:,i],old[i],new[-2:],draw);u=np.concatenate((u[:,mask],a),axis=1)
        self.chemical_correction=np.concatenate((e[:,mask],b),axis=1);self.activator,self.inhibitor=u
        f=new[-2]/new[-2:].sum();first=self.target[i]*f;self.target=np.r_[self.target[mask],first,self.target[i]-first]
        self.polarity=np.concatenate((self.polarity[mask],np.repeat(self.polarity[i][None],2,axis=0)))
        children=[self.next_id,self.next_id+1];self.next_id+=2;self.ids=np.r_[self.ids[mask],children]
        self.parents=np.r_[self.parents[mask],parent,parent];self.due=np.r_[self.due[mask],self.time+self._cycle(),self.time+self._cycle()]
        self.lineage[parent].update(division=self.time,neck_at_abscission=neck,overlap_at_abscission=overlap)
        self.lineage.extend(dict(id=child,parent=parent,birth=self.time,division=None) for child in children);del self.divisions[parent]
        self.fate=np.zeros(len(self.ids));after=(u*new).sum(1)
        self.event_audits.append(dict(time=self.time,parent=parent,children=children,neck=neck,overlap=overlap,fraction=fraction,draws=draw,
            amount_relative_jump=float(abs(after/before-1).max()),target_relative_jump=abs(self.target.sum()/target_before-1),
            volume_relative_jump=abs(new.sum()/volume_before-1),effective_occupancy_absolute_jump=jump))

    def step(self):
        c=self.config;contact,_=self.contacts();delta=self.signaling_graph(contact).delta
        u,e=NumpyStepper(delta,c.dt,c.signal_beta,c.signal_da,c.signal_dh)(np.stack((self.activator,self.inhibitor)),self.chemical_correction)
        self.activator,self.inhibitor=u;self.chemical_correction=e
        if c.polarity_enabled:
            cutoff=c.graph_contact_cutoff if c.polarity_contact_cutoff==-1 else c.polarity_contact_cutoff
            self.polarity=evolve(self.polarity,exposure_cue(self.phi,self.dx),self.activator,normalized_graph(contact,cutoff).delta,
                c.dt,c.polarity_rate,c.polarity_alignment,c.polarity_decay)
        old=self.volumes();before=(u*old).sum(1);self.mechanical_step();new=self.volumes()
        u,e=numpy_scale(u,self.chemical_correction,old/new);self.activator,self.inhibitor=u;self.chemical_correction=e
        self.dilution_amount_error=float(abs((u*new).sum(1)/before-1).max())
        self.step_number+=1;self.time=self.step_number*c.dt;self._finish_divisions()
        while len(self.ids)+len(self.divisions)<c.max_cells and np.any(self.due<=self.time):self.divide(int(np.argmin(self.due)))
