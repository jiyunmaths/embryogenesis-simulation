"""Fresh development without a fate switch or prescribed identity classes."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.ndimage import laplace

from .model import Config, Simulation, occupancy
from .polarity import tension_field, flux_divergence
from .resolution import write_json, _steps

MODES=('direct','no_feedback')
FEATURES=('log_activator','log_inhibitor','polarity_magnitude','axis_ratio_minus_one','asphericity')


class AttributeSimulation(Simulation):
    """Compatibility fate arrays stay zero and never enter dynamics or exports."""
    def __init__(self,config=None,mode=None):
        if mode is None:mode='direct' if config is None or config.feedback else 'no_feedback'
        if mode not in MODES:raise ValueError('unknown attribute-development mode')
        super().__init__(config)
        self.attribute_mode=mode
        self.config.differentiation=False
        self.config.fate_rate=0.;self.config.fate_noise=0.;self.config.partition_noise=0.
        self.config.feedback=mode=='direct'
        self.fate[:]=0

    def update_fate(self,contacts,exposure):
        self.fate[:]=0

    def material_coefficients(self):
        c=self.config
        response=np.tanh(self.activator-1.)
        tension=np.full(len(self.phi),c.surface_tension,dtype=float)
        attraction=np.full((len(self.phi),len(self.phi)),c.adhesion,dtype=float)
        if c.feedback:
            tension*=1+c.fate_tension*response
            attraction*=1+c.fate_adhesion*response[:,None]*response[None,:]
        np.fill_diagonal(attraction,0)
        return tension,attraction

    def mechanical_step(self):
        # Same spatial mechanics as the core model; material coefficients now
        # depend instantaneously on measured activator, never on stored fate.
        c=self.config;phi=self.phi
        shell2=(phi*(1-phi))**2
        derivative=2*phi*(1-phi)*(1-2*phi)
        tension,adhesion=self.material_coefficients()
        attract=(adhesion@shell2.reshape(len(phi),-1)).reshape(phi.shape)
        exclude=np.sum(phi**2,axis=0)[None]-phi**2
        h=occupancy(phi);volumes=h.sum(axis=(1,2,3),dtype=np.float64)*self.dx**3
        volume_force=c.volume_stiffness*(self.target-volumes)/self.target
        updated=np.empty_like(phi);self.volume_projection_max=0.
        centers=np.einsum('nijk,dijk->nd',h,self.xyz)*self.dx**3/volumes[:,None]
        for i in range(len(phi)):
            if c.feedback and c.polarity_enabled and np.linalg.norm(self.polarity[i])>1e-10:
                gamma=tension_field(self.xyz-centers[i,:,None,None,None],self.polarity[i],tension[i],c.polarity_tension,c.interface_width)
                force=c.interface_width**2*flux_divergence(phi[i],gamma,self.dx)-gamma*derivative[i]
            else:
                diffusion=laplace(phi[i],mode='nearest')/self.dx**2
                force=tension[i]*(c.interface_width**2*diffusion-derivative[i])
            force+=volume_force[i]*6*phi[i]*(1-phi[i])
            force-=c.repulsion*phi[i]*exclude[i]
            force+=derivative[i]*attract[i]
            event=self.divisions.get(int(self.ids[i]))
            if event is not None:
                axial,radial=self._division_coordinates(i,event)
                progress=np.clip((self.time-event['start'])/c.cytokinesis_duration,0,1)
                ramp=progress**2*(3-2*progress);ring_radius=event['radius']*(1-ramp)
                band=np.exp(-.5*(axial/c.interface_width)**2)
                outside=.5*(1+np.tanh((radial-ring_radius)/c.interface_width))
                force-=c.ring_strength*ramp*band*outside*6*phi[i]*(1-phi[i])
            updated[i]=phi[i]+c.dt*force
        self.clipped_fraction=float(np.mean((updated<0)|(updated>1)))
        self.phi=np.clip(updated,0,1)
        for i,cell_id in enumerate(self.ids):
            event=self.divisions.get(int(cell_id))
            if event is not None:self.phi[i]=self._project_volume(self.phi[i],event['volume'])

    def metrics(self):
        result=super().metrics()
        for key in ('fate_a','fate_b','uncommitted','fate_separation'):result.pop(key,None)
        result['attribute_mode']=self.attribute_mode
        return result

    def surfaces(self,max_points=200,*,include_mesh=True):
        cells=super().surfaces(max_points,include_mesh=include_mesh)
        for cell in cells:cell.pop('fate',None)
        return cells

    def checkpoint(self,path):
        if np.any(self.fate!=0):raise ValueError('unused compatibility fate must remain zero')
        super().checkpoint(path)
        with np.load(path,allow_pickle=False) as f:payload={key:f[key].copy() for key in f.files}
        payload['attribute_mode']=np.array(self.attribute_mode)
        payload['attribute_schema']=np.array(1)
        np.savez_compressed(path,**payload)

    @classmethod
    def restore(cls,path):
        with np.load(path,allow_pickle=False) as f:
            if 'attribute_schema' not in f or int(f['attribute_schema'])!=1:raise ValueError('requires an attribute-development checkpoint')
            mode=str(f['attribute_mode'])
        sim=super().restore(path);sim.attribute_mode=mode
        if np.any(sim.fate!=0) or sim.config.differentiation or sim.config.feedback!=(mode=='direct'):
            raise ValueError('inconsistent attribute checkpoint')
        return sim


def attributes(sim):
    h=occupancy(sim.phi).reshape(len(sim.phi),-1).astype(float)
    xyz=sim.xyz.reshape(3,-1);weight=h.sum(axis=1)
    centers=h@xyz.T/weight[:,None]
    eigen=[]
    for i in range(len(h)):
        offset=xyz-centers[i,:,None]
        covariance=(offset*h[i])@offset.T/weight[i]
        eigen.append(np.linalg.eigvalsh(covariance))
    eigen=np.array(eigen)
    ratio=np.sqrt(eigen[:,-1]/np.maximum(eigen[:,0],1e-30))
    asphericity=1.5*np.sum((eigen-eigen.mean(axis=1)[:,None])**2,axis=1)/np.sum(eigen,axis=1)**2
    contacts,exposure=sim.contacts();base_tension,attraction=sim.material_coefficients()
    vectors=np.column_stack([np.log(sim.activator),np.log(sim.inhibitor),np.linalg.norm(sim.polarity,axis=1),ratio-1,asphericity])
    cells=[]
    for i in range(len(h)):
        cells.append({'id':int(sim.ids[i]),'parent':int(sim.parents[i]),'attributes':dict(zip(FEATURES,vectors[i].tolist())),
                      'context':{'center':centers[i].tolist(),'exposure':float(exposure[i]),'volume':float(weight[i]*sim.dx**3),
                                 'target_volume':float(sim.target[i]),'raw_contact_sum':float(contacts[i].sum()),'dividing':int(sim.ids[i]) in sim.divisions},
                      'derived_material':{'baseline_tension':float(base_tension[i]),'attraction_to_cells':attraction[i].tolist()}})
    return {'metrics':sim.metrics(),'cells':cells,'cell_order':sim.ids.tolist(),
            'attribute_mean':vectors.mean(axis=0).tolist(),'attribute_std':vectors.std(axis=0).tolist()}


def prepare(output,until=90.,grid=72,dt=.0075,interval=.6):
    output=Path(output)
    if output.exists():raise FileExistsError('choose a fresh output directory')
    config=Config(grid=grid,extent=2.24,dt=dt,steps=_steps(until,dt),save_every=_steps(interval,dt),differentiation=False,fate_rate=0.,fate_noise=0.,partition_noise=0.)
    config.validate()
    if config.steps%config.save_every:raise ValueError('observations must align')
    p={'config':asdict(config),'modes':list(MODES),'until':until,'interval':interval,'late_duration':15.,'features':list(FEATURES),
       'coupling':{'response':'tanh(activator-1)','tension_contrast':config.fate_tension,'attraction_contrast':config.fate_adhesion,
                   'law':'gamma0*(1+c_gamma*response_i); A0*(1+c_A*response_i*response_j). Existing polarity tension retained only in direct mode.'},
       'criteria':{'volume_max':.05,'boundary_max':.01,'radius_min':4.,'max_clipping':0.},
       'source_sha256':{str(Path(__file__).with_name(f).resolve()):hashlib.sha256(Path(__file__).with_name(f).read_bytes()).hexdigest() for f in ['attribute_development.py','model.py','signaling.py','transport.py','polarity.py']},
       'scope':'Fresh analytic zygote, seed 7, 16-cell cap. No downstream fate drift, fate noise, fate labels, or prescribed cluster count. Direct continuous activity-dependent mechanics versus no regulatory mechanical feedback. Same initial geometry and random streams; later division histories may differ. Single-seed exploratory phenotype assay, not verified emergent identities. Historical config names fate_tension/fate_adhesion store material contrast numbers only.'}
    output.mkdir(parents=True);write_json(output/'protocol.json',p);write_json(output/'status.json',{'state':'prepared','protocol_sha256':hashlib.sha256((output/'protocol.json').read_bytes()).hexdigest()});return p


def run_arm(output,mode):
    output=Path(output);p=json.loads((output/'protocol.json').read_text());path=output/mode;path.mkdir()
    sim=AttributeSimulation(Config(**p['config']),mode);history=[];started=time.monotonic()
    max_volume=0.;min_radius=float('inf');max_clip=0.;max_boundary=0.
    for step in range(sim.config.steps+1):
        volume=sim.volumes();max_volume=max(max_volume,float(np.max(abs(volume/sim.target-1))))
        min_radius=min(min_radius,float(np.min((3*volume/(4*np.pi))**(1/3))/sim.dx));max_clip=max(max_clip,sim.clipped_fraction)
        if step%sim.config.save_every==0:
            row=attributes(sim);history.append(row);max_boundary=max(max_boundary,row['metrics']['boundary_occupancy'])
            write_json(path/'history.json',history);write_json(path/'status.json',{'state':'running','time':sim.time})
            if step%(10*sim.config.save_every)==0 or step==sim.config.steps:
                sim.checkpoint(path/f'state-{sim.time:g}.npz')
                write_json(path/f'surface-{sim.time:g}.json',{'metrics':row['metrics'],'cells':sim.surfaces()})
                print(f'{mode}: t={sim.time:g}, elapsed={time.monotonic()-started:.1f}s',flush=True)
        if step<sim.config.steps:sim.step()
    sim.checkpoint(path/'final_state.npz');write_json(path/'lineage.json',sim.lineage)
    result={'max_volume_error':max_volume,'min_radius':min_radius,'max_clipping':max_clip,'max_sampled_boundary':max_boundary,
            'final_cells':len(sim.phi),'active_divisions':len(sim.divisions),'elapsed_seconds':time.monotonic()-started}
    write_json(path/'analysis.json',result);write_json(path/'status.json',{'state':'completed'});return result


def compare(output):
    output=Path(output);p=json.loads((output/'protocol.json').read_text());result={};c=p['criteria']
    for mode in MODES:
        history=json.loads((output/mode/'history.json').read_text());audit=json.loads((output/mode/'analysis.json').read_text())
        late=[r for r in history if r['metrics']['time']>=p['until']-p['late_duration']-1e-10]
        stable_ids=all(r['cell_order']==late[0]['cell_order'] for r in late)
        ranges=None
        if stable_ids:
            data=np.array([[[cell['attributes'][name] for name in FEATURES] for cell in r['cells']] for r in late])
            ranges=np.median(np.ptp(data,axis=0),axis=0).tolist()
        result[mode]={'quality':{'volume':audit['max_volume_error']<c['volume_max'],'radius':audit['min_radius']>=c['radius_min'],
                      'clipping':audit['max_clipping']==0,'boundary':audit['max_sampled_boundary']<c['boundary_max'],
                      'completed_development':audit['final_cells']==p['config']['max_cells'] and audit['active_divisions']==0},
                      'late_ids_stable':stable_ids,'late_mean_attribute_spread':np.mean([r['attribute_std'] for r in late],axis=0).tolist(),
                      'late_median_within_cell_attribute_range':ranges,'late_mean_axis_ratio':float(np.mean([r['metrics']['axis_ratio'] for r in late]))}
    report={'protocol':p,'arms':result,'all_quality_pass':all(all(r['quality'].values()) for r in result.values()),
            'interpretation':'Attribute dispersion and stable morphology are descriptive phenotypes, not an inferred number of identities. Constitutive mechanics is assumed. Persistence under perturbation and environmental changes, multiple seeds, and resolution checks remain necessary. No clustering was performed.'}
    write_json(output/'comparison.json',report)
    lines=['# Development without prescribed identities','',report['interpretation'],'',f'Numerical quality passes: **{report["all_quality_pass"]}**','',
           '| Mode | Late mean axis ratio | Late mean attribute SDs (in documented order) |','|---|---:|---|']
    for mode,r in result.items():lines.append(f'| {mode} | {r["late_mean_axis_ratio"]:.6f} | {r["late_mean_attribute_spread"]} |')
    (output/'RESULTS.md').write_text('\n'.join(lines)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,3,figsize=(12,6),layout='constrained')
    for mode in MODES:
        history=json.loads((output/mode/'history.json').read_text())
        times=[r['metrics']['time'] for r in history]
        spread=np.array([r['attribute_std'] for r in history])
        for i,name in enumerate(FEATURES):
            axes.flat[i].plot(times,spread[:,i],label=mode)
            axes.flat[i].set(title=name.replace('_',' '),ylabel='Cell-to-cell standard deviation')
        axes.flat[5].plot(times,[r['metrics']['axis_ratio'] for r in history],label=mode)
    axes.flat[5].set(title='Aggregate geometry',ylabel='Axis ratio')
    for ax in axes.flat:ax.set_xlabel('Developmental time');ax.grid(alpha=.2)
    axes.flat[0].legend();fig.suptitle('Continuous attributes without prescribed identities; dispersion is not a cell-type count')
    fig.savefig(output/'comparison.png',dpi=160);plt.close(fig)
    return report


def run(output,parallel=True):
    output=Path(output);p=json.loads((output/'protocol.json').read_text())
    state=json.loads((output/'status.json').read_text())
    if state['state']!='prepared' or hashlib.sha256((output/'protocol.json').read_bytes()).hexdigest()!=state['protocol_sha256']:raise ValueError('requires unchanged prepared study')
    for path,digest in p['source_sha256'].items():
        if hashlib.sha256(Path(path).read_bytes()).hexdigest()!=digest:raise ValueError('source changed since preparation')
    write_json(output/'status.json',{'state':'running','modes':list(MODES)})
    try:
        if parallel:
            with ProcessPoolExecutor(max_workers=2) as pool:list(pool.map(run_arm,[output]*2,MODES))
        else:
            for mode in MODES:run_arm(output,mode)
        report=compare(output);write_json(output/'status.json',{'state':'completed','all_quality_pass':report['all_quality_pass']});return report
    except Exception as error:
        write_json(output/'status.json',{'state':'failed','error':str(error)});raise


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=['prepare','run']);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();prepare(args.output) if args.action=='prepare' else run(args.output)
