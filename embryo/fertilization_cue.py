"""Finite-duration directional zygote polarity cue; no prescribed identities."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import time
import numpy as np
from .attribute_development import AttributeSimulation
from .model import Config,occupancy
from .feedback_long import digest,save_checkpoint
from .resolution import write_json,_steps


def cue_increment(t,dt,strength,duration):
    """Exact integral of (2 strength / T) sin²(pi t/T), supported on [0,T]."""
    if duration<=0 or strength<0 or dt<=0:raise ValueError('Invalid cue parameters')
    a,b=np.clip([t,t+dt],0.,duration)
    return float(strength*((b-a)/duration-(np.sin(2*np.pi*b/duration)-np.sin(2*np.pi*a/duration))/(2*np.pi)))


class CueSimulation(AttributeSimulation):
    def __init__(self,config=None,mode=None,*,strength=0.,duration=.75,direction=(1.,0.,0.),clamp=False):
        super().__init__(config,mode)
        d=np.asarray(direction,float)
        if d.shape!=(3,) or not np.isfinite(d).all() or np.linalg.norm(d)==0 or not np.isfinite([strength,duration]).all() or strength<0 or duration<=0:raise ValueError('Invalid cue')
        self.cue=dict(strength=float(strength),duration=float(duration),direction=(d/np.linalg.norm(d)).tolist(),clamp=bool(clamp))
        self.cue_delivered=0.
        if duration>=self.config.division_interval*(1-self.config.cycle_jitter):raise ValueError('Cue must end before earliest scheduled cleavage')
        if clamp:self.config.signal_partition_noise=0.

    def mechanical_step(self):
        c=self.cue;increment=cue_increment(self.time,self.config.dt,c['strength'],c['duration'])
        if increment:
            if len(self.ids)!=1 or self.divisions:raise RuntimeError('Cue overlaps cleavage')
            self.polarity[0]+=increment*np.array(c['direction']);self.cue_delivered+=increment
            if np.linalg.norm(self.polarity[0])>1.:raise RuntimeError('Cue exceeds polarity norm bound; reduce strength')
        # Cue is added after intrinsic polarity evolution and before mechanics:
        # a first-order split update, assessed by timestep refinement.
        enabled=self.config.polarity_enabled
        # A zero-coupling control must use the same scalar mechanics kernel
        # regardless of whether its diagnostic polarity vector is nonzero.
        if self.config.polarity_tension==0:self.config.polarity_enabled=False
        try:super().mechanical_step()
        finally:self.config.polarity_enabled=enabled

    def step(self,**kwargs):
        if self.cue['clamp']:
            if kwargs:raise ValueError('Clamp is configured by the experiment')
            kwargs=dict(prescribed_signals=(np.ones(len(self.ids)),np.ones(len(self.ids))))
        super().step(**kwargs)

    def checkpoint(self,path):
        super().checkpoint(path)
        with np.load(path) as d:payload={k:d[k].copy() for k in d.files}
        payload['fertilization_cue']=np.array(json.dumps(dict(parameters=self.cue,delivered=self.cue_delivered)))
        np.savez_compressed(path,**payload)

    @classmethod
    def restore(cls,path):
        sim=super().restore(path)
        with np.load(path) as d:metadata=json.loads(str(d['fertilization_cue']))
        sim.cue=metadata['parameters'];sim.cue_delivered=metadata['delivered']
        return sim


def observe(sim):
    h=occupancy(sim.phi);union=np.minimum(h.sum(axis=0),1.);w=union.sum()
    center=np.einsum('ijk,dijk->d',union,sim.xyz)/w
    offsets=(sim.xyz-center[:,None,None,None]).reshape(3,-1)
    cov=(offsets*union.ravel())@offsets.T/w;values,vectors=np.linalg.eigh(cov)
    direction=np.array(sim.cue['direction']);gap=float((values[-1]-values[-2])/values[-1])
    m=sim.volumes();centers=sim.centers();r=centers-np.average(centers,weights=m,axis=0)
    la=np.log(sim.activator);fluct=la-np.average(la,weights=m);dipole=np.sum((m*fluct)[:,None]*r,axis=0)
    scale=float(np.sum(m*abs(fluct)*np.linalg.norm(r,axis=1)));strength=float(np.linalg.norm(dipole)/max(scale,1e-30))
    alignment=float(dipole@direction/np.linalg.norm(dipole)) if strength>.01 and np.std(la)>.001 else None
    return dict(time=sim.time,ids=sim.ids.tolist(),chemistry=[sim.activator.tolist(),sim.inhibitor.tolist()],polarity=sim.polarity.tolist(),centers=centers.tolist(),volumes=m.tolist(),centroid=center.tolist(),covariance=cov.tolist(),
        shape_axis_eigengap=gap,shape_axis_alignment=float(abs(vectors[:,-1]@direction)) if gap>.01 else None,
        chemical_dipole_alignment=alignment,chemical_dipole_strength=strength,log_activator_sd=float(np.std(la)),
        cue_delivered=sim.cue_delivered,metrics=sim.metrics())


def prepare(root):
    root=Path(root).resolve()
    if root.exists():raise FileExistsError(root)
    c=Config(grid=72,extent=2.24,dt=.00375,steps=24000,save_every=40,seed=7)
    c.validate()
    def job(key,strength=0.,direction=(1,0,0),duration=.75,polar=True,clamp=False):
        return dict(key=key,strength=strength,direction=list(direction),duration=duration,polar=polar,clamp=clamp)
    jobs=[job('baseline'),job('weak_x',.1),job('strong_x',.4),job('strong_diagonal',.4,(1,1,1)),job('strong_minus_x',.4,(-1,0,0)),job('long_x',.4,duration=1.5),job('no_polar_baseline',polar=False),job('no_polar_cue',.4,polar=False),job('clamped_baseline',clamp=True),job('clamped_cue',.4,clamp=True)]
    root.mkdir(parents=True)
    files=[Path(__file__),*[Path(__file__).with_name(f) for f in ('attribute_development.py','model.py','polarity.py','transport.py','signaling.py','feedback_long.py','resolution.py')]]
    p=dict(config=asdict(c),jobs=jobs,duration=90.,interval=.15,screen_duration=1.5,screen_dts=[.0075,.00375,.001875],
        criteria=dict(volume_max=.05,radius_min=4.,boundary_max=.01,clipping_max=0.,screen_polarity_abs_error=.002,screen_relative_covariance_error=.002,screen_centroid_abs_error=.002),
        cue='Add finite-time polarity forcing only in undivided zygote. g(t)=2/T sin²(pi t/T) for 0<t<T, zero otherwise. Strength is integrated forcing, not final polarity. Positive direction is toward nominal entry point and lowers tension there under existing polarity law.',
        inheritance='Existing division prolongation copies maternal polarity into both daughters; no cue is applied after T. Not determinant segregation or cortical flow.',
        scope='Exploratory seed-7 ten-arm development from analytic zygote to t=90, 16-cell cap. No assigned fates. Directions probe grid sensitivity but one history cannot establish rotational invariance or probability of axis selection. Cue amplitude/duration in model units, not calibrated sperm measurements.',
        controls='No cue; strengths .1/.4; opposite and diagonal orientations; durations .75/1.5 at equal integrated drive; matched polarity-tension-zero pair; matched externally uniform chemistry pair (includes zero partition noise). Clamp is an intervention, not selective molecular inhibition.',
        analysis='Centered shape covariance with axis undefined at leading eigengap <=.01; chemical dipole alignment undefined when normalized strength <=.01 or log-activator SD <=.001. Separate translation, shape magnitude, chemical contrast, cleavage axes, and orientation. No identity claim.',
        scheduling='Run short screen now on one process; full pilot waits for both previous moving studies to complete successfully. Gate prevents silent continuation after numerical failure.',
        dependencies=[str(Path(x).resolve()) for x in ('outputs/cell-response-moving-refined','outputs/cell-exchange-moving')],
        dependency_protocol_sha256={str(Path(x).resolve()):digest(Path(x)/'protocol.json') for x in ('outputs/cell-response-moving-refined','outputs/cell-exchange-moving')},
        source_sha256={str(f.resolve()):digest(f) for f in files})
    write_json(root/'protocol.json',p);write_json(root/'status.json',dict(state='prepared',completed=0,total=len(jobs)))


def simulate(root,p,job,dt,until):
    folder=root/job['key'];folder.mkdir(exist_ok=True);ph=digest(root/'protocol.json');ck=folder/'latest_state.npz'
    if (folder/'result.json').exists():
        result=json.loads((folder/'result.json').read_text())
        if result['protocol_sha256']!=ph:raise ValueError('Result protocol mismatch')
        return result
    if ck.exists():
        sim=CueSimulation.restore(ck)
        with np.load(ck) as d:meta=json.loads(str(d['long_experiment']))
        if meta['protocol_hash']!=ph or meta['job']!=job:raise ValueError('Restart mismatch')
        audit,history=meta['audit'],meta['history']
    else:
        config=Config(**dict(p['config'],dt=dt,steps=_steps(until,dt),save_every=_steps(p['interval'],dt)))
        if not job['polar']:config.polarity_tension=0.
        sim=CueSimulation(config,'direct',strength=job['strength'],duration=job['duration'],direction=job['direction'],clamp=job['clamp'])
        audit=dict(max_volume_error=0.,min_radius=1e100,max_clipping=0.,max_boundary=0.,elapsed_seconds=0.);history=[]
    old=audit['elapsed_seconds'];clock=time.monotonic();last=history[-1]['time'] if history else -1
    try:
        while sim.step_number<=sim.config.steps:
            v=sim.volumes();audit['max_volume_error']=max(audit['max_volume_error'],float(abs(v/sim.target-1).max()));audit['min_radius']=min(audit['min_radius'],float(((3*v/(4*np.pi))**(1/3)/sim.dx).min()));audit['max_clipping']=max(audit['max_clipping'],sim.clipped_fraction)
            if not np.isfinite(sim.activator).all() or not np.isfinite(sim.inhibitor).all() or np.any(sim.activator<=0) or np.any(sim.inhibitor<=0):raise RuntimeError('Chemical positivity/finiteness failure')
            if audit['max_volume_error']>=.05 or audit['min_radius']<4 or audit['max_clipping']>0:raise RuntimeError('Numerical quality failure')
            if sim.step_number%sim.config.save_every==0 and sim.time>last+1e-9:
                row=observe(sim);history.append(row);last=sim.time;audit['max_boundary']=max(audit['max_boundary'],row['metrics']['boundary_occupancy'])
                if audit['max_boundary']>=.01:raise RuntimeError('Boundary quality failure')
                audit['elapsed_seconds']=old+time.monotonic()-clock
                write_json(folder/'history.json',history);write_json(folder/'status.json',dict(state='running',time=sim.time,audit=audit))
                if sim.step_number%_steps(1.5,dt)==0 or sim.step_number==sim.config.steps:save_checkpoint(sim,ck,audit,history,ph,job)
            if sim.step_number==sim.config.steps:break
            sim.step()
        result=dict(protocol_sha256=ph,job=job,dt=dt,until=until,audit=audit,quality_pass=True,final=history[-1],lineage=sim.lineage,completed_cleavage=bool(len(sim.ids)==sim.config.max_cells and not sim.divisions))
        write_json(folder/'result.json',result);write_json(folder/'status.json',dict(state='completed',time=sim.time));return result
    except Exception as exc:
        write_json(folder/'status.json',dict(state='failed',time=sim.time,error=str(exc)));raise


def screen(root,p):
    results=[]
    for dt in p['screen_dts']:
        job=dict(p['jobs'][2],key=f'screen-dt-{dt:g}')
        results.append(simulate(root,p,job,dt,p['screen_duration']))
    finest=results[-1]['final'];rows=[]
    for result in results[:-1]:
        final=result['final']
        rows.append(dict(dt=result['dt'],polarity_error=float(np.max(abs(np.array(final['polarity'])-finest['polarity']))),covariance_error=float(np.linalg.norm(np.array(final['covariance'])-finest['covariance'])/np.linalg.norm(finest['covariance'])),centroid_error=float(np.linalg.norm(np.array(final['centroid'])-finest['centroid']))))
    passed=all(r['polarity_error']<.002 and r['covariance_error']<.002 and r['centroid_error']<.002 for r in rows)
    passed=passed and all(abs(r['final']['cue_delivered']-.4)<1e-12 and len(r['final']['ids'])==1 and not any(x.get('division_start') is not None for x in r['lineage']) for r in results)
    report=dict(passed=passed,comparisons=rows,scope='Pre-cleavage strong-x cue screen only; does not establish developmental or spatial convergence.')
    write_json(root/'screen.json',report);return passed


def run(root):
    root=Path(root);p=json.loads((root/'protocol.json').read_text())
    def verify():
        for f,h in p['source_sha256'].items():
            if digest(f)!=h:raise ValueError('Source changed: '+f)
        for d,h in p['dependency_protocol_sha256'].items():
            if digest(Path(d)/'protocol.json')!=h:raise ValueError('Dependency changed')
    verify()
    try:
        write_json(root/'status.json',dict(state='screening',completed=0,total=10))
        if not screen(root,p):
            write_json(root/'status.json',dict(state='blocked',reason='Cue timestep screen failed'));return
        write_json(root/'status.json',dict(state='queued',completed=0,total=10,reason='Cue screen passed; waiting for prior moving studies'))
        while True:
            statuses=[json.loads((Path(d)/'status.json').read_text()) for d in p['dependencies']]
            if any(s['state'] in ('failed','blocked') for s in statuses):
                write_json(root/'status.json',dict(state='blocked',reason='Prior moving study failed'));return
            if all(s['state']=='completed' for s in statuses):
                ref=json.loads((Path(p['dependencies'][0])/'refinement.json').read_text())
                if not ref['passed']:
                    write_json(root/'status.json',dict(state='blocked',reason='Response refinement failed'));return
                break
            time.sleep(30)
        verify();results=[]
        for job in p['jobs']:
            write_json(root/'status.json',dict(state='running',completed=len(results),total=10,current=job['key'],workers=1))
            results.append(simulate(root,p,job,p['config']['dt'],p['duration']))
            write_json(root/'results.json',results)
        summary={r['job']['key']:dict(late_mean_log_activator_sd=float(np.mean([x['log_activator_sd'] for x in json.loads((root/r['job']['key']/'history.json').read_text()) if x['time']>=75])),completed_cleavage=r['completed_cleavage'],final=r['final']) for r in results}
        write_json(root/'comparison.json',dict(arms=summary,scope=p['scope'],warning='Descriptive single-history pilot; full developmental timestep and spatial convergence and multi-history orientation tests remain outstanding.'))
        write_json(root/'status.json',dict(state='completed',completed=10,total=10,all_completed_cleavage=all(r['completed_cleavage'] for r in results)))
    except Exception as exc:
        write_json(root/'status.json',dict(state='failed',error=str(exc)));raise


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('command',choices=('prepare','run'));parser.add_argument('--output',type=Path,default=Path('outputs/fertilization-cue'))
    a=parser.parse_args();prepare(a.output) if a.command=='prepare' else run(a.output)
