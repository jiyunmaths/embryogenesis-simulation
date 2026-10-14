"""Prospective single-cleavage qualification; no fresh scientific histories."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
from multiprocessing import get_context
from pathlib import Path
import hashlib,json,time,fcntl
import numpy as np
from scipy.integrate import quad
from scipy.ndimage import label
from scipy.special import expit
import torch
from .model import Config,occupancy
from .attribute_development import AttributeSimulation
from .development_reference import ReferenceDevelopment
from .development_gpu import DevelopmentGpu,METHOD,SCHEMA,library
from .compensated_moving import numpy_scale,validate_correction,MovingSimulation
from .compensated_chemistry import NumpyStepper
from .cell_response_moving import observe
from .native_mechanics import NativeSimulation
from .resolution import _steps

LIMITS=dict(chemical_log_max=1e-5,polarity_abs_max=1e-5,relative_axis_max=1e-4,
    relative_volume_max=1e-5,relative_transport_max=1e-5,final_phi_abs_max=2e-5,
    ancestor_chemical_log_max=1e-5,whole_amount_relative_max=1e-5,
    event_occupancy_jump_max=1e-6,event_amount_jump_max=2e-14,event_target_jump_max=2e-14,
    effective_occupancy_absolute_jump_max=1e-14,dilution_amount_error_max=2e-14,
    native_prefix_phi_max=2e-5,native_prefix_chemical_max=1e-5,
    numerical_chemical_max=1e-11,volume_error_max=.05,radius_min=4.,boundary_max=.01,
    maximum_unmatched_observations=1)
TIME_LIMITS=dict(chemical_log_max=.01,polarity_abs_max=.01,relative_axis_max=.01,
    relative_volume_max=.005,relative_transport_max=.01,final_phi_abs_max=.02,
    ancestor_chemical_log_max=.01,whole_amount_relative_max=.005,maximum_unmatched_observations=1)


def digest(f):
    h=hashlib.sha256()
    with Path(f).open('rb') as x:
        for b in iter(lambda:x.read(1024*1024),b''):h.update(b)
    return h.hexdigest()


def read(f):return json.loads(Path(f).read_text())
def write(f,x):Path(f).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def pin(files):return {str(Path(f).resolve()):digest(f) for f in files}


def save_cpu(sim,path):
    sim.checkpoint(path)
    with np.load(path,allow_pickle=False) as z:out={k:z[k].copy() for k in z.files}
    out.update(phase_carry=sim.phase_carry,chemical_correction=sim.chemical_correction,
        chemical_method=np.array(METHOD),development_schema=np.array(SCHEMA),event_audits=np.array(json.dumps(sim.event_audits)),
        rounding=np.zeros(2,dtype=np.int32),precision_arm=np.array('phase_carry'))
    np.savez_compressed(path,**out)


def reference(path):
    sim=AttributeSimulation.restore(path);sim.__class__=ReferenceDevelopment
    with np.load(path,allow_pickle=False) as z:
        sim.phase_carry=z['phase_carry'].copy();sim.chemical_correction=z['chemical_correction'].copy();sim.event_audits=json.loads(str(z['event_audits']))
    sim.dilution_amount_error=0.;return sim


def verify(p):
    for file,h in {**p['source_sha256'],**p['input_sha256']}.items():
        if digest(file)!=h:raise ValueError('Pinned qualification evidence changed: '+file)


def prepare(root):
    root=Path(root).resolve()
    if root.exists():raise FileExistsError('Fresh qualification output required')
    root.mkdir(parents=True);jobs=[];files=[]
    for axis,direction in (('axial',[1,0,0]),('oblique',[1,2,3])):
        for branch in ('off','tension_adhesion','polarity','full'):
            for level,dt in (('coarse',.00375),('fine',.001875)):
                cfg=Config(grid=72,extent=2.24,dt=dt,steps=_steps(3.,dt),save_every=_steps(.03,dt),
                    max_cells=2,division_interval=100.,signal_dh=.55,signal_partition_noise=.001,
                    feedback=branch!='off',fate_tension=.25 if branch in ('tension_adhesion','full') else 0.,
                    fate_adhesion=.35 if branch in ('tension_adhesion','full') else 0.,
                    polarity_tension=.35 if branch in ('polarity','full') else 0.)
                sim=ReferenceDevelopment(cfg);r=np.sqrt(np.sum(sim.xyz**2,axis=0))
                sim.phi=expit(np.sqrt(2)*(.4-r)/cfg.interface_width)[None].astype(np.float32)
                sim.phase_carry=np.zeros_like(sim.phi,dtype=float)
                v=quad(lambda x:4*np.pi*x*x*occupancy(expit(np.sqrt(2)*(.4-x)/cfg.interface_width)),0,.4+30*cfg.interface_width,epsabs=1e-12,epsrel=1e-12)[0]
                sim.target[:]=v;sim.initial_volume=v;sim.activator[:]=1.03;sim.inhibitor[:]=.98
                if axis=='oblique':sim.polarity[:]=[.08,.03,-.02]
                name=axis+'_'+branch+'_'+level;folder=root/name;folder.mkdir();source=folder/'source.npz';save_cpu(sim,source);files.append(source)
                jobs.append(dict(key=name,axis=axis,direction=direction,branch=branch,level=level,dt=dt,source=str(source),config=asdict(sim.config)))
    accepted=Path('output/studies/moving-pulse-time-amount-response');ap=read(accepted/'protocol.json')
    job=next(j for j in ap['jobs'] if j['level']=='fine' and j['dose'] is None)
    mature=accepted/job['key']/'latest_state.npz';files.extend([accepted/'protocol.json',accepted/'summary.json',mature])
    names=('qualify_development_gpu.py','development_gpu.py','development_reference.py','gpu_development_kernels.cu',
        'gpu_backend.py','gpu_kernels.cu','gpu_precision_control.py','gpu_precision_kernels.cu','cuda_mechanics.cu','cuda_spatial_bench.cu',
        'compensated_chemistry.py','compensated_moving.py','attribute_development.py','model.py','polarity.py','transport.py','signaling.py','native_mechanics.py','native_mechanics.cpp','fast_mechanics.py','fast_mechanics.c')
    sources=[Path(__file__).with_name(n) for n in names];sources.extend([Path('tests/test_development_gpu.py'),Path('docs/development_gpu_qualification.md')])
    p=dict(jobs=jobs,duration=3.,interval=.03,native_prefix=.06,active_restart_time=.45,restart_steps=8,limits=LIMITS,time_limits=TIME_LIMITS,
        event_time_difference_dt_multiple=2,method=METHOD,mature_regression_source=str(mature),mature_regression_steps=8,
        gpu_uuid='GPU-2991ab3d-dcbd-41f3-eac1-133be278f56e',source_sha256=pin(sources),input_sha256=pin(files),
        histories=0,new_scientific_histories=0,scope='Sixteen prescribed single-cleavage synthetic cases: two axes, four mechanical branches and two timesteps. Analytic small cell, zero or modest specified polarity, nonuniform-in-time but spatially uniform starting chemistry, stochastic amount-conserving daughter partition. Backend/event/restart qualification only; no fresh repeated development, spatial convergence or identity evidence.')
    write(root/'protocol.json',p);write(root/'status.json',dict(state='prepared',total=16));return p


def row(sim,is_gpu):
    h=sim.observe(sim.time) if is_gpu else observe(sim,sim.time)
    active=sim.host.divisions if is_gpu else sim.divisions
    phi=sim.phi.cpu().numpy() if is_gpu else sim.phi
    h['nondividing_components']=[int(label(f>=.5)[1]) for i,f in enumerate(phi) if int(sim.ids[i]) not in active]
    h['dividing_ids']=sorted(active);return h


def audit_cpu(sim):
    v=sim.volumes();return dict(max_volume_error=float(abs(v/sim.target-1).max()),min_radius=float(((3*v/(4*np.pi))**(1/3)/sim.dx).min()),
        max_clipping=sim.clipped_fraction,dilution_amount_error=sim.dilution_amount_error)


def quality(a,p):
    c=p['limits'];return a['max_volume_error']<c['volume_error_max'] and a['min_radius']>=c['radius_min'] and a['max_clipping']==0 and a['dilution_amount_error']<=c['dilution_amount_error_max']


def run_cpu(args):
    root,job=args;root=Path(root);p=read(root/'protocol.json');verify(p);folder=root/job['key'];sim=reference(job['source']);sim.divide(0,job['direction'])
    history=[];maxima=dict(max_volume_error=0.,min_radius=1e100,max_clipping=0.,dilution_amount_error=0.);started=time.perf_counter()
    for k in range(_steps(p['duration'],job['dt'])+1):
        a=audit_cpu(sim)
        for key,v in a.items():maxima[key]=min(maxima[key],v) if key=='min_radius' else max(maxima[key],v)
        if not quality(a,p):raise FloatingPointError('CPU development quality failure '+job['key'])
        if k%job['config']['save_every']==0:history.append(row(sim,False))
        if k<_steps(p['duration'],job['dt']):sim.step()
    write(folder/'cpu-history.json',history);save_cpu(sim,folder/'cpu-final.npz')
    # Historical native arithmetic is an additional short physical reference,
    # not a substitute for corrected full development.
    native=AttributeSimulation.restore(job['source']);native.__class__=NativeSimulation;native.native_threads=1;native.divide(0,job['direction'])
    prefix=[]
    for k in range(_steps(p['native_prefix'],job['dt'])+1):
        if k%job['config']['save_every']==0:prefix.append(observe(native,native.time))
        if k<_steps(p['native_prefix'],job['dt']):native.step()
    write(folder/'native-prefix.json',prefix);np.savez_compressed(folder/'native-prefix-phi.npz',phi=native.phi)
    result=dict(job=job,passed=True,protocol_sha256=digest(root/'protocol.json'),audit=maxima,events=sim.event_audits,seconds=time.perf_counter()-started,
        evidence_sha256=pin(folder/f for f in ('cpu-history.json','cpu-final.npz','native-prefix.json','native-prefix-phi.npz')))
    write(folder/'cpu-result.json',result);return result


def exact(a,b):
    for name in ('phi','activator','inhibitor','polarity','target','geometry','phase_carry','rounding','chemical_correction'):
        if not torch.equal(getattr(a,name),getattr(b,name)):raise RuntimeError('Development restart changed '+name)
    for name in ('ids','parents','due'):
        x,y=(a.ids,b.ids) if name=='ids' else (getattr(a.host,name),getattr(b.host,name));np.testing.assert_array_equal(x,y)
    assert a.time==b.time and a.step_number==b.step_number and a.host.next_id==b.host.next_id
    for name in ('rng','fate_rng','signal_rng'):assert getattr(a.host,name).bit_generator.state==getattr(b.host,name).bit_generator.state
    assert a.host.divisions==b.host.divisions and a.host.lineage==b.host.lineage and a.event_audits==b.event_audits


def run_gpu(root,job,p):
    folder=root/job['key'];sim=DevelopmentGpu.restore(job['source']);sim.divide(0,job['direction']);sim.checkpoint(folder/'gpu-started.npz')
    history=[];maxima=dict(max_volume_error=0.,min_radius=1e100,max_clipping=0.,dilution_amount_error=0.);started=time.perf_counter()
    restarts={};clone=None;label_name=None;remaining=0;independent_error=0.
    for k in range(_steps(p['duration'],job['dt'])+1):
        a=sim.audit()
        for key in maxima:maxima[key]=min(maxima[key],a[key]) if key=='min_radius' else max(maxima[key],a[key])
        if not quality(a,p):raise FloatingPointError('GPU development quality failure '+job['key'])
        if k%job['config']['save_every']==0:history.append(row(sim,True))
        if abs(sim.time-p['native_prefix'])<1e-10:np.savez_compressed(folder/'gpu-prefix-phi.npz',phi=sim.phi.cpu().numpy())
        if abs(sim.time-p['active_restart_time'])<1e-10 or (len(sim.ids)==2 and 'post' not in restarts and label_name!='post'):
            label_name='active' if len(sim.ids)==1 else 'post'
            if clone is not None:raise RuntimeError('Unexpected overlapping restart windows')
            cp=folder/(label_name+'-restart.npz');sim.checkpoint(cp);clone=DevelopmentGpu.restore(cp);exact(sim,clone);remaining=p['restart_steps']
        if k==_steps(p['duration'],job['dt']):break
        if k<8:
            u=torch.stack((sim.activator,sim.inhibitor)).cpu().numpy();e=sim.chemical_correction.cpu().numpy();old=sim.geometry[:,0].cpu().numpy().copy()
            expected,correction=NumpyStepper(sim.matrices()[1].cpu().numpy(),job['dt'],sim.config.signal_beta,sim.config.signal_da,sim.config.signal_dh)(u,e)
        sim.step()
        if k<8:
            expected,correction=numpy_scale(expected,correction,old/sim.geometry[:,0].cpu().numpy());actual=torch.stack((sim.activator,sim.inhibitor)).cpu().numpy()
            independent_error=max(independent_error,float(abs(np.log(actual/expected)).max()))
            if independent_error>p['limits']['numerical_chemical_max']:raise RuntimeError('Independent moving chemical check failed')
        if clone is not None:
            clone.step();exact(sim,clone);remaining-=1
            if remaining==0:restarts[label_name]=p['restart_steps'];del clone;clone=None;label_name=None
    if restarts!={'active':8,'post':8} or len(sim.event_audits)!=1 or sim.host.divisions:raise RuntimeError('Unfinished event/restart qualification')
    write(folder/'gpu-history.json',history);sim.checkpoint(folder/'gpu-final.npz')
    result=dict(job=job,passed=True,protocol_sha256=digest(root/'protocol.json'),audit=maxima,events=sim.event_audits,seconds=time.perf_counter()-started,
        exact_restart_steps=restarts,independent_chemical_steps=8,independent_chemical_log_max=independent_error,
        evidence_sha256=pin(folder/f for f in ('gpu-history.json','gpu-final.npz','gpu-started.npz','active-restart.npz','post-restart.npz','gpu-prefix-phi.npz')))
    write(folder/'gpu-result.json',result);del sim;torch.cuda.empty_cache();return result


def compare_histories(a,b):
    keys=('chemical_log_max','polarity_abs_max','relative_axis_max','relative_volume_max','relative_transport_max','ancestor_chemical_log_max','whole_amount_relative_max')
    errors=dict.fromkeys(keys,0.);unmatched=[]
    assert len(a)==len(b)
    for x,y in zip(a,b):
        assert abs(x['time']-y['time'])<1e-9
        ux,uy=np.asarray(x['chemistry']),np.asarray(y['chemistry']);vx,vy=np.asarray(x['volumes']),np.asarray(y['volumes'])
        errors['relative_axis_max']=max(errors['relative_axis_max'],abs(x['axis_ratio']/y['axis_ratio']-1))
        errors['whole_amount_relative_max']=max(errors['whole_amount_relative_max'],float(abs((ux@vx)/(uy@vy)-1).max()))
        errors['ancestor_chemical_log_max']=max(errors['ancestor_chemical_log_max'],float(abs(np.log((ux@vx/vx.sum())/(uy@vy/vy.sum()))).max()))
        if x['ids']!=y['ids']:unmatched.append(x['time']);continue
        dx,dy=np.asarray(x['delta']),np.asarray(y['delta']);den=float(np.linalg.norm(dy))
        for key,v in dict(chemical_log_max=float(abs(np.log(ux/uy)).max()),polarity_abs_max=float(abs(np.asarray(x['polarity'])-y['polarity']).max()),
            relative_volume_max=float(abs(vx/vy-1).max()),relative_transport_max=float(np.linalg.norm(dx-dy)/den) if den else float(np.linalg.norm(dx-dy))).items():errors[key]=max(errors[key],v)
    errors['maximum_unmatched_observations']=len(unmatched);return errors,unmatched


def assess(root):
    root=Path(root);p=read(root/'protocol.json');verify(p);rows=[];timesteps=[];evidence=[]
    for job in p['jobs']:
        folder=root/job['key'];reports=[read(folder/(arm+'-result.json')) for arm in ('gpu','cpu')]
        for report in reports:
            assert report['passed'] and report['protocol_sha256']==digest(root/'protocol.json')
            for f,h in report['evidence_sha256'].items():assert digest(f)==h
        a,b=[read(folder/(arm+'-history.json')) for arm in ('gpu','cpu')];error,unmatched=compare_histories(a,b)
        with np.load(folder/'gpu-final.npz') as x,np.load(folder/'cpu-final.npz') as y:
            error['final_phi_abs_max']=float(abs(x['phi']-y['phi']).max());assert np.array_equal(x['ids'],y['ids'])
        gpu_event,cpu_event=[r['events'][0] for r in reports];time_error=abs(gpu_event['time']-cpu_event['time'])
        passed=all(v<=p['limits'][k] for k,v in error.items()) and time_error<=2*job['dt']+1e-12
        for report,h in zip(reports,(a,b)):
            passed=passed and len(report['events'])==1 and all(r['boundary_occupancy']<.01 and all(c==1 for c in r['nondividing_components']) for r in h)
            for ev in report['events']:
                passed=passed and ev['daughter_components']==[1,1]
                for key,limit in [('amount_relative_jump','event_amount_jump_max'),('target_relative_jump','event_target_jump_max'),('volume_relative_jump','event_occupancy_jump_max'),('effective_occupancy_absolute_jump','effective_occupancy_absolute_jump_max')]:passed=passed and ev[key]<=p['limits'][limit]
        native=read(folder/'native-prefix.json');prefix_error,_=compare_histories(a[:len(native)],native)
        with np.load(folder/'gpu-prefix-phi.npz') as x,np.load(folder/'native-prefix-phi.npz') as y:native_phi=float(abs(x['phi']-y['phi']).max())
        passed=passed and native_phi<=p['limits']['native_prefix_phi_max'] and prefix_error['chemical_log_max']<=p['limits']['native_prefix_chemical_max']
        rows.append(dict(case=job['key'],errors=error,event_time_error=time_error,unmatched_times=unmatched,native_prefix_phi_error=native_phi,native_prefix_chemical_error=prefix_error['chemical_log_max'],passed=bool(passed)))
        evidence.extend(folder/f for f in ('gpu-result.json','cpu-result.json','gpu-history.json','cpu-history.json','gpu-final.npz','cpu-final.npz'))
    for axis in ('axial','oblique'):
        for branch in ('off','tension_adhesion','polarity','full'):
            names=[axis+'_'+branch+'_'+l for l in ('coarse','fine')];h=[read(root/n/'gpu-history.json') for n in names];err,unmatched=compare_histories(*h)
            with np.load(root/names[0]/'gpu-final.npz') as a,np.load(root/names[1]/'gpu-final.npz') as b:err['final_phi_abs_max']=float(abs(a['phi']-b['phi']).max())
            times=[read(root/n/'gpu-result.json')['events'][0]['time'] for n in names];difference=abs(times[0]-times[1])
            passed=all(v<=p['time_limits'][k] for k,v in err.items()) and difference<=.0075+1e-12
            timesteps.append(dict(axis=axis,branch=branch,errors=err,unmatched_times=unmatched,event_time_error=difference,passed=bool(passed)))
    result=dict(passed=all(r['passed'] for r in rows+timesteps),independent_review_pending=True,protocol_sha256=digest(root/'protocol.json'),backend_rows=rows,timestep_rows=timesteps,new_scientific_histories=0,evidence_sha256=pin(evidence),scope=p['scope'])
    write(root/'summary.json',result);return result


def mature_regression(root,p):
    path=p['mature_regression_source'];a=MovingSimulation.restore(path);b=DevelopmentGpu(AttributeSimulation.restore(path),_restore=True)
    for name in ('phase_carry','rounding','chemical_correction'):getattr(b,name).copy_(getattr(a,name))
    for _ in range(p['mature_regression_steps']):
        a.step();b.step()
        for name in ('phi','activator','inhibitor','polarity','target','geometry','phase_carry','rounding','chemical_correction'):
            if not torch.equal(getattr(a,name),getattr(b,name)):raise RuntimeError('Mature regression differs: '+name)
    write(root/'mature-regression.json',dict(passed=True,exact_steps=8,source=path,source_sha256=digest(path)));del a,b;torch.cuda.empty_cache()


def _run(root,workers):
    root=Path(root).resolve();p=read(root/'protocol.json');verify(p);torch.set_num_threads(1);started=time.perf_counter()
    try:
        lib=library();write(root/'runtime.json',dict(gpu=torch.cuda.get_device_name(0),torch=torch.__version__,binary=lib.binary_path,binary_sha256=digest(lib.binary_path)))
        mature_regression(root,p)
        with ProcessPoolExecutor(max_workers=workers,mp_context=get_context('spawn')) as pool:
            futures=[pool.submit(run_cpu,(str(root),j)) for j in p['jobs']]
            for number,job in enumerate(p['jobs'],1):
                for future in futures:
                    if future.done():future.result()
                write(root/'status.json',dict(state='running',case=job['key'],completed_gpu=number-1,total=16));r=run_gpu(root,job,p)
                print(job['key'],r['events'][0]['time'],r['audit'],flush=True)
            for future in futures:future.result()
        r=assess(root);r['seconds']=time.perf_counter()-started;write(root/'timing.json',dict(seconds=r['seconds']));write(root/'status.json',dict(state='completed' if r['passed'] else 'completed_with_unresolved_checks',passed=r['passed'],independent_review_pending=True,total=16))
        return r
    except Exception as e:write(root/'status.json',dict(state='failed',error=str(e)));raise


def run(root,workers=4):
    root=Path(root).resolve()
    if not 1<=workers<=4:raise ValueError('Use one to four CPU reference workers')
    with (root/'coordinator.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        if read(root/'status.json')['state']!='prepared':raise ValueError('Qualification requires an unrun prepared directory; preserve previous attempts')
        return _run(root,workers)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('prepare','run','assess'));parser.add_argument('--output',type=Path,default=Path('output/studies/development-gpu-cleavage-qualification'));parser.add_argument('--workers',type=int,default=4);args=parser.parse_args()
    result=prepare(args.output) if args.action=='prepare' else run(args.output,args.workers) if args.action=='run' else assess(args.output)
    print({k:v for k,v in result.items() if k in ('passed','seconds','new_scientific_histories')})
