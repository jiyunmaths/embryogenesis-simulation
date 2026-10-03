"""Short saved-reference replay including GPU spatial reductions and contacts.

This is a transfer-heavy precision assay, not a production GPU integrator.
"""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from .cuda_benchmark import CudaMechanicalPilot
from .cuda_spatial_bench import library
from .cell_response import pulse
from .cell_response_moving import observe
from .feedback_long import digest
from .resolution import write_json

_LIB=None


class SpatialPrecisionPilot(CudaMechanicalPilot):
    def contacts(self):
        global _LIB
        if _LIB is None:_LIB=library()
        n=len(self.ids);geom=np.empty((n,7));contact=np.empty((n,n),np.float32)
        attraction=np.empty(self.phi.shape);times=np.empty(4,np.float32)
        _,adhesion=self.material_coefficients()
        code=_LIB.spatial_benchmark(self.phi,adhesion,n,self.config.grid,self.dx,self.config.extent,
                                   1,geom,contact,attraction,times)
        if code:raise RuntimeError(_LIB.cuda_mechanics_error().decode())
        # Keep the diagnostic exposure scalar from the original implementation.
        _,exposure=super().contacts()
        self._spatial_phi=self.phi;self._spatial_geometry=geom
        graph=contact.astype(float)*self.dx**3;np.fill_diagonal(graph,0)
        return graph,exposure

    def native_exposure(self):
        if getattr(self,'_spatial_phi',None) is not self.phi:self.contacts()
        return self._spatial_geometry[:,4:].copy()

    def volumes(self):
        if getattr(self,'_spatial_phi',None) is self.phi:return self._spatial_geometry[:,0].copy()
        return super().volumes()

    def centers(self):
        if getattr(self,'_spatial_phi',None) is self.phi:return self._spatial_geometry[:,1:4].copy()
        return super().centers()


def run(output,steps=160):
    if steps<40 or steps%40:raise ValueError('Use a positive multiple of 40 steps')
    root=Path(output);root.mkdir(parents=True,exist_ok=False)
    source=Path('outputs/cell-exchange-response-moving/unexchanged')
    jobs=[('unexchanged_control',None),('unexchanged_cell-27_factor-0.9',27)]
    criteria=dict(chemical_max_log=1e-5,polarity_max_abs=1e-5,volume_max_relative=1e-5,axis_max_relative=1e-4)
    records=[];paths={}
    files=[source/'source.npz']+[source/key/'history.json' for key,_ in jobs]
    write_json(root/'protocol.json',dict(steps=steps,dt=.00375,precisions=[64,32],criteria=criteria,
        normalized_response_limit=.001,input_sha256={str(f.resolve()):digest(f) for f in files},
        source_sha256={str(p.resolve()):digest(p) for p in [Path(__file__),Path(__file__).with_name('cuda_spatial_bench.py'),Path(__file__).with_name('cuda_spatial_bench.cu'),Path(__file__).with_name('cuda_benchmark.py'),Path(__file__).with_name('cuda_mechanics.cu')]},
        scope='Short prefix only, GPU contact/reduction/cue arithmetic plus GPU mechanics. Chemistry CPU, post-mechanical volumes CPU until next contact refresh. No cleavage; not a fully resident GPU timestep or full-horizon scientific validation.'))
    for precision in (64,32):
        for key,cell in jobs:
            sim=SpatialPrecisionPilot.restore(source/'source.npz');sim.precision=precision;sim.native_threads=1
            start_time=sim.time
            if cell is not None:sim.activator,sim.inhibitor=pulse(np.array([sim.activator,sim.inhibitor]),list(sim.ids).index(cell),.9)
            reference=json.loads((source/key/'history.json').read_text());history=[]
            errors=dict.fromkeys(criteria,0.);max_clipping=0.;max_volume=0.;started=time.perf_counter()
            for step in range(steps+1):
                if step%40==0:
                    row=observe(sim,sim.time-start_time);ref=reference[step//40]
                    if abs(row['time']-ref['time'])>1e-8 or row['ids']!=ref['ids']:raise ValueError('Reference mismatch')
                    current=dict(chemical_max_log=float(abs(np.log(np.array(row['chemistry'])/ref['chemistry'])).max()),polarity_max_abs=float(abs(np.array(row['polarity'])-ref['polarity']).max()),volume_max_relative=float(abs(np.array(row['volumes'])/ref['volumes']-1).max()),axis_max_relative=abs(row['axis_ratio']/ref['axis_ratio']-1))
                    for k,v in current.items():errors[k]=max(errors[k],v)
                    history.append(row)
                    status=dict(state='running',precision=precision,job=key,steps=step,total_steps=steps,errors=errors)
                    write_json(root/'status.json',status);print(json.dumps(status),flush=True)
                if step==steps:break
                sim.step();max_clipping=max(max_clipping,sim.clipped_fraction)
                max_volume=max(max_volume,float(abs(sim.volumes()/sim.target-1).max()))
                if not np.isfinite(sim.phi).all() or not np.isfinite(sim.activator).all() or np.any(sim.activator<=0) or np.any(sim.inhibitor<=0):raise RuntimeError('Invalid dynamics')
            paths[(precision,key)]=np.array([r['chemistry'] for r in history])
            record=dict(precision=precision,job=key,steps=steps,model_duration=sim.time-start_time,errors=errors,
                max_clipping=max_clipping,max_volume_error=max_volume,wall_seconds=time.perf_counter()-started,
                passed=all(errors[k]<=v for k,v in criteria.items()) and max_clipping==0 and max_volume<.05)
            records.append(record);write_json(root/f'{key}-{precision}-history.json',history);write_json(root/'runs.json',records)
    response=[]
    refs={key:np.array([r['chemistry'] for r in json.loads((source/key/'history.json').read_text())[:steps//40+1]]) for key,_ in jobs}
    ref_response=np.log(refs[jobs[1][0]]/refs[jobs[0][0]])/abs(np.log(.9))
    for precision in (64,32):
        value=np.log(paths[(precision,jobs[1][0])]/paths[(precision,jobs[0][0])])/abs(np.log(.9))
        response.append(dict(precision=precision,max_normalized_response_error=float(abs(value-ref_response).max())))
    passed=all(r['passed'] for r in records) and all(r['max_normalized_response_error']<=.001 for r in response)
    write_json(root/'comparison.json',dict(passed=passed,runs=records,response_errors=response,scope='Short prefix only; does not measure full recovery or long-term stability.'))
    write_json(root/'status.json',dict(state='completed',passed=passed))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True);p.add_argument('--steps',type=int,default=160);a=p.parse_args();run(a.output,a.steps)
