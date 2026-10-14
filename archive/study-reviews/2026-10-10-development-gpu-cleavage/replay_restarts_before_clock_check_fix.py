"""Replay each saved live-furrow state through the complete remaining event.

This independently repeats checkpoint continuation using the pinned candidate
solver; it is a restart gate, not an independent physical model.
"""
from pathlib import Path
import hashlib,json,time
import numpy as np
import torch
from embryo.development_gpu import DevelopmentGpu

ROOT=Path('output/studies/development-gpu-cleavage-qualification-v2')
HERE=Path(__file__).parent
def sha(f):return hashlib.sha256(Path(f).read_bytes()).hexdigest()


def main():
    p=json.loads((ROOT/'protocol.json').read_text());torch.set_num_threads(1);rows=[];watched={}
    for f,h in {**p['source_sha256'],**p['input_sha256']}.items():assert sha(f)==h;watched[f]=h
    for job in p['jobs']:
        folder=ROOT/job['key'];files=[folder/f for f in ('active-restart.npz','gpu-final.npz','gpu-result.json')]
        for f in files:watched[str(f.resolve())]=sha(f)
        result=json.loads(files[2].read_text());assert result['passed'] and result['exact_restart_steps']=={'active':8,'post':8}
        sim=DevelopmentGpu.restore(files[0]);assert sim.time==.45 and len(sim.ids)==1 and len(sim.host.divisions)==1
        before=sim.step_number;started=time.perf_counter()
        while sim.time<p['duration']-1e-12:sim.step()
        expected=DevelopmentGpu.restore(files[1])
        for name in ('phi','activator','inhibitor','polarity','target','geometry','phase_carry','rounding','chemical_correction'):
            assert torch.equal(getattr(sim,name),getattr(expected,name)),(job['key'],name)
        np.testing.assert_array_equal(sim.ids,expected.ids)
        for name in ('parents','due'):np.testing.assert_array_equal(getattr(sim.host,name),getattr(expected.host,name))
        assert sim.time==expected.time and sim.step_number==expected.step_number and sim.host.next_id==expected.host.next_id
        assert sim.host.divisions==expected.host.divisions and sim.host.lineage==expected.host.lineage and sim.event_audits==expected.event_audits
        for name in ('rng','fate_rng','signal_rng'):assert getattr(sim.host,name).bit_generator.state==getattr(expected.host,name).bit_generator.state
        rows.append(dict(case=job['key'],passed=True,start=.45,end=sim.time,replayed_steps=sim.step_number-before,seconds=time.perf_counter()-started,full_state_bitwise_equal=True))
        print(job['key'],'full live-furrow replay exact',flush=True);del sim,expected;torch.cuda.empty_cache()
    for f,h in watched.items():assert sha(f)==h,f
    (HERE/'restart-replay.json').write_text(json.dumps(dict(passed=True,cases=16,protocol_sha256=sha(ROOT/'protocol.json'),rows=rows,evidence_sha256=watched,script_sha256=sha(__file__)),indent=2)+'\n')


if __name__=='__main__':main()
