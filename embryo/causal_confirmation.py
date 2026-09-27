"""Independent finer-time confirmation of the fixed-geometry causal screen."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

from .causal_signaling import simulate
from .model import Simulation
from .resolution import write_json


def run(source,output):
    source,output=Path(source).resolve(),Path(output)
    if output.exists():raise FileExistsError('choose a fresh output directory')
    old=json.loads((source/'protocol.json').read_text());checkpoint=Path(old['checkpoint'])
    if hashlib.sha256(checkpoint.read_bytes()).hexdigest()!=old['checkpoint_sha256']:raise ValueError('source geometry changed')
    protocol={'source':str(source),'source_report_sha256':hashlib.sha256((source/'comparison.json').read_bytes()).hexdigest(),
              'checkpoint_sha256':old['checkpoint_sha256'],'arms':['full','no_self_activation','no_transport','equal_diffusion'],
              'dts':[.00125,.000625],'signal_rms_max':.01,'fate_absolute_max':.05,
              'seeds':old['seeds'],'until':old['until'],'interval':old['interval'],
              'scope':'New finer-time comparisons after original fate-timing check failed. Original result unchanged; same 20 seeds, geometry, perturbations, criteria and reaction equations.'}
    output.mkdir(parents=True);write_json(output/'protocol.json',protocol);write_json(output/'status.json',{'state':'running'})
    sim=Simulation.restore(checkpoint);graph=sim.signaling_graph();rows=[]
    try:
        for arm in protocol['arms']:
            reports=[]
            for dt in protocol['dts']:
                report=simulate(graph,sim.config,old['seeds'],arm,dt,old['until'],old['interval'],old['perturbation_rms'],old['regulator_ceiling'])
                write_json(output/f'{arm}-{dt:g}.json',report);reports.append(report)
            coarse,fine=reports
            signal=max(float(np.max(np.sqrt(np.mean((np.array(a[k])-b[k])**2,axis=1)))) for a,b in zip(coarse['history'],fine['history']) for k in ('activator','inhibitor'))
            fate=max(float(np.max(abs(np.array(a['fate'])-b['fate']))) for a,b in zip(coarse['history'],fine['history']))
            row={'arm':arm,'signal_rms_difference':signal,'fate_max_difference':fate,
                 'all_complete':all(x['completed'] for r in reports for x in r['outcomes']),
                 'final_counts_match':all(a['final_fate_a']==b['final_fate_a'] and a['final_fate_b']==b['final_fate_b'] for a,b in zip(coarse['outcomes'],fine['outcomes'])),
                 'fine_persistent_signals':sum(x['persistent_contrast_above_0_1'] for x in fine['outcomes']),
                 'fine_both_fates':sum(x['both_fates_at_end'] for x in fine['outcomes'])}
            rows.append(row);print(row,flush=True)
            write_json(output/'progress.json',rows)
        unchanged=hashlib.sha256((source/'comparison.json').read_bytes()).hexdigest()==protocol['source_report_sha256'] and hashlib.sha256(checkpoint.read_bytes()).hexdigest()==protocol['checkpoint_sha256']
        passed=unchanged and all(r['all_complete'] and r['final_counts_match'] and r['signal_rms_difference']<protocol['signal_rms_max'] and r['fate_max_difference']<protocol['fate_absolute_max'] for r in rows)
        result={'protocol':protocol,'rows':rows,'source_unchanged':unchanged,'passed':passed}
        write_json(output/'comparison.json',result)
        lines=['# Causal screen finer-time confirmation','',f'All confirmation checks pass: **{passed}**','',protocol['scope'],'',
               '| Arm | Signal RMS discrepancy | Maximum fate discrepancy | Both fates / 20 | Persistent signals / 20 |','|---|---:|---:|---:|---:|']
        lines += [f'| {r["arm"]} | {r["signal_rms_difference"]:.6g} | {r["fate_max_difference"]:.6g} | {r["fine_both_fates"]} | {r["fine_persistent_signals"]} |' for r in rows]
        (output/'RESULTS.md').write_text('\n'.join(lines)+'\n');write_json(output/'status.json',{'state':'completed','passed':passed})
        return result
    except Exception as error:
        write_json(output/'status.json',{'state':'failed','error':str(error)});raise


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--source',type=Path,default=Path('outputs/causal-signaling'));parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();run(args.source,args.output)
