"""Opposing signal pulses and additive-noise robustness of the isolated fate switch."""
import argparse
import hashlib
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from .fate_memory import neutral_fate
from .resolution import write_json, _steps


def pulse(initial,bias,duration,dt,rate=.8):
    state=np.array(initial,float,copy=True);direction=np.sign(state)
    for _ in range(_steps(duration,dt)):
        drift=lambda value:rate*(value-value**3-direction*bias)
        stage=state+dt*drift(state)
        state=.5*state+.5*(stage+dt*drift(stage))
    return state


def wilson(successes,total):
    z=1.959963984540054;p=successes/total;den=1+z*z/total
    mid=(p+z*z/(2*total))/den
    half=z*np.sqrt(p*(1-p)/total+z*z/(4*total*total))/den
    return [float(mid-half),float(mid+half)]


def noise_paths(sigmas,per_sign=256,until=60.,dt=.0075,interval=.6,seed=20260927,rate=.8):
    signs=np.r_[np.ones(per_sign),-np.ones(per_sign)]
    coarse=np.repeat(signs[None],len(sigmas),axis=0);fine=coarse.copy()
    sigma=np.array(sigmas)[:,None];rng=np.random.default_rng(seed)
    histories=[[coarse.copy()],[fine.copy()]];crossed_c=np.zeros_like(coarse,dtype=bool);crossed_f=crossed_c.copy()
    every=_steps(interval,dt)
    for i in range(1,_steps(until,dt)+1):
        dw=rng.normal(size=(2,len(signs)))*np.sqrt(dt/2)
        for increment in dw:
            fine += rate*(fine-fine**3)*(dt/2)+sigma*increment
            crossed_f |= fine*signs<0
        coarse += rate*(coarse-coarse**3)*dt+sigma*dw.sum(axis=0)
        crossed_c |= coarse*signs<0
        if not np.isfinite(fine).all() or not np.isfinite(coarse).all():raise FloatingPointError('nonfinite noise trajectory')
        if i%every==0:
            histories[0].append(coarse.copy());histories[1].append(fine.copy())
    return signs,np.array(histories[0]),np.array(histories[1]),crossed_c,crossed_f


def run(output):
    output=Path(output)
    if output.exists():raise FileExistsError('choose a fresh directory')
    p={'initial_fates':[-1.,1.],'biases':[0.,.1,.25,.375,.4,.5,.75], 'durations':[.6,3.,6.,12.,24.],
       'rate':.8,'gain':1.,'threshold':.55,'post_pulse_recovery':40.,'dts':[.0075,.00375],
       'sigmas':[0.,.05,.1,.2,.3,.4],'noise_until':60.,'noise_paths_per_initial_sign':256,'noise_seed':20260927,
       'criteria':{'pulse_fate_error_max':.05,'noise_sampled_rms_max':.05,'noise_final_flip_disagreement_max':.02,'noise_ever_cross_fraction_difference_max':.02},
       'source_sha256':{str(Path(__file__).with_name(n).resolve()):hashlib.sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in ['fate_robustness.py','fate_memory.py']},
       'scope':'Isolated assumed fate law, initially at ideal equilibria +/-1. Opposing prescribed input a=1-sign(f0)*bias (positive at all tested strengths); full reaction-diffusion chemistry and geometry absent. Noise is additive Ito forcing of fate under neutral activator. Shared Brownian increments couple fine/coarse checks; paths across noise amplitudes are paired, not independent biological replicates.'}
    output.mkdir(parents=True);write_json(output/'protocol.json',p);write_json(output/'status.json',{'state':'running'})
    try:
        rows=[];initial=np.array(p['initial_fates']);critical=2/(3*np.sqrt(3))
        for bias in p['biases']:
            for duration in p['durations']:
                endpoints=[pulse(initial,bias,duration,dt,p['rate']) for dt in p['dts']]
                references=[]
                for rtol,atol in [(1e-10,1e-13),(1e-12,1e-15)]:
                    solution=solve_ivp(lambda t,y:p['rate']*(y-y**3-np.sign(initial)*bias),(0,duration),initial,method='DOP853',rtol=rtol,atol=atol)
                    if not solution.success:raise RuntimeError(solution.message)
                    references.append(solution.y[:,-1])
                recovered=[neutral_fate(v,p['post_pulse_recovery'],p['rate']) for v in endpoints+references]
                label=lambda x:np.where(x>.55,1,np.where(x<-.55,-1,0))
                row={'bias':bias,'duration':duration,'fine_pulse_endpoint':endpoints[1].tolist(),'fine_recovered_fate':recovered[1].tolist(),
                     'switched_both_directions':bool(np.all(label(recovered[1])==-np.sign(initial))),
                     'max_error_vs_reference':float(max(np.max(abs(v-references[-1])) for v in endpoints)),
                     'max_recovered_error_vs_reference':float(max(np.max(abs(v-recovered[-1])) for v in recovered[:2])),
                     'reference_tightening_error':float(max(np.max(abs(references[0]-references[1])),np.max(abs(recovered[2]-recovered[3])))),
                     'all_labels_match_reference':bool(all(np.array_equal(label(v),label(recovered[-1])) for v in recovered[:3]))}
                rows.append(row)
        print('Pulse sweep complete',flush=True)
        signs,coarse,fine,cross_c,cross_f=noise_paths(p['sigmas'],p['noise_paths_per_initial_sign'],p['noise_until'],p['dts'][0],seed=p['noise_seed'],rate=p['rate'])
        times=np.linspace(0,p['noise_until'],len(fine));np.savez_compressed(output/'noise-trajectories.npz',times=times,initial_signs=signs,coarse=coarse,fine=fine)
        noise=[]
        for i,sigma in enumerate(p['sigmas']):
            flipped=fine[-1,i]*signs<-.55;coarse_flipped=coarse[-1,i]*signs<-.55
            retained=fine[-1,i]*signs>.55
            noise.append({'sigma':sigma,'paths':len(signs),'final_opposite_label':int(flipped.sum()),
                          'final_opposite_fraction':float(flipped.mean()),'wilson_95_interval':wilson(int(flipped.sum()),len(signs)),
                          'final_original_label':int(retained.sum()),'final_uncommitted':int((~flipped&~retained).sum()),
                          'ever_crossed_zero_fine':int(cross_f[i].sum()),
                          'max_sampled_paired_rms':float(np.sqrt(np.mean((coarse[:,i]-fine[:,i])**2,axis=-1)).max()),
                          'final_flip_disagreement':float(np.mean(flipped!=coarse_flipped)),
                          'ever_cross_fraction_difference':float(abs(cross_c[i].mean()-cross_f[i].mean())),
                          'positive_to_negative_fraction':float(flipped[signs>0].mean()),'negative_to_positive_fraction':float(flipped[signs<0].mean())})
        c=p['criteria'];checks={
            'pulse_reference_and_step_accuracy':all(r['max_error_vs_reference']<c['pulse_fate_error_max'] and r['max_recovered_error_vs_reference']<c['pulse_fate_error_max'] and r['reference_tightening_error']<1e-6 and r['all_labels_match_reference'] for r in rows),
            'no_subcritical_pulse_switches':all(not r['switched_both_directions'] for r in rows if r['bias']<=critical),
            'zero_noise_preserves_labels':noise[0]['final_original_label']==len(signs) and noise[0]['ever_crossed_zero_fine']==0,
            'noise_step_accuracy':all(r['max_sampled_paired_rms']<c['noise_sampled_rms_max'] and r['final_flip_disagreement']<=c['noise_final_flip_disagreement_max'] and r['ever_cross_fraction_difference']<=c['noise_ever_cross_fraction_difference_max'] for r in noise)}
        result={'protocol':p,'critical_constant_bias':critical,'checks':checks,'passed':all(checks.values()),'pulses':rows,'noise':noise}
        write_json(output/'comparison.json',result)
        lines=['# Fate reversal and noise robustness','',f'All checks pass: **{result["passed"]}**','',f'Constant-bias saddle-node threshold: {critical:.9f}.','',
               '| Opposing bias | Shortest tested pulse switching both directions |','|---|---:|']
        for bias in p['biases']:
            values=[r['duration'] for r in rows if r['bias']==bias and r['switched_both_directions']]
            lines.append(f'| {bias:g} | {min(values) if values else "None through 24"} |')
        lines+=['','| Noise amplitude | Opposite label at t=60 / 512 | Approx. 95% interval | Ever crossed zero | Uncommitted at t=60 |','|---|---:|---:|---:|---:|']
        for r in noise:lines.append(f'| {r["sigma"]:g} | {r["final_opposite_label"]} | {r["wilson_95_interval"][0]:.3f}–{r["wilson_95_interval"][1]:.3f} | {r["ever_crossed_zero_fine"]} | {r["final_uncommitted"]} |')
        lines+=['',*[f'- {k}: {v}' for k,v in checks.items()], '',
                'Pulse thresholds are conditional on an initial equilibrium and finite tested durations, not universal commitment thresholds. Noise crossings can reverse again; final opposite labels and ever crossing zero are different outcomes. Zero observed switches does not prove zero switching probability. Intervals quantify Monte Carlo sampling under this SDE, not biological uncertainty. This assay does not test the full embryo, chemical network, or robustness across developmental histories.']
        (output/'RESULTS.md').write_text('\n'.join(lines)+'\n')
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,axs=plt.subplots(1,2,figsize=(10,4),layout='constrained')
        matrix=np.array([[next(r['fine_recovered_fate'][1] for r in rows if r['bias']==b and r['duration']==d) for d in p['durations']] for b in p['biases']])
        im=axs[0].imshow(matrix,vmin=-1,vmax=1,cmap='coolwarm',aspect='auto',origin='lower')
        axs[0].set(xticks=range(len(p['durations'])),xticklabels=p['durations'],yticks=range(len(p['biases'])),yticklabels=p['biases'],xlabel='Opposing pulse duration',ylabel='Bias strength',title='Recovered fate (initial +1)');fig.colorbar(im,ax=axs[0])
        y=np.array([r['final_opposite_fraction'] for r in noise]);ci=np.array([r['wilson_95_interval'] for r in noise]).T
        axs[1].errorbar(p['sigmas'],y,yerr=np.maximum(np.array([y-ci[0],ci[1]-y]),0),fmt='o-',capsize=3)
        axs[1].set(xlabel='Additive fate noise amplitude',ylabel='Fraction with opposite label at t=60',title='512 independent noise paths per amplitude');axs[1].grid(alpha=.2)
        fig.savefig(output/'comparison.png',dpi=160);plt.close(fig)
        write_json(output/'status.json',{'state':'completed','passed':result['passed']});return result
    except Exception as error:
        write_json(output/'status.json',{'state':'failed','error':str(error)});raise


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();run(args.output)
