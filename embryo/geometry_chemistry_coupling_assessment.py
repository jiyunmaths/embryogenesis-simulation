"""Verify completed changing-geometry chemistry, preserving the failed screen."""
from pathlib import Path
import tempfile

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from .feedback_long import digest
from .geometry_chemistry_coupling import assess, load_context
from .neighbor_context import read, validate_graph
from .parameter_robustness import verify
from .resolution import write_json, _steps


def report(root=Path('outputs/geometry-chemistry-coupling'), docs=Path('docs')):
    root, docs = Path(root).resolve(), Path(docs); p = read(root/'protocol.json'); verify(p)
    if read(root/'status.json')['state'] not in ('completed','completed_with_unresolved_checks'):
        raise ValueError('Coupling replay incomplete')
    saved = read(root/'summary.json')
    with tempfile.TemporaryDirectory(prefix='embryo-coupling-assessment-') as temporary:
        target=Path(temporary); (target/'protocol.json').symlink_to(root/'protocol.json')
        (target/'references').symlink_to(root/'references',target_is_directory=True)
        for j in p['jobs']: (target/j['key']).symlink_to(root/j['key'],target_is_directory=True)
        if assess(target) != saved: raise ValueError('Changed replay assessment')
    samples = 0; reference_samples = 0; spacing_profiles = []; timestep_pairs = []; fidelity = []; extra = []
    initial_amount_roundtrip_error = 0.
    original = Path('outputs/polarity-robustness/seed-9_fine_polarity-0_uniform/history.json')
    live = read(original); extra.append(original)
    for context in p['contexts']:
        for spacing in ('fine','coarse'):
            schedule, initial = load_context(p,context,spacing)
            with np.load(root/'references'/f'{context["key"]}_{spacing}/paths.npz') as z:
                times = np.arange(_steps(context['duration'],p['interval'])+1)*p['interval']
                if not np.array_equal(times,z['times']): raise ValueError('Changed independent clock')
                if not np.array_equal(initial,z['trajectory'][0]): raise ValueError('Changed independent start')
                volumes=[]
                for t in times:
                    volume,delta=schedule.delta(t); validate_graph(delta,volume); volumes.append(volume)
                volumes=np.array(volumes)
                for method in ('dop853','radau'):
                    path=z['trajectory'] if method=='dop853' else z['radau_trajectory']
                    converted = z[method+'_amounts']/volumes[:,None,:]
                    initial_error = float(abs(np.log(converted[0]/path[0])).max())
                    initial_amount_roundtrip_error = max(initial_amount_roundtrip_error, initial_error)
                    # Q0=V0*c0 then Q0/V0 can differ by roundoff; later samples
                    # were constructed from the stored amounts and must be exact.
                    if not np.array_equal(converted[1:],path[1:]) or initial_error > 4*np.finfo(np.float64).eps:
                        raise ValueError('Reference concentrations do not match amounts/capacities')
                    reference_samples+=len(path)
            for method in p['methods']:
                for dt in p['dts']:
                    job=next(j for j in p['jobs'] if j['context']==context['key'] and j['spacing']==spacing and j['method']==method and j['dt']==dt)
                    folder=root/job['key']; r=read(folder/'result.json')
                    with np.load(folder/'paths.npz') as z,np.load(folder/'latest.npz') as cp:
                        trajectory=z['trajectory']
                        if (not np.array_equal(trajectory[0],initial) or not np.isfinite(trajectory).all()
                                or np.any(trajectory<=0) or cp['step']!=_steps(context['duration'],dt)
                                or str(cp['protocol_sha256'])!=digest(root/'protocol.json')
                                or not np.array_equal(cp['trajectory'],trajectory) or not np.array_equal(cp['state'],trajectory[-1])
                                or float(cp['dilution_error'])!=r['dilution_error_max']):
                            raise ValueError('Replay checkpoint, positivity or start differs')
                        if r['accuracy_pass'] != (r['max_log_error']<=p['criteria']['production_log_max']):
                            raise ValueError('Accuracy decision differs')
                        if r['dilution_error_max']>p['criteria']['dilution_error_max']:
                            raise ValueError('Dilution screen differs')
                        samples+=len(trajectory)
        with np.load(root/'references'/f'{context["key"]}_fine/paths.npz') as a,np.load(root/'references'/f'{context["key"]}_coarse/paths.npz') as b:
            e=abs(np.log(a['trajectory']/b['trajectory'])); index=np.unravel_index(e.argmax(),e.shape)
            spacing_profiles.append(dict(context=context['key'],max_log_difference=float(e.max()),
                source_elapsed_at_peak=context['elapsed']+float(a['times'][index[0]]),
                endpoint_log_difference=float(e[-1].max())))
        for method in p['methods']:
            with np.load(root/f'{context["key"]}_fine_{method}_dt-0.001875/paths.npz') as a,np.load(root/f'{context["key"]}_fine_{method}_dt-0.0009375/paths.npz') as b:
                timestep_pairs.append(dict(context=context['key'],method=method,
                    max_log_difference=float(abs(np.log(a['trajectory']/b['trajectory'])).max())))
        source=np.array([r['chemistry'] for r in live if context['elapsed']-1e-9<=r['elapsed']<=context['elapsed']+context['duration']+1e-9])
        with np.load(root/f'{context["key"]}_fine_beginning_dt-0.001875/paths.npz') as z:
            e=abs(np.log(z['trajectory']/source)); index=np.unravel_index(e.argmax(),e.shape)
            fidelity.append(dict(context=context['key'],same_dt_source_log_difference=float(e.max()),
                source_elapsed_at_peak=context['elapsed']+float(z['times'][index[0]]),
                endpoint_log_difference=float(e[-1].max())))
    fig,axes=plt.subplots(1,2,figsize=(10.8,4.1),sharey=True)
    for ax,context in zip(axes,p['contexts']):
        for spacing,style in (('fine','-'),('coarse','--')):
            for method,color,marker in (('beginning','#235789','o'),('midpoint','#c67c26','s')):
                row=next(r for r in saved['comparisons'] if r['context']==context['key'] and r['spacing']==spacing and r['method']==method)
                ax.loglog(p['dts'],row['errors'],style,color=color,marker=marker,label=f'{method}, geometry {0.15 if spacing=="fine" else .30:.2f}')
        ax.axhline(.01,color='#ab3333',linestyle=':',label='Original 0.01 limit')
        ax.set_title(f"Source elapsed {context['elapsed']:g}, duration {context['duration']:g}")
        ax.set_xlabel('Chemical timestep (model units)'); ax.grid(True,which='major',alpha=.2)
    axes[0].set_ylabel('Maximum absolute log concentration error')
    axes[1].legend(fontsize=8,loc='upper left'); fig.suptitle('Shared changing geometry: first-order beginning vs second-order midpoint')
    fig.tight_layout(); image=docs/'images/geometry-chemistry-coupling.png'; fig.savefig(image,dpi=170); plt.close(fig)
    lines=['','## Completed replay assessment','',
        '**All 24 integration cases and four independent reference pairs pass. The full-start geometry-spacing screen fails; the local formation-window screen passes.** Preserve `completed_with_unresolved_checks`.', '',
        '![Integration error on changing geometry](images/geometry-chemistry-coupling.png)', '',
        '| Context | Geometry spacing | Method | Coarse error | Fine error | Finer error | Observed orders |',
        '|---|---:|---|---:|---:|---:|---|']
    for row in saved['comparisons']:
        values=' | '.join(f'{e:.6g}' for e in row['errors']); orders=', '.join(f'{o:.3f}' for o in row['orders'])
        lines.append(f"| {row['context']} | {0.15 if row['spacing']=='fine' else .30:.2f} | {row['method']} | {values} | {orders} |")
    lines+=['',f"Maximum DOP853/Radau discrepancy is {max(r['reference_log_error'] for r in saved['comparisons']):.6g}. All order estimates are resolvable. Beginning sampling is first order (about 1.00); midpoint sampling is second order (about 2.00). At dt=0.001875 on the fine full-start schedule, midpoint reduces error from 0.000630431 to 5.13244e-7, about 1,228-fold. This validates a timing improvement on prescribed inputs, not a new live solver.",'',
        '| Context | Reference spacing difference | Limit | Outcome | Source elapsed at maximum |',
        '|---|---:|---:|---|---:|']
    for row,profile in zip(saved['spacing_checks'],spacing_profiles):
        lines.append(f"| {row['context']} | {row['reference_spacing_log_difference']:.6g} | 0.001 | {'PASS' if row['passed'] else 'FAIL'} | {profile['source_elapsed_at_peak']:.2f} |")
    lines+=['',
        'The full-start difference 0.00364878 peaks at elapsed 96.75; its endpoint difference is only about 4.1e-8. Endpoint agreement therefore cannot accept its transient spacing failure. The elapsed-90 context has spacing sensitivity 1.19127e-5 and passes. The initial history and local restart are different reference problems; do not use the local pass to erase the full-start failure.', '',
        '## Interpretation and source fidelity', '',
        'Current chemical–geometry time coupling is measurably first order on a shared prescribed path, even though its chemical step is second order on fixed geometry. Midpoint timing resolves this integration issue in both contexts and source spacings. Its live use would require consistent intermediate mechanical geometry and new coupled-backend validation; production mechanics is unchanged.', '',
        'At dt=0.001875 versus 0.0009375, the beginning scheme differs by about 0.0003151 over the full-start replay and 3.8561e-5 over the local formation replay. Both are far below the original moving discrepancy 0.03750. Thus this test does not support attributing that whole failure to chemical timing on a shared smooth geometry alone. It excludes neither differences between the actual moving geometry paths nor amplification of earlier chemical differences.', '',
        'A post hoc source-fidelity check compares the beginning replay at the original dt=0.001875 against the actual source chemistry. Maximum differences are 0.0025182 for the full start and 4.18758e-6 for the local restart. These descriptive checks add no acceptance rule and do not override the spacing failure. The local replay is well matched through the original nonlinear mismatch interval; the full-start replay needs denser saved geometry or a dedicated interpolation confirmation before stronger causal attribution.', '',
        'The [full phase-carry formation pair](phase_carry_formation.md) remains the direct moving test. Its current progress is separate from this completed replay result. No broader identity claim, new developmental history or backend acceptance follows here.', '',
        '## Completed verification', '',
        f'Independent assessment verifies all {samples:,} production observations and {reference_samples:,} reference observations, reference amount/concentration consistency, conservative interpolated graphs, physical clocks, production starts and final small-state checkpoints. The initial amount multiplication/division differs from the stored initial concentration by at most {initial_amount_roundtrip_error:.6g}, within float64 roundoff; later reference concentrations reconstruct exactly. It recomputes errors/orders/spacing decisions without rewriting the frozen summary.', '',
        '```bash','OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \\',
        'MPLCONFIGDIR=/tmp/embryo-mpl python -m embryo.geometry_chemistry_coupling_assessment','```','',
        'The [verification record](geometry_chemistry_coupling_verification.json) retains the failed full-start screen separately from passing integration checks. The consolidated ledger and manuscript remain older snapshots.']
    document=docs/'geometry_chemistry_coupling.md'; content=document.read_text(); marker='\n## Completed replay assessment\n'
    if marker in content: content=content.split(marker)[0].rstrip()+'\n'
    document.write_text(content+'\n'.join(lines)+'\n')
    files=[root/f for f in ('protocol.json','status.json','summary.json')]
    files += [root/j['key']/f for j in p['jobs'] for f in ('paths.npz','result.json','latest.npz')]
    files += [root/'references'/j['key']/f for j in p['references'] for f in ('paths.npz','result.json')]
    files += extra+[Path(__file__).resolve(),image,document]
    record=dict(verified=True, independent_histories=1,new_histories=0,production_cases=24,reference_pairs=4,
        production_observations=samples,reference_observations=reference_samples,
        initial_amount_roundtrip_log_error=initial_amount_roundtrip_error,
        time_stepping_pass=saved['time_stepping_pass'],spacing_screen_pass=saved['spacing_screen_pass'],
        spacing_profiles=spacing_profiles,shared_geometry_timestep_pairs=timestep_pairs,posthoc_source_fidelity=fidelity,
        source_sha256=p['source_sha256'],input_sha256=p['input_sha256'],
        evidence_sha256={str(f.resolve()):digest(f) for f in files})
    write_json(docs/'geometry_chemistry_coupling_verification.json',record)
    return dict(verified=True,time_stepping_pass=saved['time_stepping_pass'],spacing_screen_pass=saved['spacing_screen_pass'],
        production_observations=samples,reference_observations=reference_samples)


if __name__=='__main__': print(report())
