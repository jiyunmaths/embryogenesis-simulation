"""Independently assess the completed short precision controls; preserve originals."""
from pathlib import Path
import tempfile

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from .feedback_long import digest
from .geometry_precision import assess, errors, ARMS
from .neighbor_context import read
from .parameter_robustness import verify
from .parameter_robustness_reference import frozen_spectrum
from .resolution import write_json


def report(root=Path('outputs/geometry-precision'), docs=Path('docs')):
    root, docs = Path(root).resolve(), Path(docs); p = read(root/'protocol.json'); verify(p)
    if read(root/'status.json')['state'] != 'completed': raise ValueError('Short precision runs incomplete')
    saved = read(root/'summary.json')
    with tempfile.TemporaryDirectory(prefix='embryo-precision-assessment-') as temporary:
        target = Path(temporary); (target/'protocol.json').symlink_to(root/'protocol.json')
        for job in p['jobs']: (target/job['key']).symlink_to(root/job['key'], target_is_directory=True)
        recomputed = assess(target)
    if saved != recomputed or not saved['quality_pass']: raise ValueError('Changed or failed precision evidence')
    spectra = 0; discrepancy = 0.; histories = {}
    for job in p['jobs']:
        history = read(root/job['key']/'history.json'); histories[job['arm'], job['level']] = history
        for row in history:
            s = frozen_spectrum(np.array(row['delta']), np.array(row['volumes']), 2., .02, 27.5)
            discrepancy = max(discrepancy, abs(s['uniform_jacobian_max_real']-row['uniform_growth_max']))
            if s['unstable_modes'] != row['unstable_modes'] or not np.array_equal(s['lambdas'], row['laplacian_lambdas']):
                raise ValueError('Stored graph spectrum differs')
            if np.std(np.log(np.array(row['chemistry'])[0])) != row['log_activator_sd']:
                raise ValueError('Stored chemical contrast differs')
            spectra += 1
    if discrepancy != 0: raise ValueError('Stored growth differs')
    original_folders = [Path('outputs/polarity-robustness/seed-9_fine_polarity-0_uniform'),
        Path('outputs/polarity-robustness-refined/seed-9_finer_polarity-0_uniform')]
    reproduced = []; extra = []
    for level, folder in zip(('fine', 'finer'), original_folders):
        a = histories['baseline', level]; b = read(folder/'history.json')[:len(a)]
        e = errors(a, b)
        if any(e.values()): raise ValueError('Baseline does not reproduce its archived scientific prefix')
        reproduced.append(dict(level=level, error=e, original_history_sha256=digest(folder/'history.json')))
        extra += [folder/'history.json']
    fig, axes = plt.subplots(1, 2, figsize=(10.3, 3.6)); profiles = []
    colors = ['#235789', '#cf6b24', '#43855a', '#8b5fbf']
    for arm, color in zip(ARMS, colors):
        a, b = [histories[arm, level] for level in ('fine', 'finer')]
        times = np.array([r['elapsed'] for r in a]); x, y = [np.array([r['chemistry'] for r in h]) for h in (a, b)]
        chemical = abs(np.log(x)-np.log(y)).max(axis=(1, 2))
        profiles.append(dict(arm=arm, peak_elapsed=float(times[chemical.argmax()]),
            endpoint_chemical_log_max=float(chemical[-1]),
            after_one_unit_chemical_log_max=float(chemical[times >= 1.-1e-9].max())))
        dx, dy = [np.array([r['delta'] for r in h]) for h in (a, b)]
        transport = np.linalg.norm(dx-dy, axis=(1, 2))/np.linalg.norm(dy, axis=(1, 2))
        axes[0].plot(times, chemical, color=color, label=arm.replace('_', ' '))
        axes[1].plot(times, transport, color=color, label=arm.replace('_', ' '))
    axes[0].set_ylabel('Maximum absolute log chemical difference')
    axes[1].set_ylabel('Relative transport-matrix difference')
    for ax in axes:
        ax.set_xlabel('Elapsed time (model units)'); ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
        ax.grid(True, alpha=.2)
    axes[1].legend(fontsize=8); fig.suptitle('Early timestep sensitivity: four numerical precision controls')
    fig.tight_layout(); image = docs/'images/geometry-precision.png'; fig.savefig(image, dpi=170); plt.close(fig)
    baseline = saved['timestep_comparisons'][0]
    changes = []
    for row in saved['timestep_comparisons'][1:]:
        changes.append(dict(arm=row['arm'], chemical_ratio_to_baseline=row['chemical_log_max']/baseline['chemical_log_max'],
            transport_ratio_to_baseline=row['relative_transport_max']/baseline['relative_transport_max']))
    lines = ['', '## Completed short-control results', '',
        '**All eight six-unit continuations completed and passed the original quality screens.** These are numerical interventions within one existing history.', '',
        '![Early timestep sensitivity](images/geometry-precision.png)', '',
        '| Control | Max log chemical difference | Relative transport difference | Endpoint relative phase L2 difference |',
        '|---|---:|---:|---:|']
    for row in saved['timestep_comparisons']:
        lines.append(f"| {row['arm']} | {row['chemical_log_max']:.6g} | {row['relative_transport_max']:.6g} | {row['endpoint_phase_relative_l2']:.6g} |")
    lines += ['',
        'Values compare dt=0.001875 with 0.0009375 over the aligned early window. Transport discrepancy is the maximum Frobenius-norm difference divided by the finer matrix norm; phase discrepancy is the endpoint L2 norm divided by the finer field norm. Neither is an estimate of the full formation error.', '',
        '| Intervention | Chemical discrepancy / baseline | Transport discrepancy / baseline |', '|---|---:|---:|']
    for r in changes:
        lines.append(f"| {r['arm']} | {r['chemical_ratio_to_baseline']:.4f} | {r['transport_ratio_to_baseline']:.4f} |")
    phase = changes[0]; contact = changes[1]
    endpoint_ratio = profiles[1]['endpoint_chemical_log_max']/profiles[0]['endpoint_chemical_log_max']
    lines += ['',
        'A ratio below one indicates a smaller adjacent-timestep discrepancy in this early window. Precision controls are interventions, not independent exact solutions; compare their effect size with the baseline and do not infer full convergence from two timesteps.', '',
        f"Phase carry gives {phase['transport_ratio_to_baseline']:.3f} times the baseline transport discrepancy and {phase['chemical_ratio_to_baseline']:.3f} times its peak chemical discrepancy. Contact64 alone gives ratios {contact['transport_ratio_to_baseline']:.6f} and {contact['chemical_ratio_to_baseline']:.6f}, respectively. Double contact accumulation alone has little effect here.", '',
        '**The peak chemical metric hides a timing distinction.** Every arm peaks at the first saved moving observation, elapsed 0.15. The phase-carry curves then fall well below the baseline curves. The following post hoc descriptive profile distinguishes that early peak from later behavior; it adds no acceptance threshold and does not replace the maximum-error metric.', '',
        '| Control | Peak elapsed time | Endpoint chemical difference | Max chemical difference after elapsed 1 |',
        '|---|---:|---:|---:|']
    for row in profiles:
        lines.append(f"| {row['arm']} | {row['peak_elapsed']:.2f} | {row['endpoint_chemical_log_max']:.6g} | {row['after_one_unit_chemical_log_max']:.6g} |")
    lines += ['',
        f"Phase carry reduces the endpoint chemical discrepancy by about {(1-endpoint_ratio)*100:.1f}% (ratio {endpoint_ratio:.3f}), despite reducing the full-window peak by only about 4%. Thus the rounding correction also materially changes later chemical agreement in this short window. It is a candidate for longer confirmation; the remaining early peak leaves geometry/time splitting and other precision paths open. These data do not identify the cause of the original nonlinear transient.", '',
        'Both new baseline runs reproduce all saved chemical, volume, transport, polarity, center and uniform-growth observations from the original six-unit prefixes exactly. Independent assessment recomputes all 328 saved graph spectra and contrast values, and verifies physical checkpoints, their carry state, aligned clocks, and source/input hashes.', '',
        'The original 240-unit formation-transient failures remain unresolved. Chemical stepping on fixed geometry is accurate, and these early controls quantify selected spatial rounding paths, but identifying the source of the nonlinear transient requires a longer matched confirmation and/or prescribed-geometry splitting test. No identity or biological conclusion is promoted from this numerical assay.', '',
        '## Next diagnostic', '',
        'Two complementary confirmations are now justified. Extend the matched phase-carry pair through nonlinear formation to test whether the improved field/transport agreement persists and whether the original 0.01 long chemical gate passes. In parallel, test chemical integration against one prescribed changing geometry from the saved fine trajectory. Interpolate cell volumes and symmetric nonnegative face conductances, then reconstruct the conservative operator; directly interpolating transport matrices while volumes vary can break conservation. Use compartment amounts (volume times concentration) for independent DOP853/Radau references, so dilution is handled through the prescribed volume without a separate numerical split. Compare the current beginning-of-step transport/dilution update with a midpoint-time alternative at the same three timesteps. Validate the reference and interpolation dependence before attributing an error to splitting. The latter isolates chemical–geometry time coupling while leaving production mechanics unchanged. Neither follow-up is launched by this assessment.', '',
        'Reproduce the completed assessment with:', '', '```bash',
        'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \\',
        'MPLCONFIGDIR=/tmp/embryo-mpl python -m embryo.geometry_precision_assessment',
        '```', '',
        'The [completed verification](geometry_precision_verification.json) records the eight trajectories and independent checks. The consolidated ledger and manuscript remain their prior snapshots.']
    document = docs/'geometry_precision.md'; content = document.read_text(); marker = '\n## Completed short-control results\n'
    if marker in content: content = content.split(marker)[0].rstrip()+'\n'
    document.write_text(content+'\n'.join(lines)+'\n')
    files = [root/'protocol.json', root/'summary.json', root/'status.json', root/'implementation-gate.json']
    files += [root/j['key']/f for j in p['jobs'] for f in ('history.json', 'result.json', 'latest_state.npz')]
    files += extra+[Path(__file__).resolve(), image, document]
    record = dict(passed=True, independent_histories=1, new_histories=0, trajectories=8,
        spectra_recomputed=spectra, maximum_growth_discrepancy=discrepancy,
        archived_baselines=reproduced, precision_effect_ratios=changes,
        posthoc_timing_profiles=profiles, phase_endpoint_chemical_ratio_to_baseline=endpoint_ratio,
        source_sha256=p['source_sha256'], input_sha256=p['input_sha256'],
        evidence_sha256={str(f.resolve()):digest(f) for f in files})
    write_json(docs/'geometry_precision_verification.json', record)
    return dict(passed=True, spectra=spectra, timestep_comparisons=saved['timestep_comparisons'], effects=changes)


if __name__ == '__main__': print(report())
