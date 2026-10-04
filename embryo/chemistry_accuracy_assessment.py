"""Verify and plot the completed frozen-chemistry diagnostic without changing it."""
from pathlib import Path
import tempfile

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from .chemistry_accuracy import assess
from .feedback_long import digest
from .neighbor_context import read
from .parameter_robustness import verify
from .resolution import write_json


def report(root=Path('outputs/chemistry-accuracy'), docs=Path('docs')):
    root, docs = Path(root).resolve(), Path(docs)
    p = read(root/'protocol.json'); verify(p); saved = read(root/'summary.json')
    if read(root/'status.json')['state'] != 'completed': raise ValueError('Incomplete chemistry diagnosis')
    with tempfile.TemporaryDirectory(prefix='embryo-chemistry-assessment-') as temporary:
        target = Path(temporary)
        (target/'protocol.json').symlink_to(root/'protocol.json')
        (target/'references').symlink_to(root/'references', target_is_directory=True)
        for job in p['jobs']: (target/job['key']).symlink_to(root/job['key'], target_is_directory=True)
        recomputed = assess(target)
    if recomputed != saved or not saved['passed']: raise ValueError('Frozen chemistry evidence changed or failed')
    fig, axes = plt.subplots(1, 3, figsize=(11.8, 3.8), sharey=True)
    for ax, context in zip(axes, p['contexts']):
        for backend, color, label, marker in (('native_chemistry', '#235789', 'Native CPU', 'o'),
                ('torch_cpu_chemistry', '#c67c26', 'PyTorch on CPU', 'x')):
            row = next(c for c in saved['comparisons'] if c['context'] == context['key'] and c['backend'] == backend)
            ax.loglog(p['dts'], row['errors'], color=color, marker=marker, label=label)
        ax.axhline(p['criteria']['production_log_max'], color='#ab3333', linestyle=':', label='Original 0.01 limit')
        ax.axhline(row['reference_log_error'], color='gray', linestyle='--', label='DOP853–Radau difference')
        ax.set_title(f"Frozen at elapsed {context['elapsed']:g}\n{context['duration']:g} further time units")
        ax.set_xlabel('Chemical timestep (model units)'); ax.grid(True, which='major', alpha=.2)
        ax.text(.05, .90, f"Observed order ≈ {np.mean(row['orders']):.3f}", transform=ax.transAxes)
    axes[0].set_ylabel('Maximum absolute log concentration error')
    axes[-1].legend(fontsize=8, loc='lower right')
    fig.suptitle('Chemical integration converges when measured geometry is held fixed', fontsize=12)
    fig.tight_layout(); images = docs/'images'; images.mkdir(exist_ok=True, parents=True)
    fig.savefig(images/'chemistry-accuracy.png', dpi=170); plt.close(fig)
    maximum = max(r['max_log_error'] for r in saved['results'])
    refmaximum = max(c['reference_log_error'] for c in saved['comparisons'])
    lines = ['# Frozen-geometry chemical integration diagnosis', '',
        '**Completed: all 18 runs pass the unchanged 0.01 maximum absolute log-state error limit.** Chemical errors decrease at second order under timestep halving. This isolates chemical stepping from the unresolved moving-geometry formation error; it does not validate the coupled trajectory.', '',
        '![Chemical error under timestep halving](images/chemistry-accuracy.png)', '',
        '## Why this test', '',
        'The [completed history-9 moving refinement](polarity_refinement_assessment.md) retained persistent contrast at three timesteps but failed quantitative formation agreement. Its largest adjacent-pair discrepancy was 0.037503, near elapsed 117.75. Before halving again, this assay tests whether the actual chemical routines have a large error even when their geometry inputs do not change.', '',
        '## Matched design', '',
        'Use the exact chemical states, cell volumes and conservative transport matrices saved in the fine history-9 trajectory at elapsed 0, 105 and 117.75 (physical times 150, 255 and 267.75). Evolve chemistry for another 240, 60 and 24 model time units, respectively. Within each case, freeze transport and volumes: no mechanical update, polarity evolution or dilution.', '',
        'Run the unchanged native `transport.integrate_gm` and PyTorch `gpu_backend.gm_step` functions at timesteps 0.00375, 0.001875 and 0.0009375. Parameters are beta=2, D_a=0.02 and D_b=0.55. Both chemical routines use float64; the PyTorch routine runs on CPU for this small 32-ODE check. This is not a GPU arithmetic validation. The separate spatial controls use resident PyTorch GPU tensors and custom CUDA.', '',
        'Compare every 0.15 time units against independent DOP853 and Radau solutions, with relative tolerance 1e-12, absolute tolerance 1e-14 and an analytic Jacobian. The reference disagreement must remain below 1e-8. Estimate order only if the finest production error exceeds 100 times this disagreement. The declared order interval is 1.7–2.3.', '',
        'This reuses **one developmental history (9)**; the 18 numerical cases are not 18 histories. No new zygote trajectories or identities are introduced.', '',
        '## Completed results', '',
        '| Frozen elapsed time | Duration | Coarse error | Fine error | Finer error | Estimated orders |',
        '|---|---:|---:|---:|---:|---|']
    for c in saved['comparisons']:
        if c['backend'] != 'native_chemistry': continue
        context = next(x for x in p['contexts'] if x['key'] == c['context'])
        values = ' | '.join(f'{x:.4g}' for x in c['errors'])
        lines.append(f"| {context['elapsed']:g} | {context['duration']:g} | {values} | {c['orders'][0]:.4f}, {c['orders'][1]:.4f} |")
    lines += ['',
        f'Maximum production error: **{maximum:.6g}**. Maximum DOP853–Radau discrepancy: **{refmaximum:.6g}**. Maximum native/PyTorch log-state difference: **{saved["cpu_torch_log_max"]:.6g}**. Both routines pass all contexts and all timesteps; all six refinement comparisons have resolvable order estimates.', '',
        'The error measure is the largest value of `abs(log(c_numerical / c_reference))` over both chemicals, all cells and all sampled times. For small errors, 0.01 corresponds approximately to a 1% concentration difference. Each halving reduces the observed error by about four, as expected for second-order integration on a fixed operator.', '',
        '## Interpretation and next diagnostic', '',
        'These results rule out a large isolated chemical-stepping error in the three measured fixed contexts. They leave changing transport, volume dilution, first-order splitting, spatial precision and nonlinear amplification as possible contributors. They do not separate those mechanisms and do not erase the original 0.02574 and 0.03750 failed moving comparisons.', '',
        'The [short moving precision controls](geometry_precision.md) compare the unchanged baseline with accumulation of lost phase increments, double-precision contact accumulation, and their combination. These interventions retain the equations, physical parameters and starting state. Early-window improvement alone will not establish accuracy through nonlinear formation.', '',
        '## Reproduction', '',
        'The protocol, references, 18 trajectories and summary are in `outputs/chemistry-accuracy/`. Source and input hashes were pinned before execution. This assessment independently recomputes errors from saved paths, verifies clocks and hashes, and compares the result to the saved summary without rewriting it.', '',
        '```bash',
        'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \\',
        'MPLCONFIGDIR=/tmp/embryo-mpl python -m embryo.chemistry_accuracy_assessment',
        '```', '',
        'The [verification record](chemistry_accuracy_verification.json) records this completed diagnostic separately from the previous long study. The consolidated results ledger and manuscript remain earlier snapshots.']
    (docs/'chemistry_accuracy.md').write_text('\n'.join(lines)+'\n')
    files = [root/'protocol.json', root/'summary.json', root/'status.json', root/'references/results.json']
    files += [Path(f) for f in read(root/'references/results.json')['paths_sha256']]
    files += [root/j['key']/f for j in p['jobs'] for f in ('paths.npz', 'result.json')]
    files += [Path(__file__).resolve(), docs/'chemistry_accuracy.md', images/'chemistry-accuracy.png']
    write_json(docs/'chemistry_accuracy_verification.json', dict(passed=True, independently_recomputed=True,
        numerical_cases=18, independent_histories=1, new_histories=0, maximum_log_error=maximum,
        source_sha256=p['source_sha256'], input_sha256=p['input_sha256'],
        evidence_sha256={str(f.resolve()):digest(f) for f in files}))
    return dict(passed=True, cases=18, maximum_log_error=maximum)


if __name__ == '__main__': print(report())
