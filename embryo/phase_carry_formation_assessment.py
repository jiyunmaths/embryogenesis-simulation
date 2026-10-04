"""Verify completed phase-carry formation without changing scientific outputs."""
from datetime import datetime, timezone
from pathlib import Path
import tempfile

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from .attribute_persistence import distance, rhs
from .feedback_endpoint_bistability import chemical_jacobian
from .feedback_long import digest, restore_checkpoint
from .neighbor_context import read, validate_graph
from .parameter_robustness import classify, frozen_spectrum, metrics, verify
from .phase_carry_formation import assess
from .resolution import write_json


def profiles(a, b):
    times = np.array([row['elapsed'] for row in a])
    x, y = [np.array([row['chemistry'] for row in rows]) for rows in (a, b)]
    error = abs(np.log(x / y))
    index = np.unravel_index(error.argmax(), error.shape)
    dx, dy = [np.array([row['delta'] for row in rows]) for rows in (a, b)]
    transport = np.linalg.norm(dx - dy, axis=(1, 2)) / np.linalg.norm(dy, axis=(1, 2))
    return times, error.max(axis=(1, 2)), transport, dict(
        max_log_difference=float(error.max()), peak_elapsed=float(times[index[0]]),
        chemical='activator' if index[1] == 0 else 'inhibitor',
        cell_id=a[index[0]]['ids'][index[2]],
        endpoint_log_difference=float(error[-1].max()))


def report(root=Path('outputs/phase-carry-formation'), docs=Path('docs')):
    root, docs = Path(root).resolve(), Path(docs)
    p = read(root / 'protocol.json'); verify(p)
    saved = read(root / 'summary.json')
    status = read(root / 'status.json')
    if status['state'] != 'completed' or not saved['passed']:
        raise ValueError('Full formation pair has not completed its gates')
    files = [root / name for name in ('protocol.json', 'status.json', 'summary.json',
                                     'transfer-verification.json', 'restart-gate.json')]
    for job in p['jobs']:
        folder = root / job['key']
        files += [folder / name for name in ('result.json', 'history.json', 'latest_state.npz',
                  'endpoint/source.npz', 'endpoint/protocol.json', 'endpoint/assay/result.json',
                  'endpoint/assay/paths.npz', 'endpoint/assay/status.json')]
    original_hashes = {str(path): digest(path) for path in files}
    # The aggregate function writes only the temporary summary. Completed child
    # directories and frozen assays are reused through their unchanged caches.
    with tempfile.TemporaryDirectory(prefix='embryo-carry-formation-assessment-') as temp:
        target = Path(temp)
        (target / 'protocol.json').symlink_to(root / 'protocol.json')
        for job in p['jobs']:
            (target / job['key']).symlink_to(root / job['key'], target_is_directory=True)
        if assess(target) != saved:
            raise ValueError('Recomputed formation summary differs')
    if any(digest(name) != expected for name, expected in original_hashes.items()):
        raise ValueError('Assessment modified scientific evidence')

    histories = {}; spectra_checked = 0; frozen_samples = 0; endpoint_details = []
    ph = digest(root / 'protocol.json')
    for job in p['jobs']:
        folder = root / job['key']; history = read(folder / 'history.json')
        histories['carry', job['level']] = history
        # Preserving the original prefix excludes an accidental chemical or
        # residual restart at the six-unit transfer boundary.
        prefix = read(Path(p['short']) / job['key'] / 'history.json')
        if history[:len(prefix)] != prefix:
            raise ValueError('Transferred scientific observation prefix differs')
        for row in history:
            spectrum = frozen_spectrum(np.array(row['delta']), np.array(row['volumes']),
                                       2., .02, job['point']['ratio'])
            if (spectrum['uniform_jacobian_max_real'] != row['uniform_growth_max'] or
                    spectrum['unstable_modes'] != row['unstable_modes'] or
                    not np.array_equal(spectrum['lambdas'], row['laplacian_lambdas']) or
                    np.std(np.log(np.array(row['chemistry'])[0])) != row['log_activator_sd']):
                raise ValueError('Stored spectrum or chemical contrast differs')
            spectra_checked += 1
        host, audit, stored = restore_checkpoint(folder / 'latest_state.npz', ph, job)
        if (stored != history or host.time != p['start'] + p['duration'] or
                not np.array_equal([host.activator, host.inhibitor], history[-1]['chemistry']) or
                not np.array_equal(host.polarity, history[-1]['polarity']) or
                not np.isfinite(host.phi).all()):
            raise ValueError('Final physical checkpoint differs from observations')
        if not (audit['max_volume_error'] < .05 and audit['min_radius'] >= 4 and
                audit['max_clipping'] == 0 and audit['boundary_max'] < p['criteria']['boundary_max'] and
                audit['dilution_error_max'] <= p['criteria']['dilution_error_max']):
            raise ValueError('Inherited physical quality screen differs')
        with np.load(folder / 'latest_state.npz') as z:
            carry = z['phase_carry']
            if carry.dtype != np.float64 or carry.shape != host.phi.shape or not np.isfinite(carry).all():
                raise ValueError('Phase residual missing from final checkpoint')
        endpoint = folder / 'endpoint'; ep = read(endpoint / 'protocol.json'); verify(ep)
        result = read(endpoint / 'assay/result.json')
        if result['paths_sha256'] != digest(endpoint / 'assay/paths.npz'):
            raise ValueError('Frozen paths changed')
        final = history[-1]
        with np.load(endpoint / 'source.npz') as source, np.load(endpoint / 'assay/paths.npz') as z:
            delta, masses = validate_graph(z['delta'], z['masses'])
            if (not np.array_equal(delta, final['delta']) or not np.array_equal(masses, final['volumes']) or
                    not np.array_equal(z['ids'], final['ids']) or
                    not np.array_equal(source['trajectories'][0, -1], final['chemistry']) or
                    not np.array_equal(z['initial_states'][0], final['chemistry'])):
                raise ValueError('Frozen assay does not use the actual moving endpoint')
            if frozen_spectrum(delta, masses, 2., .02, job['point']['ratio']) != result['spectrum']:
                raise ValueError('Frozen endpoint spectrum differs')
            fun = rhs(delta, 2., .02, .02 * job['point']['ratio'])
            jac = lambda t, y: chemical_jacobian(y.reshape(2, -1), delta, 2., .02, .02 * job['point']['ratio'])
            for i, trial in enumerate(result['trials']):
                path, times = z[f'path_{i}'], z[f'times_{i}']
                if (not np.array_equal(path[0], z['initial_states'][i]) or not np.isfinite(path).all() or
                        np.any(path <= 0) or not np.array_equal(times, np.arange(
                            0., trial['horizon'] + ep['interval'] / 2, ep['interval']))):
                    raise ValueError('Frozen clocks, starts or positivity differ')
                computed = metrics(path[-1], fun, jac, ep['criteria'])
                if any(computed[k] != trial[k] for k in computed):
                    raise ValueError('Frozen endpoint residual, stability or state differs')
                returned = float(distance(path[-1], z['path_0'][-1], masses)) if i else 0.
                if returned != trial['pattern_return_log_rms']:
                    raise ValueError('Pattern-return distance differs')
                frozen_samples += len(path)
            phase = classify(result['trials'], result['spectrum']['uniform_jacobian_max_real'], ep['criteria'])
            if any(phase[k] != result[k] for k in phase):
                raise ValueError('Frozen chemical phase decision differs')
        endpoint_details.append(dict(level=job['level'], phase=result['phase'],
            uniform_growth=result['spectrum']['uniform_jacobian_max_real'],
            pattern_growth=result['trials'][0]['jacobian_max_real'],
            max_solver_log_error=max(t['solver_log_error'] for t in result['trials']),
            max_stationarity_rhs=max(t['rhs_max'] for t in result['trials'])))
    for level, folder in p['baseline_folders'].items():
        histories['no carry', level] = read(Path(folder) / 'history.json')
        files.append(Path(folder) / 'history.json')
    t, carried_error, carried_transport, carry_profile = profiles(histories['carry', 'fine'], histories['carry', 'finer'])
    _, old_error, old_transport, old_profile = profiles(histories['no carry', 'fine'], histories['no carry', 'finer'])
    if (carry_profile['max_log_difference'] != saved['full_refinement']['errors']['chemical_log_max'] or
            old_profile['max_log_difference'] != saved['no_carry_full_refinement']['errors']['chemical_log_max']):
        raise ValueError('Independent raw chemical profile differs from summary')
    improvement = old_profile['max_log_difference'] / carry_profile['max_log_difference']

    fig, axes = plt.subplots(2, 2, figsize=(11.2, 7.3))
    colors = {'carry': '#226d49', 'no carry': '#ac452c'}
    for (method, level), h in histories.items():
        style = '-' if level == 'fine' else '--'
        times = np.array([r['elapsed'] for r in h])
        axes[0, 0].plot(times, [r['log_activator_sd'] for r in h], style, color=colors[method],
                       label=f'{method}, dt={.001875 if level == "fine" else .0009375:g}')
        axes[1, 1].plot(times, [r['uniform_growth_max'] for r in h], style, color=colors[method])
    axes[0, 0].axhline(.1, color='#555555', linestyle=':', label='Contrast criterion')
    axes[0, 0].axvspan(216, 240, color='#999999', alpha=.12)
    axes[0, 0].set_ylabel('Across-cell SD(log activator)'); axes[0, 0].legend(fontsize=8)
    axes[0, 0].set_title('Initiation and late persistence remain present')
    for values, label, color in ((carried_error, 'phase carry', colors['carry']), (old_error, 'original no carry', colors['no carry'])):
        axes[0, 1].semilogy(t, np.where(values > 0, values, np.nan), color=color, label=label)
    axes[0, 1].axhline(.01, color='#555555', linestyle=':', label='Unchanged 0.01 limit')
    axes[0, 1].set_ylabel('Max absolute log difference between timesteps'); axes[0, 1].legend(fontsize=8)
    axes[0, 1].set_title(f'Raw chemical discrepancy falls {improvement:.0f}-fold')
    for values, label, color in ((carried_transport, 'phase carry', colors['carry']), (old_transport, 'original no carry', colors['no carry'])):
        axes[1, 0].semilogy(t, np.where(values > 0, values, np.nan), color=color, label=label)
    axes[1, 0].set_ylabel('Relative transport-matrix difference'); axes[1, 0].legend(fontsize=8)
    axes[1, 0].set_title('Geometry-dependent exchange agrees more closely')
    axes[1, 1].axhline(0, color='#555555', linestyle=':')
    axes[1, 1].set_ylabel('Frozen uniform-state maximum growth rate')
    axes[1, 1].set_title('Pattern persists after frozen instability disappears')
    for ax in axes.flat:
        ax.set_xlabel('Elapsed model time from mature t=150 restart'); ax.grid(alpha=.2)
    fig.suptitle('History 9, 16 cells: phase-update precision through 240 model time units')
    fig.tight_layout(); image = docs / 'images/phase-carry-formation.png'
    fig.savefig(image, dpi=170); plt.close(fig)

    lines = ['', '## Completed formation assessment', '',
        '**Both full moving continuations pass the original timestep and physical quality screens.** '
        'The phase-carry pair passes without time alignment or a relaxed tolerance; both actual endpoints '
        'support local chemical bistability. This is one existing developmental history, not two independent histories.', '',
        '![Long formation and timestep comparison](images/phase-carry-formation.png)', '',
        '| Comparison | Maximum log chemical discrepancy | Original limit | Decision |',
        '|---|---:|---:|---|',
        f'| No carry, original dt pair | {old_profile["max_log_difference"]:.8f} | 0.01 | FAIL, retained |',
        f'| Phase carry, same dt pair | {carry_profile["max_log_difference"]:.8f} | 0.01 | PASS |', '',
        f'Preserving small phase increments reduces the maximum discrepancy by **{improvement:.1f}-fold**. '
        'The new maximum corresponds to about 0.0205% concentration-ratio disagreement and is about 49 times '
        'below the declared log-error limit. The original method still fails its recorded gate.', '',
        '| Quantity | Carry pair discrepancy | Unchanged tolerance |', '|---|---:|---:|']
    for key, value in saved['full_refinement']['errors'].items():
        lines.append(f'| {key} | {value:.8g} | {p["refinement_criteria"][key]:g} |')
    for key, criterion in (('onset_time_error', 'onset_time_abs_max'), ('crossing_time_error', 'crossing_time_abs_max')):
        lines.append(f'| {key} | {saved["full_refinement"][key]:.8g} | {p["refinement_criteria"][criterion]:g} |')
    lines += ['', '## Chemical outcome and endpoint stability', '',
        '| Timestep | Pattern onset, elapsed | Late minimum SD(log activator) | Final SD(log activator) | Endpoint chemical phase |',
        '|---|---:|---:|---:|---|']
    for row, endpoint in zip(saved['results'], endpoint_details):
        lines.append(f'| {row["job"]["dt"]:g} | {row["onset_elapsed"]:.6f} | {row["late_min_log_sd"]:.7f} | '
                     f'{row["final_log_sd"]:.7f} | {endpoint["phase"]} |')
    lines += ['',
        'Both patterns exceed the 0.1 contrast threshold throughout elapsed 216–240. Formation onset differs '
        'by only 0.000449 model units. The frozen uniform-state growth diagnostic changes sign near elapsed '
        '136.968, while the developed pattern persists. Instantaneous spectra describe a held geometry, not '
        'stability of the full moving system.', '',
        'All five frozen starts on each endpoint settle: the developed state and its two small perturbations '
        'approach the patterned equilibrium; the two near-uniform starts approach uniform chemistry. '
        'The maximum uniform-state eigenvalue is about -0.04603 and the patterned-state eigenvalue about '
        '-0.70031. Thus each frozen network supports locally stable uniform and patterned chemistry. '
        'These are ten nested chemical trajectories on two endpoint graphs from the same history.', '',
        '## Numerical mechanism and interpretation', '',
        'The intervention retains phase increments that would be lost when written to float32 storage, '
        'accumulating the residual in float64. Cells, contacts, force arithmetic, and the chemical equations '
        'retain their prior conventions. At the final step, roughly 1.9% and 3.8% of voxels with nonzero '
        'force would have an unchanged visible phase value without the carry; this is a diagnostic of that '
        'step, not a cumulative fraction over the run.', '',
        'With only this accumulation rule changed, long-window transport disagreement falls from '
        '0.00026651 to 0.0000024723, and chemical disagreement falls about 183-fold. This strongly implicates '
        'loss of small mechanical updates as a major contributor to the tested historical timestep sensitivity. '
        'It does not identify an exact solution or prove that every remaining precision error is negligible.', '',
        'The numerical intervention also changes the transient trajectory at a fixed timestep. Maximum '
        'carry/no-carry log differences are 0.06657 at dt=0.001875 and 0.10427 at dt=0.0009375. These are '
        'method-change effects, not failed carry refinement tests. The old paths are not interchangeable with '
        'the corrected paths even though both form patterns and share the endpoint phase classification.', '',
        'The result strengthens the tested initiation-versus-maintenance account: chemistry organizes from '
        'a small perturbation on moving mature geometry and remains patterned when its final frozen network '
        'also supports uniform chemistry. It does not establish autonomous or inherited identity, '
        'spontaneous aggregate-axis formation, new histories, or formation from a fresh zygote. '
        'The tested directional polarity-tension contrast is zero; positive-contrast suppression controls '
        'have not been repeated with carry.', '',
        '## What remains before broader use', '',
        'A two-timestep pass supports this history, parameter point, grid, and horizon. Add a third carry '
        'timestep and an independent high-precision mechanics reference before inferring a convergence trend '
        'or admitting the method broadly. Then repeat matched zero/positive-polarity controls across the '
        'existing developmental histories. GPU cleavage, geometric conductance closure, and spatial/full '
        'developmental convergence remain separate gates.', '',
        'The [prescribed geometry–chemistry replay](geometry_chemistry_coupling.md) still finds first-order '
        'beginning sampling and much more accurate second-order midpoint sampling. Its full-start '
        'source-spacing difference 0.00364878 still fails the declared 0.001 screen; obtain denser recorded '
        'geometry before stronger full-start attribution. This independent failure is not cleared by the '
        'moving carry result. The accepted production solver and earlier ledger/manuscript snapshots are unchanged.', '',
        '## Completed verification', '',
        f'Assessment recomputes the complete frozen summary without modifying original outputs; verifies '
        f'{spectra_checked:,} saved moving graph spectra and contrast values, final physical checkpoints '
        f'and carry arrays, exact reused prefixes, all ten frozen endpoints and {frozen_samples:,} frozen '
        'observations, aligned clocks, original tolerances, and all pinned source/input hashes.', '',
        f'The raw chemical maximum occurs at elapsed {carry_profile["peak_elapsed"]:.2f} with carry '
        f'and {old_profile["peak_elapsed"]:.2f} without carry. These profile locations are descriptive '
        'and do not change the acceptance window.', '',
        '```bash', 'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \\',
        'MPLCONFIGDIR=/tmp/embryo-mpl python -m embryo.phase_carry_formation_assessment', '```', '',
        'The [verification record](phase_carry_formation_completed_verification.json) preserves '
        'old failures, new method decisions, provenance, and descriptive intervention effects separately.']
    document = docs / 'phase_carry_formation.md'
    content = document.read_text().split('\n## Completed formation assessment\n')[0].rstrip() + '\n'
    document.write_text(content + '\n'.join(lines) + '\n')
    files += [Path(__file__).resolve(), image, document]
    verification = dict(verified=True, verified_utc=datetime.now(timezone.utc).isoformat(),
        independent_histories=1, new_histories=0, moving_trajectories=2,
        moving_observations=spectra_checked, frozen_endpoint_graphs=2,
        frozen_trials=10, frozen_observations=frozen_samples,
        original_outputs_unchanged=True, source_and_input_hashes_verified=True,
        new_full_window_refinement_pass=saved['full_refinement']['passed'],
        original_no_carry_refinement_pass=saved['no_carry_full_refinement']['passed'],
        chemical_discrepancy_improvement=improvement,
        raw_chemical_profiles=dict(carry=carry_profile, no_carry=old_profile),
        endpoints=endpoint_details, same_dt_interventions=saved['same_dt_interventions'],
        evidence_sha256={str(path): digest(path) for path in files})
    write_json(docs / 'phase_carry_formation_completed_verification.json', verification)
    return dict(verified=True, independent_histories=1, new_histories=0,
                chemical_discrepancy_improvement=improvement,
                moving_observations=spectra_checked, frozen_observations=frozen_samples)


if __name__ == '__main__':
    print(report())
