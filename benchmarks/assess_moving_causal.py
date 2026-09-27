"""Describe the completed moving pilot without changing its acceptance criteria."""
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from embryo.resolution import write_json


def assess(output):
    output = Path(output)
    status = json.loads((output / 'status.json').read_text())
    if status['state'] != 'completed':
        raise ValueError('Assessment requires all four completed branches')
    comparison = json.loads((output / 'comparison.json').read_text())
    p = comparison['protocol']
    reports = {arm: json.loads((output / arm / 'analysis.json').read_text()) for arm in p['arms']}
    histories = {arm: report['history'] for arm, report in reports.items()}
    expected_times = p['start'] + np.arange(round(p['duration']/p['interval'])+1)*p['interval']
    initial = histories['full'][0]
    summary = {}
    for arm, history in histories.items():
        np.testing.assert_allclose([r['time'] for r in history], expected_times, atol=1e-9, rtol=0)
        for key in ('axis_ratio', 'activator', 'inhibitor', 'fate', 'covariance'):
            np.testing.assert_array_equal(history[0][key], initial[key])
        with np.load(output / arm / 'final_state.npz', allow_pickle=False) as final:
            assert str(final['causal_arm']) == arm
        late = [r for r in history if r['time'] >= expected_times[-1]-p['late_duration']-1e-10]
        summary[arm] = {
            'initial_axis_ratio': history[0]['axis_ratio'],
            'final_axis_ratio': history[-1]['axis_ratio'],
            'late_mean_axis_ratio': comparison['late_axis_ratios'][arm],
            'late_contrast_min': min(r['activator_std'] for r in late),
            'late_contrast_max': max(r['activator_std'] for r in late),
            'final_fate_a': history[-1]['fate_a'], 'final_fate_b': history[-1]['fate_b'],
            'final_uncommitted': history[-1]['uncommitted'],
            'late_fate_separation_mean': (float(np.mean([r['fate_separation'] for r in late if r['fate_separation'] is not None])) if any(r['fate_separation'] is not None for r in late) else None),
            'max_cell_volume_error': reports[arm]['max_volume_error'],
            'min_radius_grid_cells': reports[arm]['min_radius'],
            'max_clipping': reports[arm]['max_clipping'],
            'max_sampled_boundary_occupancy': max(r['boundary_occupancy'] for r in history),
            'contact_components_range': [min(r['contact_components'] for r in history), max(r['contact_components'] for r in history)],
            'late_unstable_graph_modes_range': [min(r['unstable_graph_modes'] for r in late), max(r['unstable_graph_modes'] for r in late)],
        }
    result = {'input_sha256': {str(output / arm / 'analysis.json'): hashlib.sha256((output / arm / 'analysis.json').read_bytes()).hexdigest() for arm in p['arms']},
              'trajectory_coverage_and_matched_initial_states': True,
              'arms': summary, 'original_comparison': comparison}
    write_json(output / 'assessment.json', result)
    lines = ['# Assessment of moving-geometry controls', '',
             f"All four branches completed the unchanged t={p['start']:g}–{expected_times[-1]:g} protocol. Recorded observation times and initial shape, regulator, and fate values match across arms.", '',
             '| Arm | Final axis ratio | Late mean ratio | Late activator SD range | Final A / B / uncommitted |',
             '|---|---:|---:|---:|---:|']
    for arm, r in summary.items():
        lines.append(f"| {arm} | {r['final_axis_ratio']:.6f} | {r['late_mean_axis_ratio']:.6f} | {r['late_contrast_min']:.6g}–{r['late_contrast_max']:.6g} | {r['final_fate_a']} / {r['final_fate_b']} / {r['final_uncommitted']} |")
    lines += ['', 'Initial axis ratio: '+f"{initial['axis_ratio']:.6f}.", '',
              '## Predeclared decisions', '']
    for key, value in comparison['causal_effect_checks'].items():
        lines.append(f'- {key}: {"PASS" if value else "FAIL"}')
    lines += [f"- All numerical-quality checks: {'PASS' if comparison['all_quality_pass'] else 'FAIL'}", '',
              'Full minus control late mean axis ratios:']
    for arm, value in comparison['shape_excess_over_controls'].items():
        lines.append(f'- {arm}: {value:+.6f} (predeclared elongation threshold 0.05 for no-feedback and no-self-activation comparisons).')
    lines += ['', '## Interpretation', '']
    checks = comparison['causal_effect_checks']
    if not comparison['all_quality_pass']:
        lines.append('At least one declared numerical-quality check fails. Treat the intervention differences as provisional; do not count them as a validated causal shape result.')
    elif checks['feedback_specific_shape_excess'] and checks['self_activation_specific_shape_excess']:
        lines.append('The full branch exceeds both the mechanical-feedback and self-activation controls by the predeclared elongation margin. This supports a feedback-dependent shape response within this pilot, subject to moving-time refinement and replication.')
    else:
        lines.append('This pilot does not satisfy the combined criterion for feedback-specific and self-activation-specific elongation. An elongated full branch alone is insufficient: all branches inherited the same already anisotropic starting geometry.')
    if not checks['full_persistent_signal_contrast']:
        lines.append('Full-loop contrast does not remain above 0.1 throughout the late window. The frozen-geometry pattern result therefore does not establish persistent pattern formation under these moving conditions over the tested duration.')
    else:
        lines.append('Full-loop contrast stays above 0.1 throughout the late window; this is a signaling outcome separate from shape and fate labels.')
    lines.append('Final fate counts are reported separately because the bistable switch can amplify transient regulator differences. A/B labels alone do not demonstrate sustained activator–inhibitor patterning.')
    lines += ['', '## Numerical quality', '',
              '| Arm | Max relative volume error | Min radius / dx | Max boundary occupancy | Max clipping |',
              '|---|---:|---:|---:|---:|']
    for arm, r in summary.items():
        lines.append(f"| {arm} | {r['max_cell_volume_error']:.6g} | {r['min_radius_grid_cells']:.6g} | {r['max_sampled_boundary_occupancy']:.6g} | {r['max_clipping']:.6g} |")
    lines += ['', 'These checks assess one pre-existing mature geometry and one chemical seed. Shape differences are intervention responses, not evidence of a spontaneously selected developmental axis. A failed elongation criterion does not exclude other feedback effects or effects in another regime. Frozen-graph fate accuracy does not establish moving-geometry time convergence. The developmental refinement failures and geometric transport approximation remain separate limitations.', '']
    (output / 'ASSESSMENT.md').write_text('\n'.join(lines))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)
    labels = {'full': 'Full feedback', 'no_feedback': 'No mechanical feedback', 'no_self_activation': 'No self-activation', 'no_signal_to_fate': 'No signal-to-fate'}
    for arm, history in histories.items():
        t = [r['time'] for r in history]
        axes[0,0].plot(t, [r['axis_ratio'] for r in history], label=labels[arm])
        axes[0,1].semilogy(t, [max(r['activator_std'], 1e-16) for r in history])
        axes[1,0].plot(t, [r['fate_a']/r['cells'] for r in history])
        axes[1,1].plot(t, [r['max_cell_volume_error'] for r in history])
    axes[0,0].set_ylabel('Principal axis ratio')
    axes[0,1].set_ylabel('Activator standard deviation')
    axes[0,1].axhline(p['criteria']['persistent_contrast_min'], color='gray', ls='--', lw=1)
    axes[1,0].set_ylabel('Fraction labeled A')
    axes[1,1].set_ylabel('Maximum relative cell-volume error')
    axes[1,1].axhline(p['criteria']['volume_max'], color='gray', ls='--', lw=1)
    for ax in axes.flat:
        ax.axvspan(max(expected_times[0], expected_times[-1]-p['late_duration']), expected_times[-1], color='gray', alpha=.1)
        ax.set_xlim(expected_times[0], expected_times[-1])
        ax.set_xlabel('Developmental time'); ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=8)
    fig.suptitle('Matched moving-geometry pilot: one geometry and seed; shaded late comparison window')
    fig.savefig(output / 'comparison.png', dpi=180)
    fig.savefig(output / 'comparison.pdf')
    plt.close(fig)
    return result


if __name__ == '__main__':
    assess(sys.argv[1] if len(sys.argv)>1 else 'outputs/moving-causal')
