"""Plot saved finer-clock chemistry and conductance without solver imports."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path('output/studies/fresh-zygote-formation')
HERE = Path(__file__).resolve().parent
BRANCHES = ['off', 'tension_adhesion', 'polarity', 'full']
LABELS = ['Feedback off', 'Tension + adhesion', 'Polarity tension', 'Full coupling']
COLORS = ['#2563eb', '#555b66', '#d97706', '#15803d']

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.spines.top': False,
                     'axes.spines.right': False, 'svg.fonttype': 'none'})
fig, axes = plt.subplots(2, 3, figsize=(12.8, 7.4), sharex=True)
summary = []
for column, history in enumerate([1101, 1102, 1103]):
    for branch, label, color in zip(BRANCHES, LABELS, COLORS):
        directory = ROOT / f'h{history}_{branch}_fine'
        rows = json.loads((directory / 'gpu-history.json').read_text())
        times = np.array([row['elapsed'] for row in rows])
        contrast = np.array([np.std(np.log(np.array(row['chemistry'])[0])) for row in rows])
        conductance = np.array([np.triu(np.array(row['volumes'])[:, None] * np.array(row['delta']), 1).sum()
                               for row in rows])
        axes[0, column].plot(times, contrast, color=color, lw=1.9, label=label)
        axes[1, column].plot(times, conductance / conductance[0], color=color, lw=1.9)
        recorded = json.loads((directory / 'gpu-result.json').read_text())
        primary = next(x for x in recorded['formation'] if x['threshold'] == .1)
        summary.append(dict(history=history, branch=branch, primary=primary,
            final_contrast=float(contrast[-1]), initial_contrast=float(contrast[0]),
            conductance_decline_percent=float(100 * (1 - conductance[-1] / conductance[0])),
            initial_modes=rows[0]['frozen_uniform_unstable_modes'], final_modes=rows[-1]['frozen_uniform_unstable_modes'],
            final_uniform_growth=rows[-1]['frozen_uniform_max_growth'], final_axis_ratio=rows[-1]['axis_ratio']))
    axes[0, column].set_title(f'History {history}', pad=10)
    axes[0, column].set_yscale('log')
    axes[0, column].set_ylim(.001, 2.0)
    axes[0, column].axhline(.1, color='#64748b', lw=1, ls=':', zorder=0)
    axes[0, column].axvspan(36, 60, color='#64748b', alpha=.08, zorder=0)
    axes[1, column].set_ylim(.4, 1.03)
    axes[1, column].set_xlabel('Elapsed model time since t = 18')
    for axis in axes[:, column]:
        axis.set_xlim(0, 60)
        axis.set_xticks([0, 12, 24, 36, 48, 60])
        axis.grid(axis='y', alpha=.18)
axes[0, 0].set_ylabel('Across-cell SD of ln(activator)')
axes[1, 0].set_ylabel('Total conductance / starting conductance')
fig.suptitle('Inherited developmental differences form sustained chemical patterns', fontsize=16, y=.985)
fig.legend(*axes[0, 0].get_legend_handles_labels(), loc='upper center', ncol=4,
           bbox_to_anchor=(.5, .94), frameon=False)
fig.subplots_adjust(top=.83, bottom=.15, left=.08, right=.98, hspace=.28, wspace=.25)
fig.text(.5, .025, 'Three independent histories · Four nested mechanical contexts · 16 cells per context\n'
         'Dotted: primary threshold 0.1. Shaded: late persistence window (elapsed 36–60).\n'
         'Finer clock dt = 0.001875; physical t = 18–78. Chemical contrast is not a biological identity label.',
         ha='center', fontsize=9, color='#475569')
fig.savefig(HERE / 'formation.png', dpi=180)
fig.savefig(HERE / 'formation.svg')
(HERE / 'descriptive_fine_results.json').write_text(json.dumps(summary, indent=2) + '\n')
print('Saved formation.png, formation.svg and descriptive_fine_results.json')
