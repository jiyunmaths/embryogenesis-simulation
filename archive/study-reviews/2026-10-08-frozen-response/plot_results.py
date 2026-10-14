"""Plot independently reconstructed moving/frozen comparisons; no simulation."""
from pathlib import Path
import csv
import hashlib
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import numpy as np

FOLDER = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    data = json.loads((FOLDER / 'independent_review.json').read_text())
    assert data['passed'] and data['inputs_unchanged']
    pairs = [(r['seed'], r['cell']) for r in data['baseline_descriptions']]
    labels = [f'H{seed}\ncell {cell}' for seed, cell in pairs]
    colors = ['#386d8a', '#e4a04b']
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False,
                         'axes.spines.right': False})
    fig = plt.figure(figsize=(12, 8.4), layout='constrained')
    grid = fig.add_gridspec(2, 2, height_ratios=[1, .95])
    axes = [fig.add_subplot(grid[0, i]) for i in range(2)]
    axes.append(fig.add_subplot(grid[1, :]))
    x = np.arange(len(pairs)); width = .35
    for axis, dose, title in zip(axes[:2], ('fractional', 'amount'),
            ('A   10% fractional pulse', 'B   Common removed amount')):
        records = [next(r for r in data['scientific_comparisons']
                        if (r['seed'], r['cell'], r['dose']) == (*pair, dose))
                   for pair in pairs]
        for shift, label, color, key in zip((-.5, .5), ('Moving', 'Frozen'),
                                            colors, ('moving-fine', 'fine')):
            axis.bar(x + shift*width, [r['full_effects'][key] for r in records],
                     width=width, color=color, label=label)
        axis.set_yscale('log'); axis.set_ylim(3e-5, .17)
        axis.axhline(.01, color='#666666', linestyle='--', linewidth=1, label='Screen = 0.01')
        axis.set_xticks(x, labels); axis.set_title(title, loc='left', fontweight='bold')
        axis.grid(axis='y', alpha=.16); axis.set_axisbelow(True)
        axis.set_ylabel('Full 60-unit response difference (RMS)')
        axis.legend(frameon=False, loc='upper left', ncol=3, fontsize=9)
    axis = axes[2]
    for shift, label, color, key in zip((-.5, .5), ('Moving', 'Frozen'), colors,
                                       ('moving-fine', 'fine')):
        axis.bar(x + shift*width,
                 [r[key]['final_activator_redistribution_sham']
                  for r in data['baseline_descriptions']],
                 width=width, color=color, label=label)
    axis.axhline(1., color='#666666', linestyle='--', linewidth=1)
    axis.set_yscale('log'); axis.set_ylim(.035, 1.9)
    axis.set_yticks([.05, .1, .2, .5, 1.]); axis.yaxis.set_major_formatter(ScalarFormatter())
    axis.set_xticks(x, labels); axis.grid(axis='y', alpha=.16); axis.set_axisbelow(True)
    axis.set_ylabel('Final activator: redistributed / sham')
    axis.set_title('C   Unperturbed chemical levels at t=570', loc='left', fontweight='bold')
    axis.text(.99, .97, 'Cells 20 / 30: frozen activator is 13.6 / 17.5-fold lower',
              transform=axis.transAxes, ha='right', va='top', fontsize=9)
    axis.legend(frameon=False, loc='upper left', ncol=2)
    fig.suptitle('Delayed response differences persist when physical evolution stops\n'
                 'Actual t=510 starts; fine dt=0.001875; six recipients nested in three histories',
                 fontsize=14, fontweight='bold')
    for suffix in ('png', 'svg'):
        fig.savefig(FOLDER / f'frozen-response.{suffix}', dpi=180,
                    metadata={'Description': 'From an independent raw-array reconstruction; no simulation or post hoc screen.'})
    plt.close(fig)
    rows = []
    for seed, cell in pairs:
        r = dict(history=seed, cell=cell)
        for dose in ('fractional', 'amount'):
            item = next(t for t in data['scientific_comparisons']
                        if (t['seed'], t['cell'], t['dose']) == (seed, cell, dose))
            r.update({f'{dose}_moving':item['full_effects']['moving-fine'],
                      f'{dose}_frozen':item['full_effects']['fine'],
                      f'{dose}_above_screen':item['frozen_detected']})
        rows.append(r)
    with (FOLDER / 'recipient-summary.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    report = dict(input_sha256=sha(FOLDER / 'independent_review.json'),
                  script_sha256=sha(Path(__file__).resolve()),
                  outputs_sha256={name:sha(FOLDER / name) for name in
                                  ('frozen-response.png','frozen-response.svg','recipient-summary.csv')},
                  matplotlib=matplotlib.__version__,
                  scope='Fine dt full60 target two-species own-control-normalized response differences, '
                        'plus final unperturbed activator ratios. Six recipients are nested in three histories; '
                        'no statistical error bars or independently replicated dose interpretation.')
    (FOLDER / 'figure-provenance.json').write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
