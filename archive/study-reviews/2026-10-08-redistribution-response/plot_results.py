"""Plot independently reconstructed data; no simulation or model helper is run."""
from pathlib import Path
import hashlib
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

folder = Path(__file__).resolve().parent
audit_path = folder / 'response_reconstruction.json'
data_path = folder / 'plot-data.json'
audit = json.loads(audit_path.read_text())
data = json.loads(data_path.read_text())
assert audit['independent_response_reconstruction_passed']
assert hashlib.sha256(data_path.read_bytes()).hexdigest() == audit['plot_data_sha256']
assert data['protocol_sha256'] == audit['protocol_sha256']
rows = audit['recipients']
labels = [f"{r['history']} / {r['cell']}" for r in rows]
x = np.arange(len(rows))
blue, orange = '#2469a4', '#c36b2c'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                     'axes.titlesize': 11, 'axes.labelsize': 10,
                     'svg.fonttype': 'none'})
fig, axes = plt.subplots(2, 2, figsize=(11.4, 8.5), layout='constrained')
ax = axes[0, 0]
for offset, dose, color, label in [(-.18, 'fractional', blue, '10% fractional depletion'),
                                   (.18, 'amount', orange, 'Common removed amount')]:
    values = [r['by_timestep']['fine'][dose]['normalized_response_shift_rms'] for r in rows]
    ax.bar(x + offset, values, width=.32, color=color, label=label)
ax.axhline(audit['effect_screen'], color='#333333', linestyle='--', linewidth=1,
           label='Prespecified screen: 0.01')
ax.set_yscale('log')
ax.set_ylim(2e-5, .5)
ax.set_xticks(x, labels)
ax.set_xlabel('History / recipient')
ax.set_ylabel('Two-species normalized response RMS')
ax.set_title('A  Delayed response effects, full 60-unit window', loc='left')
ax.legend(fontsize=8, loc='upper left', frameon=False)
ax.grid(axis='y', alpha=.15)

ax = axes[0, 1]
ratios = [r['fine_baseline']['final_species_ratios'][0] for r in rows]
ax.scatter(x, ratios, s=50, color=blue, zorder=3)
ax.axhline(1, color='#333333', linestyle='--', linewidth=1)
ax.set_yscale('log')
ax.set_ylim(.04, 1.7)
ax.set_xlim(-.5, 5.65)
ax.set_xticks(x, labels)
ax.set_xlabel('History / recipient')
ax.set_ylabel('Unperturbed activator: redistributed / sham')
ax.set_title('B  Baseline levels at t = 570', loc='left')
for i, value in enumerate(ratios):
    label = f'{value:.3f}' if value > .1 else f'{1/value:.1f}-fold lower'
    ax.annotate(label, (i, value), xytext=(0, 10), textcoords='offset points',
                ha='center', fontsize=8)
ax.grid(axis='y', alpha=.15)

for ax, target, panel in [(axes[1, 0], 19, 'C'), (axes[1, 1], 20, 'D')]:
    record = next(r for r in data['data'] if r['history'] == 8 and r['cell'] == target)
    times = np.array(record['times'])
    for context, color, context_label in [('sham', blue, 'Sham'),
                                           ('redistributed', orange, 'Redistributed')]:
        wave = np.array(record['doses']['fractional'][context]['normalized_two_species_response'])
        assert wave.shape == (401, 2)
        for species, linestyle, name in [(0, '-', 'activator'), (1, '--', 'inhibitor')]:
            ax.plot(times, wave[:, species], color=color, linestyle=linestyle, linewidth=1.6,
                    label=f'{context_label}, {name}')
    ax.axhline(0, color='#777777', linewidth=.7)
    ax.set_xlim(0, 8)
    ax.set_xlabel('Model time after t = 510 pulse')
    ax.set_ylabel('Log pulse / own control, normalized')
    ax.set_title(f'{panel}  History 8, cell {target}: fractional pulse', loc='left')
    ax.legend(fontsize=8, loc='lower right', frameon=False)
    ax.grid(alpha=.15)

fig.suptitle('Initial chemical allocation affects maintained state and selective later responses\n'
             'Three existing mature histories; doses and recipients are nested comparisons', fontsize=13)
fig.savefig(folder / 'redistribution-response.png', dpi=180)
fig.savefig(folder / 'redistribution-response.svg')
plt.close(fig)
record = {'protocol_sha256': audit['protocol_sha256'],
          'response_reconstruction_sha256': hashlib.sha256(audit_path.read_bytes()).hexdigest(),
          'plot_data_sha256': hashlib.sha256(data_path.read_bytes()).hexdigest(),
          'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
          'files_sha256': {name: hashlib.sha256((folder / name).read_bytes()).hexdigest()
                           for name in ('redistribution-response.png', 'redistribution-response.svg')},
          'scope': 'Fine dt=0.001875. Panel A uses both species over 60 units; '
                   'panel B uses own-control final activator ratios; panels C/D display '
                   'the first eight units of fractional-pulse two-species responses. '
                   'All effects are computed over the full prespecified window. '
                   'Original amounts matched at t450, not at challenge time t510.'}
(folder / 'figure-provenance.json').write_text(json.dumps(record, indent=2) + '\n')
print('Created PNG and editable SVG from independent reconstruction data.')
