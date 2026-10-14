"""Plot accepted saved-array measures; no simulation or new acceptance screen."""
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


review_path = ROOT / 'independent_review.json'
wave_path = ROOT / 'waveforms.npz'
review = json.loads(review_path.read_text())
assert review['passed'] and review['inputs_unchanged']
assert sha(wave_path) == review['waveform_sha256']
waves = np.load(wave_path)
effects = {(r['seed'], r['cell'], r['matching_reference'], r['physics']): r
           for r in review['scientific_comparisons']}
baselines = {(r['seed'], r['cell'], r['matching_reference'], r['physics']): r
             for r in review['baseline_descriptions']}
cells = [(7, 22), (7, 23), (8, 19), (8, 20), (9, 29), (9, 30)]
labels = [f'H{s}\ncell {c}' for s, c in cells]
colors = ['#2367a1', '#df852e']
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False,
                     'axes.spines.right': False, 'svg.fonttype': 'none'})
fig, axes = plt.subplots(2, 2, figsize=(12, 9.5))
fig.subplots_adjust(left=.09, right=.98, bottom=.13, top=.87,
                    hspace=.42, wspace=.29)
fig.suptitle('Surrounding chemistry affects cells after their initial state is matched',
             fontsize=16, x=.5, y=.965)
fig.text(.5, .92, 'Same initial target activator + inhibitor; shared transport and volumes; immediate 10% activator depletion',
         ha='center', fontsize=10)
x = np.arange(len(cells))
for ax, reference, panel in zip(axes[0], ('S', 'R'), ('A', 'B')):
    for i, physics in enumerate(('S', 'R')):
        xx = x + (i-.5)*.32
        values = [effects[s, c, reference, physics]['full_effects']['fine'] for s, c in cells]
        old = [effects[s, c, reference, physics]['parent_unmatched_fractional_rms'] for s, c in cells]
        ax.bar(xx, values, width=.28, color=colors[i],
               label=('Sham' if physics == 'S' else 'Redistributed')+' physics')
        ax.scatter(xx, old, marker='o', s=34, facecolors='white', edgecolors='#555555',
                   zorder=3, label='Prior unmatched target' if i == 0 else None)
    ax.axhline(.01, color='#7c3146', ls='--', lw=1.3)
    ax.text(5.55, .011, 'screen 0.01', ha='right', color='#7c3146', fontsize=9)
    ax.set_yscale('log')
    ax.set_ylim(1e-5, .5)
    ax.set_xticks(x, labels)
    ax.set_ylabel('Two-species response-difference RMS\n(normalized; full 60-unit window)')
    name = 'sham' if reference == 'S' else 'redistributed'
    ax.set_title(f'{panel}  Target matched to its measured {name} state', loc='left', fontsize=10.5)
    ax.legend(fontsize=8, loc='upper left', frameon=False, ncol=2)
    ax.grid(axis='y', alpha=.12)

ax = axes[1, 0]
for i, physics in enumerate(('S', 'R')):
    for reference, marker in [('S', 'o'), ('R', 'x')]:
        values = [baselines[s, c, reference, physics]['final_species_ratios'][0] for s, c in cells]
        ax.scatter(x + (i-.5)*.2, values, color=colors[i], marker=marker, s=46,
                   label=('Sham' if physics == 'S' else 'Redistributed')+f' physics, match {reference}')
ax.axhline(1, color='#777777', ls='--', lw=1)
ax.text(5.48, 1.13, 'All initial target ratios = 1', ha='right', fontsize=9)
ax.set_yscale('log')
ax.set_ylim(.035, 2.2)
ax.set_xticks(x, labels)
ax.set_ylabel('Final unpulsed target activator ratio\n(redistributed / sham surrounding chemistry)')
ax.set_title('C  Local concentration differences reappear', loc='left', fontsize=10.5)
ax.legend(fontsize=8, loc='lower left', frameon=False)
ax.grid(axis='y', alpha=.12)

ax = axes[1, 1]
times = waves['times'] - waves['times'][0]
for color, (seed, cell) in zip(colors, [(8, 20), (9, 30)]):
    for reference, style in [('S', '-'), ('R', '--')]:
        key = f'{seed}_{cell}_{reference}_fractional_chemical_sham_physics'
        ax.plot(times, waves[key][:, 0], color=color, ls=style, lw=1.8,
                label=f'H{seed} cell {cell}, match {reference}')
ax.axhline(0, color='#777777', lw=.8)
ax.set_xlim(0, 12)
ax.set_xlabel('Time after simultaneous matching and pulse')
ax.set_ylabel('Signed activator response difference\n(redistributed − sham; normalized)')
ax.set_title('D  Larger effects occur during the reset response', loc='left', fontsize=10.5)
ax.legend(fontsize=8, frameon=False)
ax.grid(alpha=.12)
fig.text(.09, .064, 'S/R matching: both target species set to the measured sham/redistributed pair. Panel D holds sham physics fixed.', fontsize=9)
fig.text(.09, .039, 'Three existing mature histories; recipients are nested within histories. Surrounding total chemical amounts remain unmatched.', fontsize=9)
fig.text(.09, .014, 'Endpoint agreement between matching references is descriptive. This immediate-reset assay does not establish settled response traits.', fontsize=9)
for suffix in ('png', 'svg'):
    fig.savefig(ROOT / f'recipient-state-response.{suffix}', dpi=180)
plt.close(fig)

csv_path = ROOT / 'recipient-summary.csv'
fields = ['history', 'cell', 'matching_reference', 'shared_physics', 'fine_response_rms',
          'above_primary_screen', 'parent_unmatched_fractional_rms',
          'initial_activator_ratio', 'initial_inhibitor_ratio', 'final_activator_ratio',
          'final_inhibitor_ratio', 'late_target_log_rms', 'late_network_log_rms']
with csv_path.open('w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()
    for (seed, cell, reference, physics), effect in effects.items():
        baseline = baselines[seed, cell, reference, physics]
        writer.writerow(dict(history=seed, cell=cell, matching_reference=reference,
            shared_physics=physics, fine_response_rms=effect['full_effects']['fine'],
            above_primary_screen=effect['above_screen'],
            parent_unmatched_fractional_rms=effect['parent_unmatched_fractional_rms'],
            initial_activator_ratio=baseline['initial_species_ratios'][0],
            initial_inhibitor_ratio=baseline['initial_species_ratios'][1],
            final_activator_ratio=baseline['final_species_ratios'][0],
            final_inhibitor_ratio=baseline['final_species_ratios'][1],
            late_target_log_rms=baseline['late_target_log_rms'],
            late_network_log_rms=baseline['late_network_log_rms']))
files = [review_path, wave_path, Path(__file__).resolve(), csv_path,
         ROOT/'recipient-state-response.png', ROOT/'recipient-state-response.svg']
(ROOT/'figure-provenance.json').write_text(json.dumps(dict(
    passed=True, source='Accepted raw-array audit and its saved fine response waveforms.',
    no_new_integrations=True, no_new_acceptance_threshold=True,
    endpoint_panel='All four reference/physical-context combinations shown; overlapping endpoint markers are retained.',
    files_sha256={str(p): sha(p) for p in files}), indent=2)+'\n')
print('Saved PNG/SVG, 24-row CSV and figure provenance.')
