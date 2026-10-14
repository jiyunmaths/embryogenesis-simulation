"""Plot descriptive fine-clock observations without making an acceptance claim."""
from pathlib import Path
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path('output/studies/fresh-zygote-mature-clock-fix')
HERE = Path(__file__).resolve().parent
COLORS = {'off': '#616161', 'tension_adhesion': '#0072b2',
          'polarity': '#e69f00', 'full': '#009e73'}
LABELS = {'off': 'Feedback off', 'tension_adhesion': 'Tension + adhesion',
          'polarity': 'Polarity', 'full': 'Full combination'}


def main():
    fig, axes = plt.subplots(2, 3, figsize=(12, 6.8), sharex=True, sharey='row')
    for col, history in enumerate((1101, 1102, 1103)):
        for branch, color in COLORS.items():
            folder = ROOT / f'h{history}_{branch}_fine'
            rows = json.loads((folder / 'gpu-history.json').read_text())
            network = [json.loads(line) for line in (folder / 'network.jsonl').read_text().splitlines()]
            t = [x['time'] for x in rows]
            sd = [np.log(np.asarray(x['chemistry'])[0]).std() for x in rows]
            g = np.array([np.triu(x['conductance'], 1).sum() for x in network])
            axes[0, col].plot(t, sd, color=color, label=LABELS[branch], linewidth=1.9)
            axes[1, col].plot([x['time'] for x in network], g / g[0], color=color, linewidth=1.9)
        axes[0, col].set_title(f'History {history}', fontsize=12)
        axes[1, col].set_xlabel('Physical model time')
        for ax in axes[:, col]:
            ax.set_xlim(18, 24)
            ax.set_xticks([18, 20, 22, 24])
            ax.grid(alpha=.2)
            ax.spines[['top', 'right']].set_visible(False)
    axes[0, 0].set_ylabel('SD of log activator')
    axes[1, 0].set_ylabel('Summed conductance / starting value')
    fig.suptitle('Short mature continuation: chemistry grows while contacts remodel', fontsize=15, y=.99)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(.5, .947), ncol=4, frameon=False)
    fig.text(.5, .018, 'Three existing developmental histories; four nested branches each. '
             'Fine dt = 0.001875; t = 18–24. This window does not establish persistence or cell identity.',
             ha='center', fontsize=9)
    fig.tight_layout(rect=(0, .045, 1, .89))
    for suffix in ('png', 'svg'):
        destination = HERE / f'mature-short-window.{suffix}'
        if destination.exists():
            raise FileExistsError(destination)
        fig.savefig(destination, dpi=170)
    plt.close(fig)


if __name__ == '__main__':
    main()
