"""Frozen-geometry sensitivity of conservative signaling to contact filtering."""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path

import numpy as np

from .model import Simulation
from .resolution import write_json
from .signaling import stability
from .transport import contact_transport, transport_graph


def scan(sim, cutoffs):
    if sim.config.signal_transport != 'conservative':
        raise ValueError('requires conservative transport')
    contacts = sim.contacts()[0]
    volumes, centers = sim.volumes(), sim.centers()
    rows = []
    baseline = transport_graph(contact_transport(contacts, volumes, centers, sim.config.interface_width, .02))
    base = stability(baseline, sim.config.signal_beta, sim.config.signal_da, sim.config.signal_dh)
    for cutoff in cutoffs:
        graph = transport_graph(contact_transport(contacts, volumes, centers, sim.config.interface_width, cutoff))
        result = stability(graph, sim.config.signal_beta, sim.config.signal_da, sim.config.signal_dh)
        scale = np.linalg.norm(baseline.weights)
        rates = np.asarray(result['growth_rates'])
        rows.append({'cutoff': cutoff, **result,
                     'conductance_relative_change': float(np.linalg.norm(graph.weights-baseline.weights)/scale) if scale else 0.,
                     'spectrum_relative_l2_change': float(np.linalg.norm(graph.eigenvalues-baseline.eigenvalues)/max(np.linalg.norm(baseline.eigenvalues),1e-30)),
                     'maximum_growth_rate_change': float(np.max(abs(rates-np.asarray(base['growth_rates'])))),
                     'component_count_unchanged': result['components'] == base['components'],
                     'unstable_mode_count_unchanged': len(result['unstable_modes']) == len(base['unstable_modes']),
                     'constant_residual': float(np.max(abs(graph.delta @ np.ones(len(volumes))))),
                     'amount_conservation_residual': float(np.max(abs(volumes @ graph.delta)))})
    return rows


def assess(rows, criteria):
    local = [r for r in rows if r['cutoff'] in (.01,.02,.04)]
    return {'local_connectivity_unchanged': all(r['component_count_unchanged'] for r in local),
            'local_unstable_mode_count_unchanged': all(r['unstable_mode_count_unchanged'] for r in local),
            'local_conductance_change_below_5_percent': all(r['conductance_relative_change'] < criteria['conductance_relative_max'] for r in local),
            'local_spectrum_change_below_5_percent': all(r['spectrum_relative_l2_change'] < criteria['spectrum_relative_max'] for r in local),
            'local_growth_change_below_0_01': all(r['maximum_growth_rate_change'] < criteria['growth_absolute_max'] for r in local),
            'constant_and_amount_residuals_below_1e_10': all(max(r['constant_residual'],r['amount_conservation_residual']) < 1e-10 for r in rows)}


def run(checkpoints, output):
    output = Path(output)
    if output.exists():
        raise FileExistsError('choose a fresh output directory')
    if not checkpoints:
        raise ValueError('provide at least one checkpoint')
    sources = [Path(p).resolve() for p in checkpoints]
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    protocol = {'source_sha256': hashes, 'baseline_cutoff': .02,
                'cutoffs': [0.,.005,.01,.02,.04,.08], 'acceptance_neighborhood': [.01,.02,.04],
                'criteria': {'conductance_relative_max': .05, 'spectrum_relative_max': .05, 'growth_absolute_max': .01},
                'scope': 'Frozen geometry and homogeneous-equilibrium linearization only. Extreme cutoffs are stress diagnostics, not acceptance cases. Sorted eigenvalue comparisons do not track eigenvector identity. Does not test evolving feedback, polarity sensitivity, or validate physical contact areas. Coarse developmental snapshots retain their cell-resolution limitation.'}
    output.mkdir(parents=True)
    write_json(output/'protocol.json', protocol)
    write_json(output/'status.json', {'state':'running'})
    try:
        cases = []
        for path in sources:
            sim = Simulation.restore(path)
            rows = scan(sim, protocol['cutoffs'])
            cases.append({'source': str(path), 'time': sim.time, 'config': asdict(sim.config),
                          'rows': rows, 'checks': assess(rows,protocol['criteria'])})
            print(f'{path}: {cases[-1]["checks"]}', flush=True)
        unchanged = all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in hashes.items())
        result = {'protocol':protocol, 'cases':cases, 'sources_unchanged':unchanged,
                  'passed':unchanged and all(all(c['checks'].values()) for c in cases)}
        write_json(output/'comparison.json', result)
        lines = ['# Frozen contact-cutoff sensitivity', '', f'All declared checks pass: **{result["passed"]}**', '',
                 '| State | Time | Cutoff | Edges | Components | Unstable modes | Conductance change | Max growth change |',
                 '|---|---:|---:|---:|---:|---:|---:|---:|']
        for c in cases:
            for r in c['rows']:
                lines.append(f'| {Path(c["source"]).parent.name} | {c["time"]:.3g} | {r["cutoff"]} | {r["edges"]} | {r["components"]} | {len(r["unstable_modes"])} | {r["conductance_relative_change"]:.3%} | {r["maximum_growth_rate_change"]:.5g} |')
        lines += ['', protocol['scope'], '', 'Failed local checks:']
        lines += [f'- t={c["time"]:g}, {c["source"]}: {key}' for c in cases for key,value in c['checks'].items() if not value]
        (output/'RESULTS.md').write_text('\n'.join(lines)+'\n')
        write_json(output/'status.json', {'state':'completed', 'all_checks_pass':result['passed']})
        return result
    except Exception as error:
        write_json(output/'status.json', {'state':'failed', 'error':str(error)})
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--checkpoints', nargs='+', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    run(args.checkpoints, args.output)


if __name__ == '__main__':
    main()
