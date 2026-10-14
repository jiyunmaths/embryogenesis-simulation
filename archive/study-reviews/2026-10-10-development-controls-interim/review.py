"""Read-only interim checks; no GPU work, field reconstruction or acceptance."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import os

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
from embryo import qualify_repeated_development as q
from embryo import development_space_refinement as space

ROOT = Path('output/studies/development-control-qualification')
OUT = Path(__file__).parent


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    protocol = read(ROOT / 'protocol.json')
    protocol_hash = sha(ROOT / 'protocol.json')
    sources_verified = all(sha(name) == value for name, value in protocol['source_sha256'].items())
    rows = []
    reports = {}
    for job in protocol['jobs']:
        if job['reused']:
            continue
        folder = ROOT / job['key']
        result = folder / 'gpu-result.json'
        if not result.exists():
            progress = folder / 'gpu-progress.json'
            rows.append(dict(case=job['key'], state='active' if progress.exists() else 'pending',
                             progress=read(progress) if progress.exists() else None))
            continue
        report = read(result)
        history = read(folder / 'gpu-history.json')
        small_pins = {name: value for name, value in report['evidence_sha256'].items()
                      if Path(name).suffix in ('.json', '.jsonl')}
        verified = all(sha(name) == value for name, value in small_pins.items())
        observation_pass = all(
            x['boundary_occupancy'] < protocol['physical_limits']['boundary_max']
            and all(n == 1 for n in x['nondividing_components'])
            and all(value > 0 for species in x['chemistry'] for value in species)
            for x in history)
        physical_pass = q.path_quality(report['audit'], protocol['physical_limits'])
        events_pass = all(q.event_passed(e, protocol['physical_limits']) for e in report['events'])
        provenance_pass = report['job'] == job and report['protocol_sha256'] == protocol_hash
        assert all((verified, observation_pass, physical_pass, events_pass, provenance_pass))
        last = history[-1]
        rows.append(dict(case=job['key'], state='completed', reported_passed=report['passed'],
                         seconds=report['seconds'], audit=report['audit'], events=len(report['events']),
                         cells=report['final_cells'], time=report['final_time'], history_rows=len(history),
                         runtime_summary_rechecked=True, small_evidence_hashes_verified=verified,
                         whole_checkpoint_hashes_rechecked=False,
                         restart_steps=report['exact_restart_steps'],
                         independent_chemical_log_max=report['independent_chemical_log_max'],
                         endpoint={key: last[key] for key in
                                   ('log_activator_sd', 'axis_ratio', 'boundary_occupancy', 'dividing_ids')}))
        reports[job['key']] = report

    diagnostics = []
    for branch in ('tension_adhesion', 'polarity', 'full'):
        pairs = [('space', pair) for pair in protocol['space_pairs'][branch]]
        pairs.append(('time', protocol['time_pairs'][branch]))
        for kind, pair in pairs:
            if not all(name in reports for name in pair):
                continue
            errors, unmatched = q.compare_observations(
                *[read(ROOT / name / 'gpu-history.json') for name in pair])
            limits = protocol['space_limits' if kind == 'space' else 'time_limits']
            maturity = abs(max(e['time'] for e in reports[pair[0]]['events'])
                           - max(e['time'] for e in reports[pair[1]]['events']))
            failures = {key: value for key, value in errors.items() if value > limits[key]}
            if maturity > limits['maturity_time_difference_max']:
                failures['maturity_time_difference_max'] = maturity
            diagnostics.append(dict(kind=kind, cases=pair, observation_errors=errors,
                                    unmatched_times=unmatched, maturity_difference=maturity,
                                    observation_only_failures=failures,
                                    dense_spectrum=space.dense_pair(ROOT, *pair),
                                    field_errors_assessed=False,
                                    full_comparison_accepted=False))
    result = dict(utc=datetime.now(timezone.utc).isoformat(), protocol_sha256=protocol_hash,
                  source_pins_verified=sources_verified, source_pins_count=len(protocol['source_sha256']),
                  status=read(ROOT / 'status.json'), paths=rows, partial_diagnostics=diagnostics,
                  production_completed=sum(r['state'] == 'completed' for r in rows), production_total=12,
                  full_replays_completed=len(list(ROOT.glob('replays/*/stage-*/result.json'))),
                  refinement_reports_present=[p.name for p in ROOT.glob('*-refinement.json')],
                  qualification_accepted=False, new_independent_scientific_histories=0,
                  note='Interim read-only observation checks exclude field reconstruction and full checkpoint hash review.')
    target = OUT / ('interim-' + datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '.json')
    target.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(dict(output=str(target), utc=result['utc'], completed=result['production_completed'],
                         status=result['status'], source_pins_verified=sources_verified,
                         diagnostics=[dict(cases=d['cases'], spectrum=d['observation_errors']['relative_spectrum_max'],
                                           dense_peak=d['dense_spectrum']['peak'],
                                           failures=d['observation_only_failures']) for d in diagnostics]), indent=2))


if __name__ == '__main__':
    main()
