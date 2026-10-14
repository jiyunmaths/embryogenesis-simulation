"""Read-only interim reconstruction; no scientific result is accepted here."""
from pathlib import Path
from datetime import datetime, timezone, timedelta
from zoneinfo import ZoneInfo
import hashlib
import json
import statistics
import shutil
import numpy as np

ROOT = Path('output/studies/moving-recipient-amount-response')
HERE = Path(__file__).parent
read = lambda p: json.loads(Path(p).read_text())
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
p = read(ROOT/'protocol.json')
protocol_sha = sha(ROOT/'protocol.json')
snapshot = read(ROOT/'status.json')
completed, histories, evidence, progress = [], {}, {}, {}
for job in p['jobs']:
    folder = ROOT/job['key']
    status_file = folder/'status.json'
    progress[job['key']] = read(status_file).get('elapsed', 0.) if status_file.exists() else 0.
    result_file = folder/'result.json'
    if not result_file.exists():
        continue
    result = read(result_file)
    assert result['job'] == job and result['protocol_sha256'] == protocol_sha
    assert result['quality_pass'] and result['method'] == p['method']
    for file, key in ((folder/'history.json', 'history_sha256'),
                      (folder/'latest_state.npz', 'checkpoint_sha256'),
                      (Path(job['source']), 'source_sha256')):
        assert sha(file) == result[key], str(file)
        evidence[str(file.resolve())] = result[key]
    evidence[str(result_file.resolve())] = sha(result_file)
    h = read(folder/'history.json')
    u = np.asarray([r['chemistry'] for r in h])
    v = np.asarray([r['volumes'] for r in h])
    delta = np.asarray([r['delta'] for r in h])
    assert len(h) == 401 and h[0]['elapsed'] == 0. and h[-1]['elapsed'] == 60.
    assert np.isfinite(u).all() and (u > 0).all() and (v > 0).all()
    assert np.max(abs(np.std(np.log(u[:, 0]), axis=1)-[r['log_activator_sd'] for r in h])) < 1e-12
    assert np.max(abs(delta.sum(axis=2))) < 1e-12
    assert np.max(abs(np.einsum('ti,tij->tj', v, delta))) < 1e-12
    histories[job['key']] = h
    completed.append(result)
    progress[job['key']] = 60.

reports = {}
for label in ('pilot', 'full'):
    rows = []
    for c in p['contexts']:
        file = ROOT/'pairs'/c['key']/('wait-'+label+'.json')
        if file.exists():
            r = read(file)
            assert r['passed'] and r['protocol_sha256'] == protocol_sha
            for asset, expected in r['evidence_sha256'].items():
                assert sha(asset) == expected, asset
            evidence[str(file.resolve())] = sha(file)
            rows.append(r)
    errors = {key: max(r['mechanical'][0]['errors'][key] for r in rows)
              for key in p['mechanical_limits']} if rows else {}
    reports[label] = {'passed': len(rows), 'expected': 12, 'maximum_errors': errors}

allocation_pairs = []
for pair in p['pairs']:
    jobs = {C: next(j for j in p['jobs'] if j['context'] == key and j['stage'] == 'wait' and j['level'] == 'fine')
            for C, key in pair['contexts'].items()}
    row = {'history': pair['seed'], 'cell': pair['cell'], 'both_fine_waits_complete': all(j['key'] in histories for j in jobs.values())}
    if row['both_fine_waits_complete']:
        hs = {C: histories[j['key']] for C, j in jobs.items()}
        i = hs['S'][0]['ids'].index(pair['cell'])
        initial = [np.asarray(h[0]['chemistry'])[:, i] for h in hs.values()]
        assert np.array_equal(*initial)
        final = {C: np.asarray(h[-1]['chemistry'])[:, i] for C, h in hs.items()}
        row.update(initial_recipient_pair=initial[0].tolist(),
                   final_recipient_pairs={C: x.tolist() for C, x in final.items()},
                   final_R_over_S=(final['R']/final['S']).tolist(),
                   final_activator_fold_separation=float(max(final['R'][0]/final['S'][0], final['S'][0]/final['R'][0])))
        budgets = {C: (np.asarray(h[0]['chemistry'])@np.asarray(h[0]['volumes'])).tolist() for C, h in hs.items()}
        assert np.max(abs(np.asarray(budgets['S'])/budgets['R']-1)) < 2e-13
        row['initial_whole_species_amounts'] = budgets
        row['final_whole_species_amounts'] = {C: (np.asarray(h[-1]['chemistry'])@np.asarray(h[-1]['volumes'])).tolist() for C, h in hs.items()}
        row['total_conductance_fractional_change'] = {}
        for C, h in hs.items():
            sums = []
            for observation in (h[0], h[-1]):
                g = np.asarray(observation['volumes'])[:, None]*np.asarray(observation['delta'])
                np.fill_diagonal(g, 0.)
                sums.append(float(g.sum()/2))
            row['total_conductance_fractional_change'][C] = sums[1]/sums[0]-1
        row['final_R_over_S_axis_ratio'] = hs['R'][-1]['axis_ratio']/hs['S'][-1]['axis_ratio']
    allocation_pairs.append(row)

rates = {}
for level in ('coarse', 'fine'):
    sample = [r['audit']['wall_seconds']/60 for r in completed if r['job']['level'] == level]
    assert sample
    rates[level] = {'samples': len(sample), 'median_seconds_per_unit': statistics.median(sample),
                    'range_seconds_per_unit': [min(sample), max(sample)]}
remaining = sum(max(0., j['duration']-progress[j['key']])*rates[j['level']]['median_seconds_per_unit'] for j in p['jobs'])
qualification = read(ROOT/'qualification.json')
prefix_pending = sum(not (ROOT/j['key']/'prefix/result.json').exists() for j in p['jobs'])
prefix_seconds = prefix_pending*qualification['seconds']/24
launch = read(Path('archive/study-launches/2026-10-08-moving-recipient-amount-response/launch-verification.json'))
now = datetime.now(timezone.utc)
wall_observed = sum(r['audit']['wall_seconds'] for r in completed)
for job in p['jobs']:
    if job['key'] not in histories and (ROOT/job['key']/'status.json').exists():
        wall_observed += read(ROOT/job['key']/'status.json').get('audit', {}).get('wall_seconds', 0.)
elapsed = (now-datetime.fromisoformat(launch['utc'])).total_seconds()
overhead_factor = max(1., elapsed/max(wall_observed, 1.))
central = remaining*overhead_factor+prefix_seconds
low, high = central*.9, central*1.15
audits = [r['audit'] for r in completed]
report = dict(utc=now.isoformat(), local=now.astimezone(ZoneInfo('America/New_York')).isoformat(),
              protocol_sha256=protocol_sha, root_status_at_start=snapshot,
              completed_paths=len(completed), total_paths=96,
              completed_wait_paths=sum(r['job']['stage']=='wait' for r in completed), expected_wait_paths=24,
              completed_response_paths=sum(r['job']['stage']=='response' for r in completed), expected_response_paths=72,
              histories=3, new_histories=0, paired_wait_checks=reports,
              full_window_contrast_retained=sum(r['whole_window_retained'] for r in completed),
              final_log_activator_sd_range=[min(r['final_log_sd'] for r in completed), max(r['final_log_sd'] for r in completed)],
              quality=dict(max_volume_error=max(a['max_volume_error'] for a in audits),
                           minimum_radius_grid_intervals=min(a['min_radius'] for a in audits),
                           max_boundary_occupancy=max(a['boundary_max'] for a in audits),
                           max_dilution_amount_error=max(a['dilution_error_max'] for a in audits),
                           max_clipping=max(a['max_clipping'] for a in audits)),
              paired_allocation_state_descriptions=allocation_pairs,
              timing=dict(elapsed_full_run_hours=elapsed/3600, measured_path_rates=rates,
                          unqualified_starts=prefix_pending, estimated_prefix_hours=prefix_seconds/3600,
                          observed_overhead_factor=overhead_factor, remaining_hours=central/3600,
                          remaining_hours_range=[low/3600, high/3600],
                          predicted_completion_local=(now+timedelta(seconds=central)).astimezone(ZoneInfo('America/New_York')).isoformat(),
                          assumptions='Response paths run at observed wait-path speed; prefix checks scale from 24 completed contexts. Range is a planning allowance, not a confidence interval. Gates may stop progression; final independent scientific review is additional.'),
              free_storage_gib=shutil.disk_usage(ROOT).free/2**30,
              evidence_sha256=evidence, script_sha256=sha(__file__),
              scope='Interim saved-data review, not scientific acceptance. State separation does not establish response separation, autonomous identity, inheritance or new development.')
file = HERE/('interim-'+now.strftime('%Y%m%dT%H%M%SZ')+'.json')
file.write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps({k: report[k] for k in ('utc', 'completed_paths', 'completed_wait_paths', 'completed_response_paths',
                                        'paired_wait_checks', 'full_window_contrast_retained', 'final_log_activator_sd_range',
                                        'quality', 'paired_allocation_state_descriptions', 'timing', 'free_storage_gib')}, indent=2))
print('report', file)
