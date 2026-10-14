"""Independent raw-history response reconstruction; no project helpers imported.

Read-only with respect to the workspace. Calculations precede comparison with
saved reports, and all scientific read paths receive before/after hashes.
"""
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path('/home/sji/Documents/Programming/Python/embryogenesis-simulation')
STUDY = ROOT/'outputs/phase-carry-redistribution-response'
OUTPUT = Path('/tmp/redistribution-response-independent-reconstruction.json')
PLOT_OUTPUT = Path('/tmp/redistribution-response-independent-plot-data.json')
HASHES = {}

def sha(path):
    path = Path(path)
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for data in iter(lambda: handle.read(2**20), b''):
            h.update(data)
    return h.hexdigest()

def checked_read(path):
    path = Path(path).resolve()
    content = path.read_bytes()
    expected = hashlib.sha256(content).hexdigest()
    if str(path) in HASHES:
        assert HASHES[str(path)] == expected, str(path)
    HASHES[str(path)] = expected
    return json.loads(content)

def record_text(path):
    path = Path(path).resolve()
    HASHES[str(path)] = sha(path)
    return path.read_text()

def integral(y, t):
    y, t = np.asarray(y), np.asarray(t)
    return float(np.sum(np.diff(t)*(y[1:]+y[:-1])/2))

def rms(first, second, t):
    return float(np.sqrt(integral(np.mean((first-second)**2, axis=1), t)/(t[-1]-t[0])))

def recover(t, deviation, threshold):
    # Independently scan every candidate and its observed suffix.
    for k, value in enumerate(t):
        if value <= t[-1]-24. and np.all(deviation[k:] <= threshold):
            return float(value)
    return None

def close(actual, expected, label):
    if actual is None or expected is None:
        assert actual is expected, label
    elif isinstance(actual, bool):
        assert actual == expected, label
    else:
        assert np.isclose(actual, expected, rtol=5e-11, atol=5e-12), (label, actual, expected)

p = checked_read(STUDY/'protocol.json')
status = checked_read(STUDY/'status.json')
assert status['state'] == 'completed' and status['passed'] and status['completed'] == 66
record_text(ROOT/'docs/phase_carry_redistribution_response.md')
for name in ('phase_carry_redistribution_response.py', 'phase_carry_delayed_response.py',
             'cell_response.py', 'attribute_persistence.py', 'cell_response_exchange.py',
             'cell_response_moving_refinement.py', 'cell_exchange_response_moving.py'):
    record_text(ROOT/'embryo'/name)

jobs = {(j['seed'], j['level'], j['background'], j['target'], j['dose']): j for j in p['jobs']}
assert len(jobs) == 66 and len(p['reused']) == 18
raw = {}
initial_dose_rows = []
max_dose_error = 0.
for key, job in jobs.items():
    reused = p['reused'].get(job['key'])
    physical = (Path(p['reuse_parent'])/reused['job']['key']) if reused else STUDY/job['key']
    history = checked_read(physical/'history.json')
    result = checked_read(physical/'result.json')
    assert result['quality_pass'] and result['history_sha256'] == HASHES[str((physical/'history.json').resolve())]
    t = np.asarray([r['elapsed'] for r in history])
    assert len(t) == 401 and np.allclose(t, np.arange(401)*.15, rtol=0, atol=1e-9)
    assert np.allclose([r['time'] for r in history], t+510., rtol=0, atol=1e-9)
    ids = history[0]['ids']
    assert len(ids) == len(set(ids)) == 16 and all(r['ids'] == ids for r in history)
    values, volumes = np.array([r['chemistry'] for r in history]), np.array([r['volumes'] for r in history])
    assert values.shape == (401, 2, 16) and volumes.shape == (401, 16)
    assert np.isfinite(values).all() and np.all(values > 0) and np.isfinite(volumes).all() and np.all(volumes > 0)
    # Independently ensure raw trajectory equals the paired saved observations.
    for label, number in (('pilot', 41), ('full', 401)):
        paired = checked_read(STUDY/'pairs'/f'seed-{job["seed"]}'/f'{label}_{job["key"]}.json')
        assert paired == history[:number]
    # Source fields are partially read, never copied or rewritten.
    with np.load(job['source']) as source:
        assert np.array_equal(source['ids'], ids)
        assert np.array_equal(values[0], [source['activator'], source['inhibitor']])
        assert source['activator'].dtype == np.dtype('float64')
        assert source['inhibitor'].dtype == np.dtype('float64')
        metadata = json.loads(str(source['metadata']))
        assert metadata['time'] == 510. and metadata['config']['dt'] == job['dt']
    raw[key] = dict(times=t, ids=ids, values=values, volumes=volumes)

responses, factors = {}, {}
for key, job in jobs.items():
    if job['target'] is None:
        continue
    seed, level, background, target, dose = key
    own = raw[seed, level, background, None, None]
    trial = raw[key]
    i = trial['ids'].index(target)
    assert own['ids'] == trial['ids'] and np.array_equal(own['volumes'][0], trial['volumes'][0])
    factor = float(trial['values'][0, 0, i]/own['values'][0, 0, i])
    expected = own['values'][0].copy()
    expected[0, i] *= job['factor']
    assert np.array_equal(trial['values'][0], expected), job['key']
    close(factor, job['factor'], 'observed factor '+job['key'])
    amplitude = float(abs(np.log(factor)))
    close(amplitude, job['initial_log_amplitude'], 'actual initial log amplitude')
    mass = own['volumes'][0, i]
    target_amount = float(own['values'][0, 0, i]*mass)
    removed = float((own['values'][0, 0, i]-trial['values'][0, 0, i])*mass)
    sham = raw[seed, level, 'sham', None, None]
    red = raw[seed, level, f'redistributed-cell-{target}', None, None]
    amounts = [float(bg['values'][0, 0, i]*bg['volumes'][0, i]) for bg in (sham, red)]
    requested = .1*target_amount if dose == 'fractional' else .1*min(amounts)
    error = abs(removed/requested-1.)
    assert error <= 2e-13, (job['key'], error)
    max_dose_error = max(max_dose_error, error)
    close(removed, job['removed_activator_amount'], 'saved removed amount')
    close(requested, job['requested_removed_amount'], 'saved requested dose')
    close(target_amount, job['initial_target_activator_amount'], 'saved target amount')
    weights = own['volumes'][0]
    log_change = np.log(trial['values']/own['values'])
    # Network distance historically uses a difference of logs; keep this exact
    # mathematical definition separate from target ratios for recovery.
    network = np.sqrt(np.sum((np.log(trial['values'])-np.log(own['values']))**2*weights, axis=(1, 2))/(2*weights.sum()))
    target_change = np.max(abs(log_change[:, :, i]), axis=1)
    other = np.arange(16) != i
    other_rms = np.sqrt(np.sum(log_change[:, :, other]**2*weights[other], axis=(1, 2))/(2*weights[other].sum()))
    metrics = dict(target_activator_peak_gain=float(abs(log_change[:, 0, i]).max()/amplitude),
        target_inhibitor_peak_gain=float(abs(log_change[:, 1, i]).max()/amplitude),
        target_activator_log_auc_per_log_pulse=integral(abs(log_change[:, 0, i]), trial['times'])/amplitude,
        other_cells_peak_log_rms_per_log_pulse=float(other_rms.max()/amplitude),
        target_recovery_time=recover(trial['times'], target_change, .1*amplitude),
        network_recovery_time=recover(trial['times'], network, .1*network[0]),
        final_log_rms_from_control=float(network[-1]),
        late_max_log_rms_from_control=float(network[trial['times'] >= 36.].max()))
    responses[key] = dict(wave=log_change[:, :, i]/amplitude, log_change=log_change,
        amplitude=amplitude, metrics=metrics)
    factors[key] = factor
    initial_dose_rows.append(dict(job=job['key'], history=seed, level=level,
        background=background, recipient=target, dose=dose, observed_factor=factor,
        actual_log_amplitude=amplitude, measured_removed_amount=removed,
        requested_removed_amount=requested, relative_dose_error=error,
        target_starting_amount=target_amount, starting_volume=float(mass),
        starting_global_species_amounts=(own['values'][0]@weights).tolist(), metrics=metrics))

recipients, comparisons, history_counts, plot_rows = [], [], [], []
all_gate_passed = True
effect_discrepancy = saved_effect_error = normalized_error_max = auc_error_max = recovery_error_max = 0.
for seed in p['histories']:
    report = checked_read(STUDY/'pairs'/f'seed-{seed}'/'full-refinement.json')
    assert report['horizon'] == 60. and report['protocol_sha256'] == HASHES[str((STUDY/'protocol.json').resolve())]
    counts = {d: 0 for d in ('fractional', 'amount')}
    both = 0
    for target in p['selections'][str(seed)]['selected_ids']:
        item = dict(history=seed, cell=target, by_timestep={})
        for level in ('coarse', 'fine'):
            saved = next(r for r in report['responses'][level] if r['cell'] == target)
            effects = {}
            for dose in ('fractional', 'amount'):
                a = responses[seed, level, 'sham', target, dose]
                b = responses[seed, level, f'redistributed-cell-{target}', target, dose]
                t = raw[seed, level, 'sham', target, dose]['times']
                effect = rms(b['wave'], a['wave'], t)
                raw_effect = rms(b['log_change'][:, :, raw[seed, level, 'sham', target, dose]['ids'].index(target)],
                                 a['log_change'][:, :, raw[seed, level, 'sham', target, dose]['ids'].index(target)], t)
                detected = bool(effect > p['effect_min'])
                original = saved['differences'][dose]
                close(effect, original['normalized_response_shift_rms'], 'saved primary effect')
                close(raw_effect, original['raw_log_response_shift_rms'], 'saved raw effect')
                close(detected, original['response_shift_detected'], 'saved primary decision')
                saved_effect_error = max(saved_effect_error, abs(effect-original['normalized_response_shift_rms']))
                effects[dose] = dict(normalized_response_shift_rms=effect,
                    raw_log_response_shift_rms=raw_effect, response_shift_detected=detected,
                    signed_margin_above_screen=effect-p['effect_min'])
                for background, response in (('sham', a), (f'redistributed-cell-{target}', b)):
                    saved_metrics = saved['contexts'][background]['doses'][dose]['metrics']
                    for metric, value in response['metrics'].items():
                        close(value, saved_metrics[metric], 'saved '+metric)
            for background in ('sham', f'redistributed-cell-{target}'):
                waves = [responses[seed, level, background, target, d]['wave'] for d in ('fractional', 'amount')]
                close(rms(*waves, t), saved['contexts'][background]['dose_dependence_normalized_rms'], 'within-background dose dependence')
            item['by_timestep'][level] = effects
        agreement = all(item['by_timestep']['coarse'][d]['response_shift_detected'] == item['by_timestep']['fine'][d]['response_shift_detected'] for d in counts)
        assert agreement
        item['same_effect_decisions'] = agreement
        item['shifted_at_both_doses'] = all(item['by_timestep']['fine'][d]['response_shift_detected'] for d in counts)
        both += int(item['shifted_at_both_doses'])
        for dose in counts:
            counts[dose] += int(item['by_timestep']['fine'][dose]['response_shift_detected'])
            discrepancy = abs(item['by_timestep']['coarse'][dose]['normalized_response_shift_rms']-item['by_timestep']['fine'][dose]['normalized_response_shift_rms'])
            effect_discrepancy = max(effect_discrepancy, discrepancy)
            assert discrepancy <= .01
            for background in ('sham', f'redistributed-cell-{target}'):
                a, b = [responses[seed, level, background, target, dose] for level in ('coarse', 'fine')]
                close(factors[seed, 'coarse', background, target, dose], factors[seed, 'fine', background, target, dose], 'paired factors')
                error = float(np.max(abs(a['log_change']-b['log_change']))/a['amplitude'])
                metric = 'target_activator_log_auc_per_log_pulse'
                auc_error = abs(b['metrics'][metric]/a['metrics'][metric]-1.)
                recovery_errors = {}
                for metric in ('target_recovery_time', 'network_recovery_time'):
                    x, y = [r['metrics'][metric] for r in (a, b)]
                    recovery_errors[metric] = 0. if x is None and y is None else None if (x is None or y is None) else abs(x-y)
                lim = p['response_limits']
                passed = (error <= lim['max_normalized_response_error'] and auc_error <= lim['max_relative_auc_error'] and all(x is not None and x <= lim['max_recovery_time_error']+1e-9 for x in recovery_errors.values()))
                row = next(r for r in report['responses']['rows'] if (r['context'], r['cell'], r['dose']) == (background, target, dose))
                close(error, row['normalized_response_error'], 'saved refinement waveform error')
                close(auc_error, row['relative_auc_error'], 'saved AUC error')
                for metric, value in recovery_errors.items(): close(value, row['recovery_time_errors'][metric], 'saved recovery difference')
                close(bool(passed), row['passed'], 'saved response refinement decision')
                all_gate_passed &= bool(passed)
                normalized_error_max = max(normalized_error_max, error)
                auc_error_max = max(auc_error_max, auc_error)
                recovery_error_max = max(recovery_error_max, *[x or 0. for x in recovery_errors.values()])
                comparisons.append(dict(history=seed, cell=target, background=background, dose=dose,
                    normalized_response_error=error, relative_auc_error=auc_error,
                    recovery_time_errors=recovery_errors, passed=bool(passed)))
        sham, red = [raw[seed, 'fine', background, None, None] for background in ('sham', f'redistributed-cell-{target}')]
        i = sham['ids'].index(target)
        late = sham['times'] >= 36.
        late_log_difference = np.log(red['values'][late, :, i]/sham['values'][late, :, i])
        baseline_late_rms = float(np.sqrt(integral(np.mean(late_log_difference**2, axis=1), sham['times'][late])/24.))
        item['fine_baseline'] = dict(initial_species_ratios=(red['values'][0, :, i]/sham['values'][0, :, i]).tolist(),
            final_species_ratios=(red['values'][-1, :, i]/sham['values'][-1, :, i]).tolist(),
            late_two_species_state_log_rms=baseline_late_rms,
            initial_target_species_amounts={name:(bg['values'][0, :, i]*bg['volumes'][0, i]).tolist() for name, bg in (('sham', sham), ('redistributed', red))},
            initial_aggregate_species_amounts={name:(bg['values'][0]@bg['volumes'][0]).tolist() for name, bg in (('sham', sham), ('redistributed', red))},
            initial_aggregate_amount_ratios=((red['values'][0]@red['volumes'][0])/(sham['values'][0]@sham['volumes'][0])).tolist(),
            final_aggregate_amount_ratios=((red['values'][-1]@red['volumes'][-1])/(sham['values'][-1]@sham['volumes'][-1])).tolist())
        plot_rows.append(dict(history=seed, cell=target, fine_baseline=item['fine_baseline'],
            times=sham['times'].tolist(), fine_baseline_species_ratios=(red['values'][:, :, i]/sham['values'][:, :, i]).tolist(),
            doses={dose:{name:dict(normalized_two_species_response=responses[seed, 'fine', background, target, dose]['wave'].tolist(),
                observed_factor=factors[seed, 'fine', background, target, dose],
                actual_log_amplitude=responses[seed, 'fine', background, target, dose]['amplitude'],
                recovery_metrics=responses[seed, 'fine', background, target, dose]['metrics'])
                for name, background in (('sham', 'sham'), ('redistributed', f'redistributed-cell-{target}'))}
                for dose in ('fractional', 'amount')}))
        recipients.append(item)
    history_counts.append(dict(history=seed, recipients=2, fractional_above_screen=counts['fractional'],
        amount_above_screen=counts['amount'], both_doses_above_screen=both))

for path, expected in HASHES.items():
    assert sha(path) == expected, ('Read evidence changed during reconstruction', path)

result = dict(created_utc=datetime.now(timezone.utc).isoformat(),
    independent_response_reconstruction_passed=all_gate_passed,
    scope='Response/dose/raw-history audit only, not full mechanical/provenance acceptance. No project numerical helpers called. Histories are the three existing independent units; recipients, doses and timestep paths are nested.',
    root=str(STUDY), protocol_sha256=HASHES[str((STUDY/'protocol.json').resolve())],
    script_sha256=sha(__file__), independent_histories=3, new_histories=0,
    raw_logical_paths_checked=len(raw), normalized_target_effects=24, pulse_metric_records=48,
    full_response_refinement_rows=len(comparisons), effect_screen=p['effect_min'],
    max_relative_pulse_amount_error=max_dose_error,
    max_difference_from_saved_primary_effect=saved_effect_error,
    max_coarse_fine_primary_effect_discrepancy=effect_discrepancy,
    max_full_network_normalized_response_error=normalized_error_max,
    max_relative_target_auc_error=auc_error_max,
    max_recovery_time_discrepancy=recovery_error_max,
    all_48_target_recovered=all(r['metrics']['target_recovery_time'] is not None for r in initial_dose_rows),
    all_48_network_recovered=all(r['metrics']['network_recovery_time'] is not None for r in initial_dose_rows),
    history_outcomes=history_counts, recipients=recipients, pulse_doses_and_metrics=initial_dose_rows,
    response_refinement=comparisons, evidence_sha256=HASHES,
    source_npz_scope='Prepared source chemistry/IDs/time/dtype/config checked by partial ZIP reads. Full source and checkpoint hashes and physical field/carry equality require the separate provenance/mechanics audit.')
PLOT_OUTPUT.write_text(json.dumps(dict(scope='Fine dt=0.001875; two species column order activator, inhibitor; target responses relative to each own control, normalized by actual log pulse amplitude; baseline ratios redistribution/sham; three existing histories.',
    protocol_sha256=result['protocol_sha256'], script_sha256=result['script_sha256'],
    data=plot_rows, evidence_sha256=HASHES), separators=(',', ':'))+'\n')
result['plot_data'] = str(PLOT_OUTPUT)
result['plot_data_sha256'] = sha(PLOT_OUTPUT)
OUTPUT.write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ('pulse_doses_and_metrics', 'response_refinement', 'evidence_sha256')}, indent=2))
print('Saved independent reconstruction:', OUTPUT)
