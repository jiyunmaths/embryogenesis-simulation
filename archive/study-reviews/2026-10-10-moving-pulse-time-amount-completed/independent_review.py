"""Reconstruct completed moving evidence without project helpers or stepping.

This checks saved arrays and committed runtime evidence. It does not recreate
unsaved GPU restart branches or establish full-horizon native equivalence.
"""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import math
import numpy as np

ROOT = Path('output/studies/moving-pulse-time-amount-response')
HERE = Path(__file__).parent
PINS = {}


def digest(file):
    h = hashlib.sha256()
    with Path(file).open('rb') as handle:
        for block in iter(lambda: handle.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def pin(file, expected=None):
    name = str(Path(file).resolve())
    value = PINS.get(name)
    if value is None:
        value = digest(name)
        PINS[name] = value
    assert expected is None or value == expected, name
    return Path(name)


def read(file):
    return json.loads(pin(file).read_text())


def arrays(file, keys):
    with np.load(pin(file), allow_pickle=False) as z:
        return {key: z[key].copy() for key in keys}


def close(actual, expected):
    if isinstance(actual, dict):
        for key, value in actual.items():
            close(value, expected[key])
    elif isinstance(actual, list):
        assert len(actual) == len(expected)
        for a, b in zip(actual, expected):
            close(a, b)
    elif actual is None or isinstance(actual, (bool, str)):
        assert actual == expected, (actual, expected)
    else:
        assert math.isclose(actual, expected, rel_tol=1e-11, abs_tol=1e-12), (actual, expected)


def errors(first, reference):
    u, v = [np.asarray([x['chemistry'] for x in h]) for h in (first, reference)]
    a, b = [np.asarray([x['delta'] for x in h]) for h in (first, reference)]
    return dict(
        chemical_log_max=float(abs(np.log(u/v)).max()),
        polarity_abs_max=float(abs(np.asarray([x['polarity'] for x in first])-np.asarray([x['polarity'] for x in reference])).max()),
        relative_axis_max=float(abs(np.asarray([x['axis_ratio'] for x in first])/np.asarray([x['axis_ratio'] for x in reference])-1).max()),
        relative_volume_max=float(abs(np.asarray([x['volumes'] for x in first])/np.asarray([x['volumes'] for x in reference])-1).max()),
        relative_transport_max=float((np.linalg.norm(a-b, axis=(1,2))/np.maximum(np.linalg.norm(b, axis=(1,2)), 1e-300)).max()))


def check_history(h, job, horizon, source):
    times = np.arange(round(horizon/.15)+1)*.15
    assert len(h) == len(times)
    np.testing.assert_allclose([r['elapsed'] for r in h], times, rtol=0, atol=1e-9)
    np.testing.assert_array_equal(h[0]['chemistry'], source['initial'])
    u = np.asarray([r['chemistry'] for r in h])
    v = np.asarray([r['volumes'] for r in h])
    d = np.asarray([r['delta'] for r in h])
    assert np.isfinite(u).all() and (u > 0).all() and np.isfinite(v).all() and (v > 0).all()
    assert np.isfinite(d).all()
    g = v[:,:,None]*d
    assert abs(d.sum(axis=2)).max() < 1e-11
    assert abs(np.einsum('ti,tij->tj',v,d)).max() < 1e-12
    assert abs(g-g.transpose(0,2,1)).max() < 1e-12
    off = d.copy()
    off[:,np.arange(16),np.arange(16)] = 0
    assert off.min() >= 0 and np.diagonal(d,axis1=1,axis2=2).max() <= 0
    sd = np.std(np.log(u[:,0]), axis=1)
    np.testing.assert_allclose(sd,[r['log_activator_sd'] for r in h],rtol=0,atol=1e-12)
    for r,t in zip(h,times):
        np.testing.assert_array_equal(r['ids'],source['ids'])
        assert abs(r['time']-job['start']-t) < 1e-9
        assert r['boundary_occupancy'] < .01
        assert np.isfinite(r['polarity']).all() and np.isfinite(r['centers']).all() and r['axis_ratio'] > 0
    return u,v,d,sd


def contrast(h, horizon):
    t = np.asarray([r['elapsed'] for r in h])
    sd = np.std(np.log(np.asarray([r['chemistry'][0] for r in h])),axis=1)
    late = t >= horizon-min(24.,horizon)-1e-9
    assert sd.min() > .1  # The completed dataset never loses contrast.
    return dict(final_log_sd=float(sd[-1]),late_min_log_sd=float(sd[late].min()),persistent_contrast=True,
                whole_window_min_log_sd=float(sd.min()),whole_window_retained=True,first_loss=None)


def recovery(t, y, threshold):
    suffix = np.logical_and.accumulate((y <= threshold)[::-1])[::-1]
    indices = np.flatnonzero(suffix & (t <= t[-1]-24.))
    return float(t[indices[0]]) if len(indices) else None


def response(h, ctrl, cell, factor):
    u,b = [np.asarray([r['chemistry'] for r in x]) for x in (h,ctrl)]
    t = np.asarray([r['elapsed'] for r in h]);v = np.asarray(h[0]['volumes'])
    i = h[0]['ids'].index(cell);amp = abs(np.log(factor));lc = np.log(u/b)
    net = np.sqrt(np.sum((np.log(u)-np.log(b))**2*v,axis=(1,2))/(2*v.sum()))
    others = np.arange(16) != i
    other = np.sqrt(np.sum(lc[:,:,others]**2*v[others],axis=(1,2))/(2*v[others].sum()))
    metrics = dict(target_activator_peak_gain=float(abs(lc[:,0,i]).max()/amp),
        target_inhibitor_peak_gain=float(abs(lc[:,1,i]).max()/amp),
        target_activator_log_auc_per_log_pulse=float(np.trapezoid(abs(lc[:,0,i]),t)/amp),
        other_cells_peak_log_rms_per_log_pulse=float(other.max()/amp),
        target_recovery_time=recovery(t,abs(lc[:,:,i]).max(axis=1),.1*amp),
        network_recovery_time=recovery(t,net,.1*net[0]),
        final_log_rms_from_control=float(net[-1]),late_max_log_rms_from_control=float(net[t>=t[-1]-24.].max()))
    return lc/amp, lc[:,:,i]/amp, metrics


def main():
    p = read(ROOT/'protocol.json'); ph = digest(ROOT/'protocol.json')
    summary = read(ROOT/'summary.json'); status = read(ROOT/'status.json'); timing = read(ROOT/'timing.json')
    assert status['state'] == 'completed' and status['passed'] and summary['passed']
    assert summary['protocol_sha256'] == ph and summary['total_paths'] == 48
    assert p['histories'] == [7,8,9] and p['new_histories'] == 0
    assert len(p['jobs']) == 48 and len(p['contexts']) == 12 and len(p['pairs']) == 6
    assert p['primary_reference'] == p['budget_reference'] == p['initial_physics'] == 'S'
    assert p['start'] == 570. and p['wait'] == 0. and p['response_duration'] == 60.
    for collection in ('source_sha256','input_sha256'):
        for file,expected in p[collection].items(): pin(file,expected)
    for file,expected in summary['evidence_sha256'].items(): pin(file,expected)
    pre = read(p['preparation_audit']); qa = read(p['qualification_audit'])
    q = read(ROOT/'qualification.json'); pilot = read(ROOT/'pilot-gate.json')
    assert all(x['passed'] and x['protocol_sha256'] == ph for x in (pre,qa,q,pilot))
    assert pre['total_paths'] == qa['actual_start_gates'] == q['jobs'] == 48
    assert pilot['contexts'] == 12 and pilot['contrast_pairs'] == 6
    for file,expected in pre['inputs_sha256'].items(): pin(file,expected)
    for file,expected in pilot['evidence_sha256'].items(): pin(file,expected)

    # The user's introduction rewrite is an authorized non-input editorial
    # exception, with original bytes preserved. Do not refresh experiment pins.
    guards = read('archive/study-document-snapshots/2026-10-09-before-moving-pulse-time-amount-response/manifest.json')
    editorial = read('archive/study-document-snapshots/2026-10-09-authorized-manuscript-introduction-edit/manifest.json')
    changed = {x['path']:x for x in editorial['changed_noninput_manuscript_files']}
    exceptions = []
    for file,expected in guards['preserved_ledger_manuscript_sha256'].items():
        actual = digest(file)
        if actual != expected:
            note = changed[file]; assert note['before'] == expected and note['after'] == actual
            assert str(Path(file).resolve()) not in p['source_sha256'] and str(Path(file).resolve()) not in p['input_sha256']
            original = Path(editorial['prior_manuscript_archive'])/file
            pin(original,expected); pin(file,actual)
            exceptions.append(dict(path=file,before=expected,after=actual,original=str(original),authorization=editorial['authorization']))
        else: pin(file,expected)
    assert len(exceptions) == 5
    before = read('archive/study-document-snapshots/2026-10-10-before-completed-moving-pulse-time-amount-review/manifest.json')
    for file,expected in before['preserved_ledger_manuscript_sha256'].items(): pin(file,expected)

    histories,results,sources,prefix_rows = {},[],{},[]
    def job_for(context,level,dose=None):
        return next(j for j in p['jobs'] if (j['context'],j['level'],j['dose']) == (context,level,dose))
    for number,job in enumerate(p['jobs'],1):
        folder = ROOT/job['key']; r = read(folder/'result.json'); h = read(folder/'history.json')
        assert r['job'] == job and r['protocol_sha256'] == ph and r['method'] == p['method'] and r['quality_pass']
        for file,key in ((folder/'history.json','history_sha256'),(folder/'latest_state.npz','checkpoint_sha256'),(job['source'],'source_sha256')): pin(file,r[key])
        source = arrays(job['source'],('initial','correction','ids','volumes','delta')); sources[job['key']] = source
        assert np.all(source['correction'] == 0) and source['correction'].dtype == np.float64
        u,v,d,sd = check_history(h,job,60.,source); histories[job['key']] = h
        cp = arrays(folder/'latest_state.npz',('metadata','activator','inhibitor','chemical_correction','chemical_method','phase_carry','rounding','phi','polarity','ids','parents','target','fate','long_experiment','precision_arm','attribute_mode'))
        meta = json.loads(str(cp['metadata'])); commit = json.loads(str(cp['long_experiment']))
        assert meta['config'] == p['accepted_configs'][job['context']+'_'+job['level']]
        assert abs(meta['time']-630.) < 1e-9 and meta['step_number'] == round(meta['time']/job['dt'])
        assert not meta['divisions'] and commit['protocol_hash'] == ph and commit['job'] == job
        assert commit['history'] == h and commit['audit'] == r['audit']
        np.testing.assert_array_equal(u[-1],[cp['activator'],cp['inhibitor']]); np.testing.assert_array_equal(h[-1]['polarity'],cp['polarity'])
        np.testing.assert_array_equal(source['ids'],cp['ids'])
        e = cp['chemical_correction']; last = np.asarray([cp['activator'],cp['inhibitor']])
        assert e.dtype == np.float64 and e.shape == last.shape and np.isfinite(e).all() and np.all(abs(e)<=2*np.spacing(last))
        assert str(cp['chemical_method']) == p['method'] and str(cp['precision_arm']) == 'phase_carry' and str(cp['attribute_mode']) == 'direct'
        assert cp['phi'].dtype == np.float32 and cp['phase_carry'].dtype == np.float64 and cp['phi'].shape == cp['phase_carry'].shape == (16,72,72,72)
        assert np.isfinite(cp['phase_carry']).all() and np.isfinite(cp['phi']).all() and cp['phi'].min() >= 0 and cp['phi'].max() <= 1
        assert cp['rounding'].dtype == np.int32 and cp['rounding'].shape == (2,) and np.all(cp['rounding'] >= 0) and np.all(cp['fate'] == 0)
        audit = r['audit']; assert audit['max_volume_error'] < .05 and audit['min_radius'] >= 4 and audit['max_clipping'] == 0
        assert audit['boundary_max'] < .01 and audit['dilution_error_max'] <= p['amount_error_max']
        assert abs(v/cp['target']-1).max() <= audit['max_volume_error']+1e-12
        assert max(x['boundary_occupancy'] for x in h) <= audit['boundary_max']+1e-16
        close(contrast(h,60.),r); results.append(r); del cp

        pref = folder/'prefix'; qr = read(pref/'result.json')
        assert qr['passed'] and qr['job'] == job and qr['protocol_sha256'] == ph and qr['method'] == p['method']
        pin(job['source'],qr['source_sha256'])
        for file,expected in qr['evidence_sha256'].items(): pin(file,expected)
        gh,ch = read(pref/'gpu-history.json'),read(pref/'cpu-history.json')
        for hh in (gh,ch): check_history(hh,job,.6,source)
        measured = errors(gh,ch); close(measured,qr['errors'])
        for key,value in measured.items(): assert value <= p['prefix_limits'][key]
        a,b = [arrays(pref/f'{arm}-phi.npz',('phi',))['phi'] for arm in ('gpu','cpu')]
        phi_error = float(abs(a-b).max()); close(phi_error,qr['final_phi_abs_max'])
        assert phi_error <= p['prefix_limits']['final_phi_abs_max']; del a,b
        assert qr['independent_chemical_steps'] == 8 and qr['restart_exact_steps'] == 4
        assert qr['cpu_chemical_log_max'] <= 1e-11 and qr['amount_accounting_max'] <= 2e-14
        cp = arrays(pref/'restart.npz',('metadata','activator','inhibitor','chemical_correction','chemical_method'))
        cu = np.asarray([cp['activator'],cp['inhibitor']]); ce = cp['chemical_correction']
        assert ce.dtype == np.float64 and ce.shape == cu.shape and np.isfinite(ce).all() and np.all(abs(ce)<=2*np.spacing(cu))
        assert str(cp['chemical_method']) == p['method']
        assert abs(json.loads(str(cp['metadata']))['time']-570.-8*job['dt']) < 1e-9
        for key in ('volumes','polarity','centers','delta'): np.testing.assert_array_equal(gh[0][key],h[0][key])
        prefix_rows.append(dict(job=job['key'],errors=measured,final_phi_abs_max=phi_error,runtime_eight_step_chemical_error=qr['cpu_chemical_log_max'],runtime_dilution_amount_error=qr['amount_accounting_max'],passed=True))
        if number%12 == 0: print(f'Verified {number}/48 completed paths and saved prefix arrays',flush=True)

    parent = read(Path(p['parent'])/'protocol.json'); states = {}; preparations = []
    for c in p['contexts']:
        old = next(x for x in parent['contexts'] if x['key'] == c['key'])
        assert all(c[k] == old[k] for k in old if k != 'source') and c['accepted_context_source'] == old['source']
        z = arrays(c['source'],('initial','delta','volumes','ids')); oldz = arrays(c['accepted_context_source'],tuple(z)); states[c['key']] = z
        for key in z: np.testing.assert_array_equal(z[key],oldz[key])
        physical = Path(c['physical_source']); physical_h = read(physical.with_name('history.json'))
        assert physical_h[-1]['time'] == 570.
        for key in ('volumes','delta','ids'): np.testing.assert_array_equal(z[key],physical_h[-1][key])
        cp = arrays(physical,('polarity','ids','parents','target','metadata','rounding','phase_carry'))
        assert json.loads(str(cp['metadata']))['time'] == 570. and np.isfinite(cp['phase_carry']).all()
        for dose in (None,'fractional'):
            job = job_for(c['key'],'fine',dose); file = Path(job['source']); m = read(file.with_suffix('.json'))
            pin(file,m['paths_sha256']); pin(physical,m['physical_source_sha256']); pin(m['parent_source'],m['parent_source_sha256'])
            assert m['physical_source'] == str(physical)
            trial = next(t for t in parent['trials'] if t['context'] == c['key'] and t['target'] == job['target'])
            assert m['parent_trial'] == trial['key'] and m['parent_source'] == trial['source']
            s = sources[job['key']]; expected = z['initial'].copy(); i = z['ids'].tolist().index(c['recipient'])
            if dose: expected[0,i] *= .9
            np.testing.assert_array_equal(s['initial'],expected); np.testing.assert_array_equal(s['initial'],arrays(m['parent_source'],('initial',))['initial'])
            close(1. if dose is None else .9,m['record']['factor'])
            removed = float((z['initial'][0,i]-s['initial'][0,i])*z['volumes'][i]); close(removed,m['record']['removed_activator_amount'])
            for level in ('coarse','fine'):
                h = histories[job_for(c['key'],level,dose)['key']]
                for key in ('volumes','delta','ids','polarity','centers','axis_ratio'): np.testing.assert_array_equal(h[0][key],physical_h[-1][key])
                lastcp = arrays(ROOT/job_for(c['key'],level,dose)['key']/'latest_state.npz',('parents','target'))
                for key in lastcp: np.testing.assert_array_equal(lastcp[key],cp[key])
            preparations.append(dict(context=c['key'],dose=dose,exact_accepted_initial_preserved=True,explicit_zero_correction_verified=True,removed_activator_amount=removed))
        del cp

    mechanical_rows,response_rows = [],[]
    for c in p['contexts']:
        for label,horizon in (('pilot',6.),('full',60.)):
            record = read(ROOT/'pairs'/c['key']/f'response-{label}.json')
            assert record['passed'] and record['protocol_sha256'] == ph
            for file,expected in record['evidence_sha256'].items(): pin(file,expected)
            hs = {}
            for n,dose in enumerate((None,'fractional')):
                paths,fields = [],[]
                for level in ('coarse','fine'):
                    job = job_for(c['key'],level,dose); folder = ROOT/job['key']
                    h = read(folder/('pilot-history.json' if label == 'pilot' else 'history.json'))
                    check_history(h,job,horizon,sources[job['key']]); assert h == histories[job['key']][:len(h)]
                    hs[level,dose] = h; paths.append(h)
                    fields.append(arrays(folder/('pilot-phi.npz' if label == 'pilot' else 'latest_state.npz'),('phi',))['phi'])
                measured = errors(*paths); measured['final_phi_abs_max'] = float(abs(fields[0]-fields[1]).max())
                close(measured,record['mechanical'][n]['errors'])
                for key,value in measured.items(): assert value <= p['mechanical_limits'][key]
                close([contrast(h,horizon) for h in paths],record['mechanical'][n]['contrast_summaries'])
                assert record['mechanical'][n]['same_contrast_outcomes'] and record['mechanical'][n]['loss_times'] == [None,None]
                mechanical_rows.append(dict(context=c['key'],label=label,dose=dose,errors=measured,passed=True)); del fields
            data = [response(hs[level,'fractional'],hs[level,None],c['recipient'],.9) for level in ('coarse','fine')]
            normalized = float(abs(data[0][0]-data[1][0]).max())
            auc = abs(data[1][2]['target_activator_log_auc_per_log_pulse']/data[0][2]['target_activator_log_auc_per_log_pulse']-1)
            recover = {}
            for key in ('target_recovery_time','network_recovery_time'):
                a,b = [x[2][key] for x in data]
                recover[key] = 0. if a is None and b is None else None if a is None or b is None else abs(a-b)
            measured = dict(normalized_response_error=normalized,relative_auc_error=auc,recovery_time_errors=recover)
            close(measured,record['response'][0])
            assert normalized <= .01 and auc <= .02 and all(x is not None and x<=.3+1e-9 for x in recover.values())
            response_rows.append(dict(context=c['key'],label=label,**measured,passed=True))
        print('Reconstructed pilot/full gates:',c['key'],flush=True)

    comparisons,recoveries,allocations,archived = [],[],[],{}
    frozen_review = read('archive/study-reviews/2026-10-09-pulse-time-amount-completed/independent_review.json')
    assert frozen_review['passed'] and frozen_review['accepted_at_tested_scope']
    frozen_arrays = pin('archive/study-reviews/2026-10-09-pulse-time-amount-completed/independent_review.npz')
    frozen_comparisons = []
    for pair in p['pairs']:
        hs = {C:histories[job_for(context,'fine')['key']] for C,context in pair['contexts'].items()}
        i = hs['S'][0]['ids'].index(pair['cell']); v = np.asarray(hs['S'][0]['volumes']); others = np.arange(16) != i
        amounts,conductance,axis,initial_errors = {},{},{},[]
        for dose in (None,'fractional'):
            a,b = [sources[job_for(pair['contexts'][C],'fine',dose)['key']] for C in ('S','R')]
            for key in ('volumes','delta','ids'): np.testing.assert_array_equal(a[key],b[key])
            np.testing.assert_array_equal(a['initial'][:,i],b['initial'][:,i])
            for x,y in ((a['initial']@v,b['initial']@v),(a['initial'][:,others]@v[others],b['initial'][:,others]@v[others])): initial_errors.append(float(abs(x/y-1).max()))
        assert max(initial_errors) <= 2e-13
        stored = next(a for a in summary['amount_evolution'] if (a['history'],a['cell']) == (pair['seed'],pair['cell']))
        for C,h in hs.items():
            amounts[C] = np.asarray([np.asarray(row['chemistry'])@np.asarray(row['volumes']) for row in h])
            surrounding = np.asarray([np.asarray(row['chemistry'])[:,others]@np.asarray(row['volumes'])[others] for row in h])
            norm = np.asarray([float(np.linalg.norm(row['delta'])) for row in h])
            np.testing.assert_array_equal(amounts[C],stored['whole_species_amounts'][C]); np.testing.assert_array_equal(surrounding,stored['surrounding_species_amounts'][C])
            np.testing.assert_array_equal(norm,stored['operator_frobenius_norm'][C]); np.testing.assert_array_equal([row['elapsed'] for row in h],stored['elapsed'])
            gs=[]
            for row in (h[0],h[-1]):
                g = np.asarray(row['volumes'])[:,None]*np.asarray(row['delta']); np.fill_diagonal(g,0.); gs.append(float(g.sum()/2))
            conductance[C] = gs[1]/gs[0]-1; axis[C] = h[-1]['axis_ratio']/h[0]['axis_ratio']-1
            archived[f"history_{pair['seed']}_cell_{pair['cell']}_{C}_amounts"] = amounts[C]
        allocations.append(dict(history=pair['seed'],cell=pair['cell'],initial_amount_relative_error=max(initial_errors),final_global_amount_R_over_S=(amounts['R'][-1]/amounts['S'][-1]).tolist(),maximum_later_global_amount_relative_separation=float(abs(amounts['R']/amounts['S']-1).max()),total_conductance_fractional_change=conductance,axis_ratio_fractional_change=axis))
        for label,horizon in (('pilot',6.),('full',60.)):
            full = read(ROOT/'contrasts'/f"history-{pair['seed']}_cell-{pair['cell']}_{label}.json")
            assert full['passed'] and full['protocol_sha256'] == ph
            for file,expected in full['evidence_sha256'].items(): pin(file,expected)
            waves,metrics = {},{}
            for level in ('coarse','fine'):
                signals=[]; metrics[level]=[]
                for C,context in pair['contexts'].items():
                    job = job_for(context,level,'fractional'); control = job_for(context,level)
                    h,b = [read(ROOT/j['key']/('pilot-history.json' if label == 'pilot' else 'history.json')) for j in (job,control)]
                    lc,w,m = response(h,b,pair['cell'],.9); signals.append(w); metrics[level].append(m)
                    if label == 'full':
                        recoveries.append(dict(history=pair['seed'],cell=pair['cell'],allocation=C,level=level,**m))
                        archived[f"history_{pair['seed']}_cell_{pair['cell']}_{C}_{level}_response"] = w
                waves[level] = signals[1]-signals[0]
            times = np.asarray([r['elapsed'] for r in h]); effects = {level:float(np.sqrt(np.trapezoid(np.mean(w*w,axis=1),times)/horizon)) for level,w in waves.items()}
            error = float(abs(waves['coarse']-waves['fine']).max()); decision = effects['fine'] > .01
            assert decision == (effects['coarse']>.01) and error <= .01 and abs(effects['fine']-effects['coarse']) <= .01
            measured = dict(seed=pair['seed'],cell=pair['cell'],dose='fractional',effects=effects,direct_signed_waveform_error=error,same_decision=True,fine_above_screen=bool(decision),metrics=metrics,passed=True)
            close(measured,full['rows'][0])
            z = arrays(ROOT/'contrasts'/f"history-{pair['seed']}_cell-{pair['cell']}_fractional_{label}.npz",('times','coarse','fine'))
            np.testing.assert_array_equal(z['times'],times)
            for level in ('coarse','fine'): np.testing.assert_array_equal(z[level],waves[level])
            if label == 'full':
                close(measured,next(r for r in summary['response_rows'] if (r['seed'],r['cell']) == (pair['seed'],pair['cell']))); comparisons.append(measured)
                fr = next(r for r in frozen_review['scientific_comparisons'] if (r['seed'],r['cell'],r['matching_reference'],r['budget_reference'],r['physics']) == (pair['seed'],pair['cell'],'S','S','S'))
                frozen_E = fr['full_effects']['fine']
                name = f"{pair['seed']}_{pair['cell']}_S_budget-S_fractional_chemical_sham_physics"
                fw = arrays(frozen_arrays,(name,))[name]
                assert fw.shape == waves['fine'].shape
                delta_wave = waves['fine']-fw
                archived[f"history_{pair['seed']}_cell_{pair['cell']}_moving_minus_frozen"] = delta_wave
                frozen_comparisons.append(dict(history=pair['seed'],cell=pair['cell'],moving_fine_E=effects['fine'],frozen_fine_E=frozen_E,relative_E_change=effects['fine']/frozen_E-1,same_classification=(effects['fine']>.01)==(frozen_E>.01),moving_minus_frozen_signed_RMS=float(np.sqrt(np.trapezoid(np.mean(delta_wave*delta_wave,axis=1),times)/60.)),moving_minus_frozen_signed_max=float(abs(delta_wave).max())))
    outcomes = [dict(history=seed,recipients=2,above_screen=sum(r['fine_above_screen'] for r in comparisons if r['seed'] == seed)) for seed in (7,8,9)]
    close(outcomes,summary['history_outcomes']); assert all(r['same_classification'] for r in frozen_comparisons)
    fine_recoveries = [r for r in recoveries if r['level'] == 'fine']
    assert len(fine_recoveries) == 12 and all(r['target_recovery_time'] is not None and r['network_recovery_time'] is not None for r in fine_recoveries)
    archived['times'] = times; np.savez_compressed(HERE/'independent_review.npz',**archived)
    with (HERE/'response-comparisons.csv').open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=['history','cell','moving_fine_E','frozen_fine_E','relative_E_change','same_classification','moving_minus_frozen_signed_RMS','moving_minus_frozen_signed_max']); writer.writeheader(); writer.writerows(frozen_comparisons)
    for file,expected in PINS.items(): assert digest(file) == expected,file
    record = dict(passed=True,accepted_at_tested_scope=True,utc=datetime.now(timezone.utc).isoformat(),protocol_sha256=ph,method=p['method'],total_new_paths=48,logical_starts=24,distinct_initial_value_problems=24,independent_histories=3,new_histories=0,
        actual_start_gates_reconstructed=48,pilot_pairs_reconstructed=12,full_pairs_reconstructed=12,pilot_signed_contrasts_reconstructed=6,full_signed_contrasts_reconstructed=6,
        source_pins_unchanged=len(p['source_sha256']),input_pins_unchanged=len(p['input_sha256']),historical_guards_checked=len(guards['preserved_ledger_manuscript_sha256']),historical_guards_unchanged=len(guards['preserved_ledger_manuscript_sha256'])-len(exceptions),authorized_editorial_exceptions=exceptions,ledger_and_current_manuscript_unchanged_by_review=True,
        prefix_rows=prefix_rows,mechanical_rows=mechanical_rows,response_refinement_rows=response_rows,preparations=preparations,response_comparisons=comparisons,history_outcomes=outcomes,allocations=allocations,frozen_comparisons=frozen_comparisons,
        recovery=dict(fine_labels=12,all_recover=True,target_range=[min(r['target_recovery_time'] for r in fine_recoveries),max(r['target_recovery_time'] for r in fine_recoveries)],network_range=[min(r['network_recovery_time'] for r in fine_recoveries),max(r['network_recovery_time'] for r in fine_recoveries)],maximum_final_network_log_RMS=max(r['final_log_rms_from_control'] for r in fine_recoveries),maximum_late_network_log_RMS=max(r['late_max_log_rms_from_control'] for r in fine_recoveries),fine_metrics=fine_recoveries),
        quality=dict(maximum_recorded_volume_error=max(r['audit']['max_volume_error'] for r in results),minimum_recorded_radius_grid_intervals=min(r['audit']['min_radius'] for r in results),maximum_recorded_boundary_occupancy=max(r['audit']['boundary_max'] for r in results),maximum_recorded_dilution_amount_error=max(r['audit']['dilution_error_max'] for r in results),clipping=0,whole_window_minimum_log_activator_SD=min(r['whole_window_min_log_sd'] for r in results)),
        coordinator_hours=timing['seconds']/3600,qualification_minutes=q['seconds']/60,evidence_sha256=PINS,script_sha256=digest(__file__),archived_arrays_sha256=digest(HERE/'independent_review.npz'),csv_sha256=digest(HERE/'response-comparisons.csv'),scope=p['scope'],
        limitations='Saved-array reconstruction plus committed eight-step independent chemistry and four-step exact restart records; unsaved branches are not independently rerun. Native references cover 0.6 units, not the entire trajectory. Per-step audit extrema are checked against saved observations. Two timesteps qualify the tested moving paths without assigning order two. No fresh development, spatial convergence, continuously clamped budgets, settled response, nearest-neighbor specificity, autonomous/inherited identity or reciprocal-reference extension. RMS changes are descriptive rather than additive causal fractions.')
    (HERE/'independent_review.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:record[k] for k in ('passed','total_new_paths','actual_start_gates_reconstructed','history_outcomes','recovery','quality','coordinator_hours')},indent=2),flush=True)


if __name__ == '__main__': main()
