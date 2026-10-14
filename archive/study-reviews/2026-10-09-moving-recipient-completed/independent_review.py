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

ROOT = Path('output/studies/moving-recipient-amount-response')
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
    HERE.mkdir(parents=True,exist_ok=True)
    p = read(ROOT/'protocol.json');ph = digest(ROOT/'protocol.json');summary = read(ROOT/'summary.json')
    status = read(ROOT/'status.json');timing = read(ROOT/'timing.json')
    assert status['state'] == 'completed' and status['passed'] and summary['passed']
    assert summary['protocol_sha256'] == ph and summary['total_paths'] == 96
    assert p['histories'] == [7,8,9] and p['new_histories'] == 0
    assert p['primary_reference'] == p['initial_physics'] == 'S'
    assert len(p['jobs']) == 96 and len(p['contexts']) == 12 and len(p['pairs']) == 6
    for collection in ('source_sha256','input_sha256'):
        for file,expected in p[collection].items():pin(file,expected)
    guard_manifest = read('archive/study-document-snapshots/2026-10-08-before-moving-recipient-amount-response/manifest.json')
    guard_differences = []
    for file,expected in guard_manifest['unchanged_guards'].items():
        actual = digest(file)
        if actual != expected:
            # Keep the failed historical export guard explicit. This file is
            # neither a numerical source nor an experiment input. Its original
            # bytes remain independently available in manuscript/main.pdf.
            assert Path(file).resolve() == Path('output/pdf/embryogenesis_o1.pdf').resolve()
            assert str(Path(file).resolve()) not in p['source_sha256']
            assert str(Path(file).resolve()) not in p['input_sha256']
            pin('manuscript/main.pdf',expected)
            guard_differences.append(dict(path=file,expected=expected,current=actual,
                original_bytes_preserved_in='manuscript/main.pdf',role='Non-input manuscript PDF export'))
            pin(file,actual)
        else:
            pin(file,expected)
    for file,expected in summary['evidence_sha256'].items():pin(file,expected)

    histories, results, sources, prefix_rows = {}, [], {}, []
    for number,job in enumerate(p['jobs'],1):
        folder = ROOT/job['key'];r = read(folder/'result.json');h = read(folder/'history.json')
        assert r['job'] == job and r['protocol_sha256'] == ph and r['method'] == p['method'] and r['quality_pass']
        for file,key in ((folder/'history.json','history_sha256'),(folder/'latest_state.npz','checkpoint_sha256'),(job['source'],'source_sha256')):
            pin(file,r[key])
        source = arrays(job['source'],('initial','ids','volumes'));sources[job['key']] = source
        u,v,d,sd = check_history(h,job,60.,source);histories[job['key']] = h
        cp = arrays(folder/'latest_state.npz',('metadata','activator','inhibitor','chemical_correction','chemical_method','phase_carry','rounding','phi','polarity','ids','parents','target','fate','long_experiment','precision_arm','attribute_mode'))
        meta = json.loads(str(cp['metadata']));commit = json.loads(str(cp['long_experiment']))
        assert meta['config'] == p['accepted_configs'][f"{job['seed']}_{job['level']}"]
        assert abs(meta['time']-job['start']-60.) < 1e-9 and meta['step_number'] == round(meta['time']/job['dt'])
        assert not meta['divisions'] and commit['protocol_hash'] == ph and commit['job'] == job
        assert commit['history'] == h and commit['audit'] == r['audit']
        np.testing.assert_array_equal(u[-1],[cp['activator'],cp['inhibitor']])
        np.testing.assert_array_equal(h[-1]['polarity'],cp['polarity'])
        np.testing.assert_array_equal(source['ids'],cp['ids'])
        e = cp['chemical_correction'];last = np.asarray([cp['activator'],cp['inhibitor']])
        assert e.dtype == np.float64 and e.shape == last.shape and np.isfinite(e).all() and np.all(abs(e)<=2*np.spacing(last))
        assert str(cp['chemical_method']) == p['method'] and str(cp['precision_arm']) == 'phase_carry' and str(cp['attribute_mode']) == 'direct'
        assert cp['phi'].dtype == np.float32 and cp['phase_carry'].dtype == np.float64 and cp['phi'].shape == cp['phase_carry'].shape == (16,72,72,72)
        assert np.isfinite(cp['phase_carry']).all() and np.isfinite(cp['phi']).all() and cp['phi'].min() >= 0 and cp['phi'].max() <= 1
        assert np.all(cp['fate']==0)
        audit = r['audit'];assert audit['max_volume_error'] < .05 and audit['min_radius'] >= 4 and audit['max_clipping'] == 0
        assert audit['boundary_max'] < .01 and audit['dilution_error_max'] <= p['amount_error_max']
        assert abs(v/cp['target']-1).max() <= audit['max_volume_error']+1e-12
        assert max(x['boundary_occupancy'] for x in h) <= audit['boundary_max']+1e-16
        close(contrast(h,60.),r);results.append(r);del cp

        pref = folder/'prefix';q = read(pref/'result.json')
        assert q['passed'] and q['job'] == job and q['protocol_sha256'] == ph and q['method'] == p['method']
        pin(job['source'],q['source_sha256'])
        for file,expected in q['evidence_sha256'].items():pin(file,expected)
        gh,ch = read(pref/'gpu-history.json'),read(pref/'cpu-history.json')
        for hh in (gh,ch):check_history(hh,job,.6,source)
        measured = errors(gh,ch);close(measured,q['errors'])
        for key,value in measured.items():assert value <= p['prefix_limits'][key]
        a,b = [arrays(pref/f'{arm}-phi.npz',('phi',))['phi'] for arm in ('gpu','cpu')]
        phi_error = float(abs(a-b).max());close(phi_error,q['final_phi_abs_max'])
        assert phi_error <= p['prefix_limits']['final_phi_abs_max'];del a,b
        assert q['independent_chemical_steps'] == 8 and q['restart_exact_steps'] == 4
        assert q['cpu_chemical_log_max'] <= 1e-11 and q['amount_accounting_max'] <= 2e-14
        cp = arrays(pref/'restart.npz',('metadata','activator','inhibitor','chemical_correction','chemical_method'))
        cu = np.asarray([cp['activator'],cp['inhibitor']]);ce = cp['chemical_correction']
        assert ce.dtype == np.float64 and ce.shape == cu.shape and np.isfinite(ce).all() and np.all(abs(ce)<=2*np.spacing(cu))
        assert str(cp['chemical_method']) == p['method']
        assert abs(json.loads(str(cp['metadata']))['time']-job['start']-8*job['dt']) < 1e-9
        for key in ('volumes','polarity','centers','delta'):
            np.testing.assert_array_equal(gh[0][key],h[0][key])
        prefix_rows.append(dict(job=job['key'],stage=job['stage'],errors=measured,final_phi_abs_max=phi_error,
            runtime_eight_step_chemical_error=q['cpu_chemical_log_max'],runtime_dilution_amount_error=q['amount_accounting_max'],passed=True))
        if number%12 == 0:print(f'Verified {number}/96 completed paths and saved prefix arrays',flush=True)

    def job_for(context,stage,level,dose=None):
        return next(j for j in p['jobs'] if (j['context'],j['stage'],j['level'],j['dose'])==(context,stage,level,dose))

    preparations = []
    for c in p['contexts']:
        fine = job_for(c['key'],'wait','fine');parent = ROOT/fine['key']/'latest_state.npz'
        cp = arrays(parent,('activator','inhibitor','chemical_correction','ids'))
        initial = np.asarray([cp['activator'],cp['inhibitor']]);i = cp['ids'].tolist().index(c['recipient'])
        volumes = np.asarray(histories[fine['key']][-1]['volumes'])
        other = next(x for x in p['contexts'] if (x['seed'],x['recipient'])==(c['seed'],c['recipient']) and x['chemical_origin']!=c['chemical_origin'])
        opposite = job_for(other['key'],'wait','fine');oh = histories[opposite['key']][-1]
        common = .1*min(initial[0,i]*volumes[i],oh['chemistry'][0][i]*oh['volumes'][i])
        doses = []
        for dose in (None,'fractional','amount'):
            job = job_for(c['key'],'response','fine',dose);file = Path(job['source']);m = read(file.with_suffix('.json'))
            pin(file,m['paths_sha256']);pin(parent,m['parent_checkpoint_sha256'])
            assert Path(m['parent_checkpoint']) == parent.resolve() and m['protocol_sha256'] == ph and m['method'] == p['method'] and m['context'] == c['key'] and m['dose'] == dose
            pin(ROOT/opposite['key']/'latest_state.npz',m['opposite_checkpoint_sha256'])
            z = arrays(file,('initial','correction','ids','volumes'));expected = initial.copy();correction = cp['chemical_correction'].copy()
            factor = 1. if dose is None else .9 if dose=='fractional' else 1-common/(initial[0,i]*volumes[i])
            if dose:expected[0,i] *= factor;correction[0,i] = 0.
            np.testing.assert_array_equal(z['initial'],expected);np.testing.assert_array_equal(z['correction'],correction)
            np.testing.assert_array_equal(z['ids'],cp['ids']);np.testing.assert_array_equal(z['volumes'],volumes)
            close(factor,m['record']['factor']);removal = float((initial[0,i]-expected[0,i])*volumes[i])
            close(removal,m['record']['removed_activator_amount'])
            if dose=='amount':assert abs(removal/common-1) <= 2e-13
            for level in ('coarse','fine'):
                hh = histories[job_for(c['key'],'response',level,dose)['key']]
                for key in ('volumes','delta','polarity','centers','axis_ratio'):
                    np.testing.assert_array_equal(hh[0][key],histories[fine['key']][-1][key])
            doses.append(dict(dose=dose,factor=float(factor),removed_amount=removal))
        preparations.append(dict(context=c['key'],actual_new_fine_start_preserved=True,correction_intervention_verified=True,doses=doses))

    mechanical_rows,response_rows = [],[]
    for c in p['contexts']:
        for stage in ('wait','response'):
            for label,horizon in (('pilot',6.),('full',60.)):
                record = read(ROOT/'pairs'/c['key']/f'{stage}-{label}.json')
                assert record['passed'] and record['protocol_sha256'] == ph
                for file,expected in record['evidence_sha256'].items():pin(file,expected)
                hs = {}
                for n,dose in enumerate((None,) if stage=='wait' else (None,'fractional','amount')):
                    paths,fields = [],[]
                    for level in ('coarse','fine'):
                        job = job_for(c['key'],stage,level,dose);folder = ROOT/job['key']
                        h = read(folder/('pilot-history.json' if label=='pilot' else 'history.json'))
                        check_history(h,job,horizon,sources[job['key']]);assert h == histories[job['key']][:len(h)]
                        hs[level,dose] = h;paths.append(h)
                        fields.append(arrays(folder/('pilot-phi.npz' if label=='pilot' else 'latest_state.npz'),('phi',))['phi'])
                    measured = errors(*paths);measured['final_phi_abs_max'] = float(abs(fields[0]-fields[1]).max())
                    close(measured,record['mechanical'][n]['errors'])
                    for key,value in measured.items():assert value <= p['mechanical_limits'][key]
                    close([contrast(h,horizon) for h in paths],record['mechanical'][n]['contrast_summaries'])
                    assert record['mechanical'][n]['same_contrast_outcomes'] and record['mechanical'][n]['loss_times'] == [None,None]
                    mechanical_rows.append(dict(context=c['key'],stage=stage,label=label,dose=dose,errors=measured,passed=True));del fields
                if stage=='response':
                    for n,dose in enumerate(('fractional','amount')):
                        job = job_for(c['key'],stage,'fine',dose);factor = read(Path(job['source']).with_suffix('.json'))['record']['factor']
                        data = [response(hs[level,dose],hs[level,None],c['recipient'],factor) for level in ('coarse','fine')]
                        normalized = float(abs(data[0][0]-data[1][0]).max())
                        auc = abs(data[1][2]['target_activator_log_auc_per_log_pulse']/data[0][2]['target_activator_log_auc_per_log_pulse']-1)
                        recover = {}
                        for key in ('target_recovery_time','network_recovery_time'):
                            a,b = [x[2][key] for x in data]
                            recover[key] = 0. if a is None and b is None else None if a is None or b is None else abs(a-b)
                        measured = dict(normalized_response_error=normalized,relative_auc_error=auc,recovery_time_errors=recover)
                        close(measured,record['response'][n])
                        assert normalized <= .01 and auc <= .02 and all(x is not None and x<=.3+1e-9 for x in recover.values())
                        response_rows.append(dict(context=c['key'],label=label,dose=dose,**measured,passed=True))
        print('Reconstructed pilot/full gates:',c['key'],flush=True)

    comparisons, recoveries, allocations, archived = [],[],[],{}
    for pair in p['pairs']:
        hs = {C:histories[job_for(context,'wait','fine')['key']] for C,context in pair['contexts'].items()}
        i = hs['S'][0]['ids'].index(pair['cell']);x,y = [np.asarray(h[0]['chemistry']) for h in hs.values()]
        v = np.asarray(hs['S'][0]['volumes']);others = np.arange(16)!=i
        np.testing.assert_array_equal(x[:,i],y[:,i]);np.testing.assert_array_equal(hs['S'][0]['delta'],hs['R'][0]['delta'])
        np.testing.assert_array_equal(v,hs['R'][0]['volumes'])
        initial_error = max(float(abs((x@v)/(y@v)-1).max()),float(abs((x[:,others]@v[others])/(y[:,others]@v[others])-1).max()))
        assert initial_error <= 2e-13
        finals = {C:np.asarray(h[-1]['chemistry'])[:,i] for C,h in hs.items()}
        ratio = finals['R']/finals['S'];conductance = {}
        for C,h in hs.items():
            sums=[]
            for row in (h[0],h[-1]):
                g = np.asarray(row['volumes'])[:,None]*np.asarray(row['delta']);np.fill_diagonal(g,0.);sums.append(float(g.sum()/2))
            conductance[C] = sums[1]/sums[0]-1
            amounts = np.asarray([np.asarray(row['chemistry'])@np.asarray(row['volumes']) for row in h])
            stored = next(a for a in summary['wait_amount_evolution'] if (a['history'],a['cell'])==(pair['seed'],pair['cell']))
            np.testing.assert_array_equal(amounts,stored['whole_species_amounts'][C])
            archived[f"history_{pair['seed']}_cell_{pair['cell']}_{C}_amounts"] = amounts
        allocations.append(dict(history=pair['seed'],cell=pair['cell'],initial_amount_relative_error=initial_error,
            final_recipient_pairs={C:value.tolist() for C,value in finals.items()},final_R_over_S=ratio.tolist(),activator_fold_separation=float(max(ratio[0],1/ratio[0])),
            total_conductance_fractional_change=conductance,final_axis_ratio_R_over_S=hs['R'][-1]['axis_ratio']/hs['S'][-1]['axis_ratio']))
        full = read(ROOT/'contrasts'/f"history-{pair['seed']}_cell-{pair['cell']}_full.json")
        assert full['passed'] and full['protocol_sha256'] == ph
        for file,expected in full['evidence_sha256'].items():pin(file,expected)
        for n,dose in enumerate(('fractional','amount')):
            waves,metrics = {},{}
            for level in ('coarse','fine'):
                signals=[];metrics[level]=[]
                for C,context in pair['contexts'].items():
                    job = job_for(context,'response',level,dose);control = job_for(context,'response',level)
                    factor = read(Path(job['source']).with_suffix('.json'))['record']['factor']
                    h,b = histories[job['key']],histories[control['key']];lc,w,m = response(h,b,pair['cell'],factor)
                    signals.append(w);metrics[level].append(m)
                    recoveries.append(dict(history=pair['seed'],cell=pair['cell'],allocation=C,dose=dose,level=level,**m))
                    archived[f"history_{pair['seed']}_cell_{pair['cell']}_{C}_{dose}_{level}_response"] = w
                waves[level] = signals[1]-signals[0]
            times = np.asarray([r['elapsed'] for r in h]);effects = {level:float(np.sqrt(np.trapezoid(np.mean(w*w,axis=1),times)/60)) for level,w in waves.items()}
            error = float(abs(waves['coarse']-waves['fine']).max());decision = effects['fine']>.01
            assert decision == (effects['coarse']>.01) and error <= .01 and abs(effects['fine']-effects['coarse']) <= .01
            measured = dict(seed=pair['seed'],cell=pair['cell'],dose=dose,effects=effects,direct_signed_waveform_error=error,same_decision=True,fine_above_screen=bool(decision),metrics=metrics,passed=True)
            close(measured,full['rows'][n]);close(measured,next(r for r in summary['response_rows'] if (r['seed'],r['cell'],r['dose'])==(pair['seed'],pair['cell'],dose)))
            z = arrays(ROOT/'contrasts'/f"history-{pair['seed']}_cell-{pair['cell']}_{dose}_full.npz",('times','coarse','fine'))
            np.testing.assert_array_equal(z['times'],times)
            for level in ('coarse','fine'):np.testing.assert_array_equal(z[level],waves[level])
            comparisons.append(measured)
    outcomes = [dict(history=seed,recipients=2,above_screen_at_both_doses=sum(all(r['fine_above_screen'] for r in comparisons if (r['seed'],r['cell'])==(seed,pair['cell'])) for pair in p['pairs'] if pair['seed']==seed)) for seed in (7,8,9)]
    close(outcomes,summary['history_outcomes'])
    fine_recoveries = [r for r in recoveries if r['level']=='fine']
    assert len(fine_recoveries)==24 and all(r['target_recovery_time'] is not None and r['network_recovery_time'] is not None for r in fine_recoveries)
    archived['times'] = times
    np.savez_compressed(HERE/'independent_review.npz',**archived)
    with (HERE/'response-comparisons.csv').open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=['history','cell','dose','coarse_E','fine_E','above_screen','signed_waveform_error'])
        writer.writeheader()
        for r in comparisons:writer.writerow(dict(history=r['seed'],cell=r['cell'],dose=r['dose'],coarse_E=r['effects']['coarse'],fine_E=r['effects']['fine'],above_screen=r['fine_above_screen'],signed_waveform_error=r['direct_signed_waveform_error']))
    # Verify files again without the digest cache before declaring acceptance.
    for file,expected in PINS.items():assert digest(file)==expected,file
    record = dict(passed=True,accepted_at_tested_scope=True,utc=datetime.now(timezone.utc).isoformat(),protocol_sha256=ph,
        method=p['method'],total_new_paths=96,wait_paths=24,response_paths=72,independent_histories=3,new_histories=0,
        actual_start_gates_reconstructed=len(prefix_rows),pilot_pairs_reconstructed=24,full_pairs_reconstructed=24,
        source_pins_unchanged=len(p['source_sha256']),input_pins_unchanged=len(p['input_sha256']),
        historical_guards_checked=len(guard_manifest['unchanged_guards']),
        historical_guards_unchanged=len(guard_manifest['unchanged_guards'])-len(guard_differences),
        historical_guards_all_unchanged=not guard_differences,historical_guard_differences=guard_differences,
        prefix_rows=prefix_rows,mechanical_rows=mechanical_rows,response_refinement_rows=response_rows,
        preparations=preparations,response_comparisons=comparisons,history_outcomes=outcomes,allocations=allocations,
        recovery=dict(fine_labels=24,all_recover=True,target_range=[min(r['target_recovery_time'] for r in fine_recoveries),max(r['target_recovery_time'] for r in fine_recoveries)],
            network_range=[min(r['network_recovery_time'] for r in fine_recoveries),max(r['network_recovery_time'] for r in fine_recoveries)],
            maximum_late_network_log_RMS=max(r['late_max_log_rms_from_control'] for r in fine_recoveries)),
        quality=dict(maximum_recorded_volume_error=max(r['audit']['max_volume_error'] for r in results),minimum_recorded_radius_grid_intervals=min(r['audit']['min_radius'] for r in results),
            maximum_recorded_boundary_occupancy=max(r['audit']['boundary_max'] for r in results),maximum_recorded_dilution_amount_error=max(r['audit']['dilution_error_max'] for r in results),clipping=0,
            whole_window_minimum_log_activator_SD=min(r['whole_window_min_log_sd'] for r in results)),
        coordinator_hours=timing['seconds']/3600,evidence_sha256=PINS,script_sha256=digest(__file__),
        archived_arrays_sha256=digest(HERE/'independent_review.npz'),csv_sha256=digest(HERE/'response-comparisons.csv'),
        scope=p['scope'],limitations='Saved-array reconstruction plus committed eight-step independent chemistry and exact restart records; unsaved branches are not independently rerun. Native references cover 0.6 units, not the entire trajectory. Two timesteps qualify the tested moving paths without assigning order two. No new development, spatial convergence, ongoing budget clamp, pulse-time local-state matching, autonomous/inherited identity or entire frozen factorial.')
    (HERE/'independent_review.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:record[k] for k in ('passed','total_new_paths','actual_start_gates_reconstructed','history_outcomes','recovery','quality','coordinator_hours')},indent=2),flush=True)


if __name__=='__main__':main()
