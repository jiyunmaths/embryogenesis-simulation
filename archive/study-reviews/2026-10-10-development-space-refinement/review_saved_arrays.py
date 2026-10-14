"""Independent spatial/time review using saved arrays; no embryo imports.

Reuse the earlier independent CDF, eigenproblem and checkpoint checks. Rebuild
every spatial primary/endpoint field and every native finer-grid time field.
Keep the production summary, older failures and reused bytes unchanged.
"""
from pathlib import Path
from datetime import datetime, timezone
import importlib.util
import json
import time
import numpy as np
from scipy.ndimage import map_coordinates

ROOT = Path('output/studies/development-space-refinement')
HERE = Path(__file__).parent
HELPER = Path('archive/study-reviews/2026-10-10-development-time-diagnosis/review_saved_arrays.py')
spec = importlib.util.spec_from_file_location('independent_time_review', HELPER)
t = importlib.util.module_from_spec(spec); spec.loader.exec_module(t)
r = t.r
read, load, pin, close = r.read, r.load, r.pin, r.close


def density(z,probe,order=3):
    """Cellwise spline on an affine common grid, then float64 occupancy.

    No phase clipping or alignment. Keep interpolated phase in native float32
    before evaluating occupancy, matching the prescribed diagnostic precision.
    """
    phi = z['phi']; grid = phi.shape[-1]
    if probe is None:
        shape = phi.shape[1:]; coordinates = None
    else:
        shape = (probe,)*3
        axis = (np.arange(probe)+.5)*grid/probe-.5
        coordinates = np.array(np.meshgrid(axis,axis,axis,indexing='ij'))
    out = np.zeros((6,)+shape,dtype=float)
    u = np.stack((z['activator'],z['inhibitor'])); lo,hi = np.inf,-np.inf
    for i,f in enumerate(phi):
        mapped = (map_coordinates(f,coordinates,order=order,mode='nearest',prefilter=True)
            if probe is not None and probe != grid else f)
        lo,hi = min(lo,float(mapped.min())),max(hi,float(mapped.max()))
        p = mapped.astype(float); h = p*p*(3-2*p)
        out[0] += h
        for j in range(2): out[1+j] += u[j,i]*h
        for j in range(3): out[3+j] += z['polarity'][i,j]*h
    return out,[lo,hi]


def dense_compare(a,b):
    rows = []; unmatched = []; topology = []; conductance = []
    assert sorted(a) == sorted(b) == [round(i*.015,9) for i in range(1201)]
    for time,x in a.items():
        y = b[time]; sx,sy = r.eig(x)[0],r.eig(y)[0]
        if len(sx) != len(sy): unmatched.append(time); continue
        er = float(np.linalg.norm(sx-sy)/(np.linalg.norm(sy) or 1.))
        rows.append(dict(time=time,relative_spectrum=er,cells=len(sx),edges=[x['edges'],y['edges']]))
        if x['ids'] == y['ids'] and x['parents'] == y['parents']:
            ga,gb = np.array(x['conductance']),np.array(y['conductance'])
            conductance.append(dict(time=time,relative_conductance=float(np.linalg.norm(ga-gb)/(np.linalg.norm(gb) or 1.))))
            if not np.array_equal(ga>0,gb>0): topology.append(time)
    return dict(peak=max(rows,key=lambda x:x['relative_spectrum']),unmatched_times=unmatched,
        topology_mismatch_times=topology,conductance_peak=max(conductance,key=lambda x:x['relative_conductance']))


def compare_histories(original,histories,reports,networks):
    a,b = original['cases']; errors,peaks,unmatched,raw,raw_unmatched = r.compare(histories[a],histories[b])
    for k,v in errors.items(): close(v,original['errors'][k],('history',a,b,k))
    assert unmatched == original['unmatched_times']
    maturity = abs(max(e['time'] for e in reports[a]['events'])-max(e['time'] for e in reports[b]['events']))
    close(maturity,original['maturity_time_error'],'maturity')
    return dict(cases=[a,b],errors=errors,peak_times=peaks,unmatched_times=unmatched,
        maturity_time_error=maturity,field_checks=[],dense=dense_compare(networks[a],networks[b]),
        production_passed=original['passed'])


def decision(row,original,limits,spread_limit=None):
    fm = max(v for f in row['field_checks'] for v in f['errors'].values())
    close(fm,original['field_max'],'field max')
    failures = {k:v for k,v in row['errors'].items() if v > limits[k]}
    if fm > limits['field_max']: failures['field_max'] = fm
    if row['maturity_time_error'] > limits['maturity_time_difference_max']:
        failures['maturity_time_error'] = row['maturity_time_error']
    spread = {k:max(f['errors'][k] for f in row['field_checks'] if f['time']==18)-
        min(f['errors'][k] for f in row['field_checks'] if f['time']==18) for k in r.FIELDS}
    for k,v in spread.items(): close(v,original['diagnostic_spread'][k],('spread',k))
    if spread_limit is not None and max(spread.values()) > spread_limit:
        failures['diagnostic_spread'] = max(spread.values())
    row.update(field_max=fm,diagnostic_spread=spread,failed_gates=failures,independently_passed=not failures)
    assert row['independently_passed'] == row['production_passed']


def main():
    started = time.perf_counter(); pin(__file__); pin(HELPER); pin(t.HELPER)
    p = read(ROOT/'protocol.json'); ph = pin(ROOT/'protocol.json')
    s = read(ROOT/'summary.json'); status = read(ROOT/'status.json'); timing = read(ROOT/'timing.json')
    assert status['state'] == 'completed_provisional' and not s['qualification_accepted']
    assert s['protocol_sha256'] == ph and s['new_scientific_histories'] == 0
    for path,expected in {**p['source_sha256'],**p['input_sha256'],**s['evidence_sha256']}.items(): pin(path,expected)
    runtime = read(ROOT/'runtime.json'); pin(runtime['binary'],runtime['binary_sha256'])
    prior_runtime = read('output/studies/development-time-diagnosis/runtime.json')
    assert runtime['binary_sha256'] == prior_runtime['binary_sha256']
    calibration = read(ROOT/'calibration.json'); cm = read(ROOT/'calibration-manifest.json')
    assert calibration['passed'] and len(calibration['rows']) == 48
    assert calibration['protocol_sha256'] == cm['protocol_sha256'] == ph
    assert cm['calibration_sha256'] == pin(ROOT/'calibration.json')
    assert max(c['relative_l2'] for c in calibration['rows']) <= p['calibration_limit']
    gate = read(p['instrumentation_gate']); assert gate['passed']
    assert gate['protocol_sha256'] == p['reuse_provenance']['protocol_sha256']
    for file,expected in gate['evidence_sha256'].items(): pin(file,expected)
    histories = {}; reports = {}; networks = {}; paths = []; lineages = {}
    for job in p['jobs']:
        key = job['key']; folder = ROOT/key; cfg = job['config']
        report = read(folder/'gpu-result.json'); reports[key] = report
        assert report['passed'] and report['failure'] is None
        if job['reused']:
            provenance = p['reuse_provenance']
            assert folder.is_symlink() and folder.resolve() == Path(provenance['folder'])
            assert pin(folder/'gpu-result.json') == provenance['report_sha256']
            assert report['job'] == provenance['job'] and report['protocol_sha256'] == provenance['protocol_sha256']
            assert report['job']['config'] == cfg
        else:
            assert not folder.is_symlink() and report['job'] == job and report['protocol_sha256'] == ph
        for file,expected in report['evidence_sha256'].items(): pin(file,expected)
        history = read(folder/'gpu-history.json'); histories[key] = history
        events = read(folder/'gpu-events.json'); assert events == report['events']
        lineage = read(folder/'gpu-lineage.json'); lineages[key] = lineage
        assert len(history) == 121 and len(events) == 15 and len(lineage) == 31
        assert report['final_time'] == 18 and report['final_cells'] == 16
        audit = report['audit']; limits = p['physical_limits']
        assert audit['max_volume_error'] < limits['volume_error_max'] and audit['min_radius'] >= limits['radius_min']
        assert audit['max_clipping'] == 0 and audit['dilution_amount_error'] <= limits['dilution_amount_error_max']
        assert report['exact_restart_steps'] == {'2':8,'4':8,'8':8}
        assert report['independent_chemical_steps'] == 8 and report['independent_chemical_log_max'] <= 1e-11
        assert 18-max(e['time'] for e in events) >= p['minimum_mature_window']
        assert len({e['parent'] for e in events}) == 15
        for e in events:
            assert e['neck'] <= .15 and e['overlap'] <= .002 and .1 < e['fraction'] < .9
            assert e['daughter_components'] == [1,1]
            for k,l in [('amount_relative_jump','event_amount_jump_max'),('target_relative_jump','event_target_jump_max'),
                ('volume_relative_jump','event_volume_jump_max'),('effective_occupancy_absolute_jump','event_occupancy_jump_max')]:
                assert e[k] <= limits[l]
            for child in e['children']: assert lineage[child]['parent'] == e['parent']
            np.testing.assert_allclose(np.array(e['chemistry_before'])@np.array(e['volumes_before']),
                np.array(e['chemistry_after'])@np.array(e['volumes_after']),rtol=2e-14,atol=1e-15)
        closure = 0.; boundary = 0.
        for i,row in enumerate(history):
            assert abs(row['time']-i*.15) < 1e-9
            n = len(row['ids']); u = np.array(row['chemistry']); v = np.array(row['volumes'])
            assert u.shape == (2,n) and np.isfinite(u).all() and u.min() > 0 and np.isfinite(v).all() and v.min() > 0
            assert all(c == 1 for c in row['nondividing_components'])
            assert len(row['nondividing_components']) == n-len(row['dividing_ids'])
            assert row['boundary_occupancy'] < limits['boundary_max']
            _,err = r.eig(row); closure = max(closure,err); boundary = max(boundary,row['boundary_occupancy'])
        trace = t.jsonl(folder/'network.jsonl'); assert len(trace) == report['network_trace_rows']
        previous = -1.; regular = {}
        for row in trace:
            assert row['time'] > previous; previous = row['time']
            close(row['time'],row['step']*cfg['dt'],'dense clock')
            _,err = t.graph(row,cfg); closure = max(closure,err)
            if abs(row['time']/.015-round(row['time']/.015)) < 1e-7: regular[round(row['time'],9)] = row
        assert len(regular) == 1201; networks[key] = regular
        pressure = t.jsonl(folder/'projection.jsonl'); seen = set(); parents = {e['parent'] for e in events}
        for row in pressure:
            assert row['parent'] in parents and (row['step'],row['parent']) not in seen
            seen.add((row['step'],row['parent']))
            close(row['time_start'],row['step']*cfg['dt'],'projection clock')
            close(row['time_end']-row['time_start'],cfg['dt'],'projection dt')
            assert row['capacity'] > 0 and np.isfinite(list(row.values())).all()
            assert abs(row['after_relative_volume_error']) <= 1e-7+1e-14 and row['correction_magnitude'] >= 0
        extrema = dict(calls=len(pressure),correction_magnitude_max=max(x['correction_magnitude'] for x in pressure),
            before_relative_volume_error_max=max(abs(x['before_relative_volume_error']) for x in pressure),
            after_relative_volume_error_max=max(abs(x['after_relative_volume_error']) for x in pressure))
        for k,v in extrema.items(): close(v,report['projection'][k],('projection extrema',k))
        paths.append(dict(case=key,reused=job['reused'],seconds=report['seconds'],audit=audit,
            boundary_max=boundary,graph_closure_max=closure,projection=extrema,network_rows=len(trace),
            recorded_runtime_event_gates_passed=True,log_activator_sd_at_18=history[-1]['log_activator_sd']))
        print('Reviewed histories, dense graphs and pressure records:',key,flush=True)
    space = [compare_histories(row,histories,reports,networks) for row in s['space_rows']]
    temporal = compare_histories(s['time_row'],histories,reports,networks)
    for actual,original in zip([row['dense'] for row in space]+[temporal['dense']],s['dense_space_rows']+[s['dense_time']]):
        close(actual['peak']['relative_spectrum'],original['peak']['relative_spectrum'],'dense spectrum')
        assert actual['peak']['time'] == original['peak']['time']
        assert actual['unmatched_times'] == original['unmatched_times']
        assert actual['topology_mismatch_times'] == original['topology_mismatch_times']
    samples = []; live = []; configs = {j['key']:j['config'] for j in p['jobs']}
    for stamp in (0,3,6,9,12,15,18):
        primary = {}; native = {}; ranges = {}
        for job in p['jobs']:
            key = job['key']; z = load(ROOT/key/f'gpu-t{stamp:02d}.npz')
            check = r.checkpoint(z,job['config'],histories[key][int(round(stamp/.15))],stamp==18)
            assert check['volume_error'] < p['physical_limits']['volume_error_max'] and check['minimum_radius'] >= p['physical_limits']['radius_min']
            assert json.loads(str(z['event_audits'])) == [e for e in reports[key]['events'] if e['time'] <= stamp+1e-12]
            samples.append(dict(case=key,**check))
            if key != 'off120_half': primary[key],ranges[key] = density(z,112,3)
            if key in p['time_pair']: native[key],_ = density(z,None)
            del z
        for i,row in enumerate(space):
            a,b = row['cases']; values = r.field_error(primary[a],primary[b])
            expected = next(f for f in s['space_rows'][i]['fields'] if f['time']==stamp and f['probe']==112 and f['order']==3)
            for k,v in values.items(): close(v,expected['errors'][k],('spatial primary',stamp,a,b,k))
            np.testing.assert_allclose([ranges[a],ranges[b]],expected['phase_ranges'],rtol=1e-6,atol=1e-9)
            row['field_checks'].append(dict(time=stamp,probe=112,order=3,errors=values,phase_ranges=[ranges[a],ranges[b]]))
        a,b = temporal['cases']; values = r.field_error(native[a],native[b])
        expected = next(f for f in s['time_row']['fields'] if f['time']==stamp)
        for k,v in values.items(): close(v,expected['errors'][k],('native temporal',stamp,k))
        temporal['field_checks'].append(dict(time=stamp,probe=None,order=3,errors=values))
        del native
        if stamp == 18:
            for probe,order in ((112,5),(128,3),(128,5)):
                fields = {}; ranges = {}
                for key in ('off72','off88','off104','off120'):
                    z = load(ROOT/key/'gpu-t18.npz'); fields[key],ranges[key] = density(z,probe,order); del z
                for i,row in enumerate(space):
                    a,b = row['cases']; values = r.field_error(fields[a],fields[b])
                    expected = next(f for f in s['space_rows'][i]['fields'] if f['time']==18 and f['probe']==probe and f['order']==order)
                    for k,v in values.items(): close(v,expected['errors'][k],('endpoint variants',probe,order,a,b,k))
                    row['field_checks'].append(dict(time=18,probe=probe,order=order,errors=values,phase_ranges=[ranges[a],ranges[b]]))
                del fields
                print('Reconstructed all spatial endpoint variants:',probe,order,flush=True)
        del primary
        print('Reconstructed all primary spatial and native temporal fields:',stamp,flush=True)
    for job in p['jobs']:
        for stage in (2,4,8):
            z = load(ROOT/job['key']/f'restart-stage-{stage}.npz'); check = r.checkpoint(z,job['config'])
            assert stage < check['cells'] < 2*stage and check['active'] > 0
            live.append(dict(case=job['key'],stage=stage,**check))
    for row,original in zip(space,s['space_rows']): decision(row,original,p['space_limits'],p['diagnostic_spread_limit'])
    decision(temporal,s['time_row'],p['time_limits'])
    keys = ('concentration_wasserstein_max','log_wasserstein_max','relative_axis_max')
    sequences = {k:[row['errors'][k] for row in space] for k in keys}
    sequences.update({k:[max(f['errors'][k] for f in row['field_checks']) for row in space] for k in r.FIELDS})
    trend = all(all(b <= a+p['trend_floor'] for a,b in zip(values,values[1:])) for values in sequences.values())
    assert trend == s['spatial_error_decay_passed']
    dense_pass = temporal['dense']['peak']['relative_spectrum'] <= p['dense_time_spectrum_limit']
    scoped = all(i['recorded_runtime_event_gates_passed'] for i in paths) and space[-1]['independently_passed'] and temporal['independently_passed'] and dense_pass and trend
    assert scoped == s['scoped_screens_passed'] and dense_pass == s['dense_time_screen_passed']
    preserved = read(p['preserved_manifest'])
    for file,expected in preserved['preserved_sha256'].items(): pin(file,expected)
    result = dict(review_completed=True,scoped_screen_accepted=bool(scoped),full_development_qualified=False,
        qualification_accepted=False,reviewed_at_utc=datetime.now(timezone.utc).isoformat(),protocol_sha256=ph,
        source_count=len(p['source_sha256']),input_count=len(p['input_sha256']),preserved_count=len(preserved['preserved_sha256']),
        paths=paths,space_rows=space,time_row=temporal,dense_time_screen_passed=dense_pass,spatial_error_decay_passed=trend,
        spatial_error_sequences=sequences,regular_checkpoints=samples,live_checkpoint_samples=live,
        independent_full_restart_replay_performed=False,
        recorded_checks_not_independently_replayed=['Every-step physical extrema and eight-step exact restarts; saved arrays inspected, full replays remain separate.',
            'Inherited instrumentation gate, initial chemical prefixes and analytic reconstruction calibration; original records and hashes verified.'],
        production_seconds=timing['seconds'],new_scientific_histories=0,new_numerical_paths=4,reused_numerical_paths=1,
        consumed_sha256=r.WATCHED,script_sha256=pin(__file__),review_seconds=time.perf_counter()-started,
        scope=p['scope'])
    target = HERE/'independent_review.json'
    if target.exists(): raise FileExistsError(target)
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('review_completed','scoped_screen_accepted','full_development_qualified','preserved_count','review_seconds')}),flush=True)


if __name__ == '__main__': main()
