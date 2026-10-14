"""Independent review of completed repeated cleavage; no embryo imports.

Recompute all history metrics and all native CPU/GPU and timestep fields.
For spatial pairs, reconstruct one prescribed-primary failed-field witness
per pair. Hash all evidence, and distinguish recorded every-step extrema
and short exact replays from sampled-array checks and full restart replay.
This review cannot accept the already failed study.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import time
import numpy as np
from scipy.linalg import eigh
from scipy.ndimage import label, map_coordinates

ROOT = Path('output/studies/repeated-development-gpu-qualification-v2')
HERE = Path(__file__).parent
WATCHED = {}
FIELDS = ('occupancy', 'activator_density', 'inhibitor_density', 'polarity_density')


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def pin(path, expected=None):
    key = str(Path(path).resolve())
    if key not in WATCHED:
        WATCHED[key] = sha(path)
    if expected is not None:
        assert WATCHED[key] == expected, key
    return WATCHED[key]


def read(path):
    pin(path)
    return json.loads(Path(path).read_text())


def load(path):
    pin(path)
    with np.load(path, allow_pickle=False) as z:
        return {k: z[k] for k in z.files}


def close(a, b, what):
    assert np.isclose(a, b, rtol=2e-7, atol=5e-12), (what, a, b)


def eig(row):
    v = np.asarray(row['volumes']); d = np.asarray(row['delta'])
    k = -v[:, None] * d
    scale = max(1., np.linalg.norm(d))
    errors = [abs(d.sum(1)).max(), abs(v @ d).max(), abs(k-k.T).max()]
    assert max(errors) / scale <= 1e-12
    off = d.copy(); np.fill_diagonal(off, 0)
    assert off.min() >= 0 and np.diag(d).max() <= 0
    values = eigh((k+k.T)/2, np.diag(v), eigvals_only=True)
    assert values.min() >= -1e-12 * scale
    np.testing.assert_allclose(values, row['spectrum'], rtol=2e-12, atol=2e-12)
    return values, float(max(errors)/scale)


def w1(a, b, v, z):
    """Integral of weighted empirical CDF difference, independent of labels."""
    a, b = np.asarray(a), np.asarray(b)
    ia, ib = np.argsort(a), np.argsort(b)
    points = np.unique(np.r_[a, b])
    ca = np.r_[0., np.cumsum(v[ia]/v.sum())]
    cb = np.r_[0., np.cumsum(z[ib]/z.sum())]
    return float(np.sum(np.diff(points) * abs(
        ca[np.searchsorted(a[ia], points[:-1], side='right')] -
        cb[np.searchsorted(b[ib], points[:-1], side='right')])))


def compare(a, b):
    maxima = {}; peaks = {}; unmatched = []; raw = {}; raw_unmatched = []
    assert len(a) == len(b) == 121
    for x, y in zip(a, b):
        assert abs(x['time']-y['time']) < 1e-9
        u, w = np.array(x['chemistry']), np.array(y['chemistry'])
        v, z = np.array(x['volumes']), np.array(y['volumes'])
        values = dict(concentration_wasserstein_max=max(w1(i,j,v,z) for i,j in zip(u,w)),
            log_wasserstein_max=max(w1(np.log(i),np.log(j),v,z) for i,j in zip(u,w)),
            contrast_abs_max=float(abs(u.std(1)-w.std(1)).max()),
            relative_axis_max=abs(x['axis_ratio']/y['axis_ratio']-1),
            relative_total_volume_max=abs(v.sum()/z.sum()-1),
            whole_amount_relative_max=float(abs((u@v)/(w@z)-1).max()),
            mean_polarity_abs_max=float(abs(v@np.array(x['polarity'])/v.sum()-z@np.array(y['polarity'])/z.sum()).max()),
            cell_count_difference_max=abs(len(v)-len(z)))
        sx, _ = eig(x); sy, _ = eig(y)
        if len(v) == len(z):
            values['relative_spectrum_max'] = float(np.linalg.norm(sx-sy)/(np.linalg.norm(sy) or 1.))
        else:
            unmatched.append(x['time'])
        for key, value in values.items():
            if key not in maxima or value > maxima[key]:
                maxima[key] = float(value); peaks[key] = x['time']
        if x['ids'] == y['ids'] and x['parents'] == y['parents']:
            d, e = np.array(x['delta']), np.array(y['delta'])
            rv = dict(chemical_log_max=float(abs(np.log(u/w)).max()),
                polarity_abs_max=float(abs(np.array(x['polarity'])-np.array(y['polarity'])).max()),
                relative_volume_max=float(abs(v/z-1).max()),
                relative_transport_max=float(np.linalg.norm(d-e)/(np.linalg.norm(e) or 1.)))
            for key, value in rv.items(): raw[key] = max(raw.get(key,0.),value)
        else:
            raw_unmatched.append(x['time'])
    maxima['unmatched_observations_max'] = len(unmatched)
    raw['unmatched_observations_max'] = len(raw_unmatched)
    return maxima, peaks, unmatched, raw, raw_unmatched


def checkpoint(z, config, history=None, final=False):
    m = json.loads(str(z['metadata'])); phi = z['phi']; n = len(z['ids'])
    assert m['config'] == config
    assert str(z['chemical_method']) == 'compensated resident progressive development v1'
    assert int(z['development_schema']) == int(z['attribute_schema']) == 1
    assert str(z['attribute_mode']) == ('direct' if config['feedback'] else 'no_feedback')
    assert phi.dtype == np.float32 and phi.shape == (n,)+(config['grid'],)*3
    assert z['phase_carry'].dtype == np.float64 and z['phase_carry'].shape == phi.shape
    assert np.isfinite(z['phase_carry']).all() and phi.min() >= 0 and phi.max() <= 1
    effective = phi.astype(float)+z['phase_carry']
    assert effective.min() >= -1e-15 and effective.max() <= 1+1e-15
    u = np.stack((z['activator'],z['inhibitor'])); correction = z['chemical_correction']
    assert u.dtype == correction.dtype == np.float64 and u.shape == correction.shape == (2,n)
    assert np.isfinite(u).all() and u.min() > 0 and np.isfinite(correction).all()
    assert np.all(abs(correction) <= 2*np.spacing(u))
    assert np.isfinite(z['polarity']).all() and z['polarity'].shape == (n,3)
    assert z['target'].shape == (n,) and z['target'].min() > 0 and np.isfinite(z['target']).all()
    assert np.all(z['fate'] == 0) and z['due'].shape == (n,)
    for cell, due in zip(z['ids'],z['due']):
        assert np.isposinf(due) if str(int(cell)) in m['divisions'] else np.isfinite(due)
    dx = 2*config['extent']/config['grid']
    h = phi*phi*(np.float32(3)-np.float32(2)*phi)
    volumes = h.sum((1,2,3),dtype=float)*dx**3
    assert abs(z['target'].sum()/m['initial_volume']-1) <= 2e-14
    if history is not None:
        assert abs(m['time']-history['time']) < 1e-9
        np.testing.assert_array_equal(z['ids'],history['ids'])
        np.testing.assert_array_equal(z['parents'],history['parents'])
        np.testing.assert_array_equal(u,history['chemistry'])
        np.testing.assert_array_equal(z['polarity'],history['polarity'])
        np.testing.assert_allclose(volumes,history['volumes'],rtol=2e-14,atol=1e-15)
        observed=[label(f>=.5)[1] for cell,f in zip(z['ids'],phi) if str(int(cell)) not in m['divisions']]
        assert observed == history['nondividing_components']
    if final:
        assert n == 16 and m['time'] == 18. and m['step_number'] == config['steps']
        assert m['divisions'] == {} and m['next_id'] == 31 and len(m['lineage']) == 31
        assert all(label(f>=.5)[1] == 1 for f in phi)
    return dict(time=m['time'],cells=n,active=len(m['divisions']),
        volume_error=float(abs(volumes/z['target']-1).max()),
        minimum_radius=float(((3*volumes/(4*np.pi))**(1/3)/dx).min()))


def density(z, probe=None):
    """Interpolate native phi cellwise, then evaluate float64 occupancy."""
    phi=z['phi']; grid=phi.shape[-1]
    if probe is not None:
        axis=(np.arange(probe)+.5)*grid/probe-.5
        coordinates=np.array(np.meshgrid(axis,axis,axis,indexing='ij'))
    shape=(probe,)*3 if probe is not None else phi.shape[1:]
    result=np.zeros((6,)+shape,dtype=float)
    u=np.stack((z['activator'],z['inhibitor']))
    for i, native in enumerate(phi):
        p=(map_coordinates(native,coordinates,order=3,mode='nearest',prefilter=True)
            if probe is not None else native).astype(float)
        h=p*p*(3-2*p)
        result[0]+=h
        for species in range(2): result[1+species]+=u[species,i]*h
        for axis in range(3): result[3+axis]+=z['polarity'][i,axis]*h
    return result


def field_error(a,b):
    norm=lambda v:float(np.sqrt(np.sum(v*v,dtype=float)))
    return dict(zip(FIELDS,[norm(a[i]-b[i])/max(norm(b[i]),1e-30) for i in range(3)] +
        [norm(a[3:]-b[3:])/max(norm(b[0]),1e-30)]))


def main():
    started=time.perf_counter();pin(__file__)
    protocol=read(ROOT/'protocol.json'); protocol_hash=pin(ROOT/'protocol.json')
    production=read(ROOT/'summary.json');status=read(ROOT/'status.json');timing=read(ROOT/'timing.json')
    assert not production['passed'] and status['state']=='completed_with_unresolved_checks'
    assert production['protocol_sha256']==protocol_hash and production['new_scientific_histories']==0
    for path,expected in {**protocol['source_sha256'],**protocol['input_sha256'],**production['evidence_sha256']}.items(): pin(path,expected)
    runtime=read(ROOT/'runtime.json');pin(runtime['binary'],runtime['binary_sha256'])
    calibration=read(ROOT/'calibration.json');cm=read(ROOT/'calibration-manifest.json')
    assert calibration['passed'] and len(calibration['rows'])==36
    assert cm['calibration_sha256']==pin(ROOT/'calibration.json')
    assert cm['protocol_sha256']==calibration['protocol_sha256']==protocol_hash
    histories={};reports={};paths=[];sampled_checkpoints=[];restart_samples=[]
    print('Reviewing 20 histories and recorded runtime/event gates',flush=True)
    for job in protocol['jobs']:
        for arm in (('gpu','cpu') if job['key'] in protocol['cpu_jobs'] else ('gpu',)):
            folder=ROOT/job['key'];key=(job['key'],arm)
            report=read(folder/(arm+'-result.json')); reports[key]=report
            assert report['job']==job and report['protocol_sha256']==protocol_hash and report['passed']
            for path,expected in report['evidence_sha256'].items():pin(path,expected)
            history=read(folder/(arm+'-history.json'));histories[key]=history
            assert len(history)==121 and report['final_time']==18 and report['final_cells']==16
            assert len(report['events'])==15
            audit=report['audit']; limits=protocol['physical_limits']
            assert audit['max_volume_error']<limits['volume_error_max'] and audit['min_radius']>=limits['radius_min']
            assert audit['max_clipping']==0 and audit['dilution_amount_error']<=limits['dilution_amount_error_max']
            for ev in report['events']:
                assert ev['neck']<=.15 and ev['overlap']<=.002 and .1<ev['fraction']<.9 and ev['daughter_components']==[1,1]
                for k,l in [('amount_relative_jump','event_amount_jump_max'),('target_relative_jump','event_target_jump_max'),
                    ('volume_relative_jump','event_volume_jump_max'),('effective_occupancy_absolute_jump','event_occupancy_jump_max')]:assert ev[k]<=limits[l]
            assert 18-max(e['time'] for e in report['events'])>=3
            lineage=read(folder/(arm+'-lineage.json'));assert len(lineage)==31 and lineage[0]['parent']==-1
            assert len(set(ev['parent'] for ev in report['events']))==15
            for ev in report['events']:
                for child in ev['children']:assert lineage[child]['parent']==ev['parent']
            graph_error=0.;boundary=0.
            for index,row in enumerate(history):
                assert abs(row['time']-index*.15)<1e-9
                n=len(row['ids']);u=np.array(row['chemistry']);v=np.array(row['volumes'])
                assert u.shape==(2,n) and np.isfinite(u).all() and u.min()>0 and np.isfinite(v).all() and v.min()>0
                assert all(c==1 for c in row['nondividing_components']) and len(row['nondividing_components'])==n-len(row['dividing_ids'])
                assert row['boundary_occupancy']<.01
                _,g=eig(row);graph_error=max(graph_error,g);boundary=max(boundary,row['boundary_occupancy'])
            assert np.allclose(np.std(np.log(np.array(history[-1]['chemistry'])[0])),history[-1]['log_activator_sd'],rtol=1e-12,atol=1e-15)
            if arm=='gpu':
                assert report['exact_restart_steps']=={'2':8,'4':8,'8':8} and report['independent_chemical_steps']==8
                assert report['independent_chemical_log_max']<=1e-11
            paths.append(dict(case=job['key'],arm=arm,recorded_runtime_event_gates_pass=True,
                normalized_transport_closure_max=graph_error,boundary_max=boundary,audit=audit,
                log_activator_sd_at_18=history[-1]['log_activator_sd']))
    groups={}
    for group,limit_key in [('backend_rows','backend_limits'),('time_rows','time_limits'),('space_rows','space_limits')]:
        rows=[]
        for original in production[group]:
            a,b=original['cases'];aa,ab=original['arms'];errors,peaks,unmatched,raw,raw_unmatched=compare(histories[a,aa],histories[b,ab])
            for k,v in errors.items():close(v,original['errors'][k],(group,a,b,k))
            assert unmatched==original['unmatched_times']
            maturity=abs(max(e['time'] for e in reports[a,aa]['events'])-max(e['time'] for e in reports[b,ab]['events']))
            close(maturity,original['maturity_time_error'],'maturity')
            failures={k:v for k,v in errors.items() if v>protocol[limit_key][k]}
            if maturity>protocol[limit_key]['maturity_time_difference_max']:failures['maturity_time_error']=maturity
            if group=='backend_rows':
                ta=read(ROOT/a/(aa+'-lineage.json'));tb=read(ROOT/b/(ab+'-lineage.json'))
                assert [(r['id'],r['parent']) for r in ta]==[(r['id'],r['parent']) for r in tb]
                ea={e['parent']:e for e in reports[a,aa]['events']};eb={e['parent']:e for e in reports[b,ab]['events']}
                assert ea.keys()==eb.keys()
                assert all(ea[k]['children']==eb[k]['children'] and ea[k]['draws']==eb[k]['draws'] for k in ea)
                event_time=max(abs(ea[k]['time']-eb[k]['time']) for k in ea)
                assert event_time<=.0075+1e-12 and original['lineage_and_partition_correspondence']
                close(event_time,original['parent_event_time_max'],'parent events')
                for k,v in raw.items():close(v,original['homologous_errors'][k],('raw',a,k))
                assert raw_unmatched==original['homologous_unmatched_times']
            rows.append(dict(cases=[a,b],arms=[aa,ab],errors=errors,peak_times=peaks,
                failed_history_gates=failures,maturity_time_error=maturity,
                field_checks=[],production_field_max=original['field_max'],production_passed=original['passed']))
        groups[group]=rows
    print('Reconstructing all native checkpoints and restart samples',flush=True)
    for branch in ('off','tension_adhesion','polarity','full'):
        for stamp in (0,3,6,9,12,15,18):
            zs={}
            for level,arm in [('coarse','gpu'),('fine','gpu'),('coarse','cpu'),('space88','gpu'),('space104','gpu')]:
                case=branch+'_'+level;z=load(ROOT/case/(arm+f'-t{stamp:02d}.npz'))
                checks=checkpoint(z,reports[case,arm]['job']['config'],histories[case,arm][int(round(stamp/.15))],stamp==18)
                assert json.loads(str(z['event_audits']))==[e for e in reports[case,arm]['events'] if e['time']<=stamp+1e-12]
                sampled_checkpoints.append(dict(case=case,arm=arm,**checks))
                if level=='coarse' or level=='fine':zs[level,arm]=z
            ds={key:density(z) for key,z in zs.items()}
            for group,keys in [('backend_rows',[('coarse','gpu'),('coarse','cpu')]),('time_rows',[('coarse','gpu'),('fine','gpu')])]:
                index=('off','tension_adhesion','polarity','full').index(branch)
                row=groups[group][index];values=field_error(ds[keys[0]],ds[keys[1]])
                expected=next(f for f in production[group][index]['fields'] if f['time']==stamp)
                for k,v in values.items():close(v,expected['errors'][k],('native field',branch,stamp,k))
                row['field_checks'].append(dict(time=stamp,errors=values))
                if group=='backend_rows':
                    assert np.array_equal(zs[keys[0]]['ids'],zs[keys[1]]['ids'])
                    assert np.array_equal(zs[keys[0]]['parents'],zs[keys[1]]['parents'])
                    phase=float(abs(zs[keys[0]]['phi']-zs[keys[1]]['phi']).max())
                    row['independent_phase_abs_max']=max(row.get('independent_phase_abs_max',0.),phase)
                    if stamp==18:
                        ma,mb=(json.loads(str(zs[k]['metadata'])) for k in keys)
                        for rng in ('rng','fate_rng','signal_rng'):assert ma[rng]==mb[rng]
                        np.testing.assert_array_equal(zs[keys[0]]['due'],zs[keys[1]]['due'])
            del zs,ds,z
        for level in ('coarse','fine','space88','space104'):
            case=branch+'_'+level
            for stage in (2,4,8):
                z=load(ROOT/case/f'restart-stage-{stage}.npz');checks=checkpoint(z,reports[case,'gpu']['job']['config'])
                assert stage<checks['cells']<2*stage and checks['active']>0
                restart_samples.append(dict(case=case,stage=stage,**checks))
        print('Native fields and checkpoint schema reviewed:',branch,flush=True)
    print('Reconstructing one failed spatial-field witness per pair',flush=True)
    for i,original in enumerate(production['space_rows']):
        witness=max([f for f in original['fields'] if f['probe']==112 and f['order']==3],key=lambda f:max(f['errors'].values()))
        stamp=witness['time'];ds=[]
        for case in original['cases']:
            z=load(ROOT/case/f'gpu-t{stamp:02d}.npz');ds.append(density(z,112))
        values=field_error(*ds)
        for k,v in values.items():close(v,witness['errors'][k],('spatial witness',i,k))
        assert max(values.values())>protocol['space_limits']['field_max']
        groups['space_rows'][i]['field_checks'].append(dict(time=stamp,probe=112,order=3,errors=values))
        del ds,z
        print('Spatial failure independently reproduced:',original['cases'],stamp,flush=True)
    for group,limit_key in [('backend_rows','backend_limits'),('time_rows','time_limits'),('space_rows','space_limits')]:
        for i,row in enumerate(groups[group]):
            field_max=max(v for f in row['field_checks'] for v in f['errors'].values())
            if group!='space_rows':close(field_max,row['production_field_max'],(group,'field max'))
            failures=dict(row['failed_history_gates'])
            if field_max>protocol[limit_key]['field_max']:failures['field_max']=field_max
            if group=='backend_rows':
                close(row['independent_phase_abs_max'],production[group][i]['homologous_errors']['phase_abs_max'],'phase max')
                assert row['independent_phase_abs_max']<=protocol['homologous_limits']['phase_abs_max']
            row['independently_checked_failures']=failures
            row['independently_checked_passed']=not failures
            assert row['independently_checked_passed']==row['production_passed']
    print('Checking preserved accepted evidence',flush=True)
    preserved=read(protocol['preserved_manifest'])
    for path,expected in preserved['preserved_sha256'].items():pin(path,expected)
    result=dict(review_completed=True,qualification_accepted=False,study_passed=False,
        reviewed_at_utc=datetime.now(timezone.utc).isoformat(),protocol_sha256=protocol_hash,
        source_count=len(protocol['source_sha256']),input_count=len(protocol['input_sha256']),
        preserved_count=len(preserved['preserved_sha256']),paths=paths,**groups,
        sampled_checkpoints=sampled_checkpoints,restart_checkpoint_samples=restart_samples,
        independent_full_restart_replay_performed=False,
        unreplicated_production_checks=['All spatial field checkpoints/endpoint reconstruction variants and spatial trend extrema; failure witnesses suffice to confirm nonacceptance.',
            'Every-step physical extrema and short bitwise restart tests are recorded runtime evidence; saved arrays sample those states.'],
        spatial_trend_production=production['spatial_trend_rows'],
        consumed_sha256=WATCHED,script_sha256=pin(__file__),new_scientific_histories=0,
        production_seconds=timing['seconds'],review_seconds=time.perf_counter()-started,
        interpretation='All 20 paths finish and four coarse backend comparisons reproduce. Four temporal and eight spatial comparisons fail unchanged gates; no developmental accuracy or fresh-history acceptance.')
    (HERE/'independent_review.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('review_completed','qualification_accepted','preserved_count','review_seconds')}),flush=True)


if __name__=='__main__':main()
