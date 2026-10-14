"""Protect independent history accounting, matched starts and hard numerical gates."""
from dataclasses import asdict
import json
from types import SimpleNamespace

import numpy as np
import pytest

from embryo import fresh_zygote_pilot as f
from embryo.model import Config
from embryo.qualify_development_gpu import save_cpu


def baseline():
    return asdict(Config(grid=120,extent=2.24,interface_width=.085,dt=.00046875,
        steps=38400,save_every=320,signal_dh=.55,signal_partition_noise=.001,
        differentiation=False,fate_rate=0.,fate_noise=0.,partition_noise=0.,
        feedback=False,fate_tension=0.,fate_adhesion=0.,polarity_tension=0.))


def test_three_new_histories_four_controls_two_exact_clocks():
    jobs=f.design(baseline())
    assert len(jobs)==24
    assert {j['seed'] for j in jobs}=={1101,1102,1103}
    assert len({(j['seed'],j['branch'],j['level']) for j in jobs})==24
    for job in jobs:
        c=job['config'];f.config_changes(baseline(),c,job['branch'])
        assert c['grid']==120 and c['steps']*c['dt']==18
        assert c['save_every']*c['dt']==.15
        assert c['signal_da']==.02 and c['signal_dh']==.55
        assert c['signal_partition_noise']==.001
    assert [j['level'] for j in jobs[:2]]==['development','development_half']


@pytest.mark.parametrize('change', [dict(signal_dh=.56),dict(grid=104),dict(seed=7),
    dict(dt=.0009375),dict(save_every=1),dict(polarity_tension=.7)])
def test_unregistered_physics_history_clock_or_coupling_is_rejected(change):
    job=f.design(baseline())[-2]
    with pytest.raises(ValueError):
        f.config_changes(baseline(),dict(job['config'],**change),job['branch'])


def source_jobs(tmp_path):
    jobs=f.design(baseline())
    for j in jobs:
        # Small analytic arrays exercise matching/schema without a CUDA allocation.
        c=Config(**dict(j['config'],grid=32));host=f.q.initial(c)
        j['config']=asdict(host.config);source=tmp_path/(j['key']+'.npz')
        save_cpu(host,source);j['source']=str(source)
    return jobs


def change_npz(path,name,value):
    with np.load(path,allow_pickle=False) as z: data={k:z[k].copy() for k in z.files}
    data[name]=np.asarray(value);np.savez_compressed(path,**data)


def test_matching_includes_every_state_array_and_all_random_states(tmp_path):
    jobs=source_jobs(tmp_path);signatures=f.match_starts(jobs)
    assert len(signatures)==24 and len(set(signatures.values()))==3
    with np.load(jobs[1]['source']) as z:due=z['due'].copy()
    change_npz(jobs[1]['source'],'due',due+.01)
    with pytest.raises(ValueError,match='not matched'):
        f.match_starts(jobs)


@pytest.mark.parametrize('field,value', [('fate',[1.]),('activator',[1.04]),
    ('phase_carry',np.ones((1,32,32,32))),('chemical_method','wrong method')])
def test_a_label_pattern_residual_or_wrong_schema_is_not_a_fresh_start(tmp_path,field,value):
    jobs=source_jobs(tmp_path);change_npz(jobs[0]['source'],field,value)
    with pytest.raises(ValueError):f.initial_signature(jobs[0]['source'],jobs[0])


def test_distinct_seed_numbers_with_copied_random_states_are_rejected(tmp_path):
    jobs=source_jobs(tmp_path)
    first=next(j for j in jobs if j['seed']==1101)
    with np.load(first['source']) as z:
        meta=json.loads(str(z['metadata']));due=z['due'].copy()
    for j in jobs:
        with np.load(j['source']) as z:current=json.loads(str(z['metadata']))
        for key in ('rng','fate_rng','signal_rng'):current[key]=meta[key]
        change_npz(j['source'],'metadata',json.dumps(current))
        change_npz(j['source'],'due',due)
    with pytest.raises(ValueError,match='distinct'):
        f.match_starts(jobs)


def test_seed_reservation_checks_original_bytes_and_used_identifiers(tmp_path):
    file=tmp_path/'old-protocol.json';file.write_text('{"seed":7}\n')
    inventory=dict(proposed_histories=list(f.HISTORIES),all_proposed_unused=True,
        used_seed_or_history_ids=[7],files=[dict(path=str(file),sha256=f.digest(file))])
    f.verify_seed_inventory(inventory)
    with pytest.raises(ValueError,match='already used'):
        f.verify_seed_inventory(dict(inventory,used_seed_or_history_ids=[1101]))
    file.write_text('{"seed":8}\n')
    with pytest.raises(ValueError,match='input changed'):
        f.verify_seed_inventory(inventory)


def test_existing_data_and_recorded_prefix_are_never_overwritten(tmp_path):
    marker=tmp_path/'preserved';marker.write_text('original')
    with pytest.raises(FileExistsError):f.prepare(tmp_path,tmp_path/'missing')
    folder=tmp_path/'case';folder.mkdir();(folder/'actual-start-prefix.json').write_text('original')
    with pytest.raises(FileExistsError):f.prefix(tmp_path,dict(key='case'),{})
    assert marker.read_text()=='original'


def test_all_twelve_unique_history_branch_pairs_are_required():
    pairs=[dict(seed=seed,branch=branch,passed=True) for seed in f.HISTORIES for branch in f.BRANCHES]
    assert f.complete_design(pairs)
    assert not f.complete_design(pairs[:-1])
    repeated=list(pairs);repeated[-1]=pairs[0];assert not f.complete_design(repeated)
    failed=list(pairs);failed[5]=dict(failed[5],passed=False);assert not f.complete_design(failed)


@pytest.mark.parametrize('stage',['prefix','path','refinement'])
def test_failed_gate_stops_next_case_without_long_window_launch(tmp_path,monkeypatch,stage):
    binary=tmp_path/'binary';binary.write_bytes(b'unchanged kernel')
    jobs=[dict(key='first',reused=False),dict(key='second',reused=False)]
    p=dict(jobs=jobs,time_pairs=[dict(seed=1101,branch='off',cases=['first','second'])],
        runtime_binary_sha256=f.digest(binary),preserved_manifest='unused')
    f.write(tmp_path/'protocol.json',p)
    monkeypatch.setattr(f,'validate_prelaunch',lambda *args:None)
    monkeypatch.setattr(f.td.gpu,'library',lambda:SimpleNamespace(binary_path=str(binary)))
    monkeypatch.setattr(f.torch.cuda,'get_device_name',lambda *args:'dedicated fixture')
    calls=[]
    def prefix(*args):
        calls.append(('prefix',args[1]['key']));return dict(passed=stage!='prefix')
    def run_path(*args):
        calls.append(('path',args[1]['key']));return dict(passed=stage!='path',seconds=1.)
    monkeypatch.setattr(f,'prefix',prefix);monkeypatch.setattr(f.td,'run_path',run_path)
    monkeypatch.setattr(f,'assess_pair',lambda *args:dict(passed=stage!='refinement'))
    with pytest.raises(RuntimeError):f.run(tmp_path)
    status=f.read(tmp_path/'status.json')
    assert status['state']=='failed' and not status['post_cleavage_started']
    if stage in ('prefix','path'):assert all(case=='first' for _,case in calls)
    else:assert len(calls)==4
    assert not (tmp_path/'summary.json').exists()
