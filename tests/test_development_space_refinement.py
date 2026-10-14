"""Guard reuse, model matching and complete finer-grid sampling before launch."""
import json
from pathlib import Path
import pytest

from embryo import development_space_refinement as s


def test_reused_and_symlinked_paths_cannot_be_executed(tmp_path,monkeypatch):
    def forbidden(*args): raise AssertionError('Original evidence would be overwritten')
    monkeypatch.setattr(s.td,'run_path',forbidden)
    with pytest.raises(ValueError,match='read-only'):
        s.execute_path(tmp_path,dict(key='old',reused=True))
    original = tmp_path/'original'; original.mkdir()
    (tmp_path/'link').symlink_to(original,target_is_directory=True)
    with pytest.raises(ValueError,match='read-only'):
        s.execute_path(tmp_path,dict(key='link',reused=False))


def test_preparation_refuses_an_existing_study(tmp_path):
    marker = tmp_path/'old-result'; marker.write_text('preserved')
    with pytest.raises(FileExistsError): s.prepare(tmp_path)
    assert marker.read_text() == 'preserved'


def test_source_matching_rejects_physical_parameter_and_clock_drift():
    base = dict(grid=72,dt=.00046875,steps=38400,save_every=320,feedback=False,
        fate_tension=0.,fate_adhesion=0.,polarity_tension=0.,signal_dh=.55)
    candidate = dict(base,grid=120,dt=.000234375,steps=76800,save_every=640)
    assert set(s.config_changes(base,candidate)) == {'grid','dt','steps','save_every'}
    with pytest.raises(ValueError,match='Physical model'): s.config_changes(base,dict(candidate,signal_dh=.56))
    with pytest.raises(ValueError,match='clock'): s.config_changes(base,dict(candidate,save_every=320))


def test_reuse_report_must_keep_its_original_protocol_and_config(tmp_path):
    original = tmp_path/'preserved'; original.mkdir(); (tmp_path/'off72').symlink_to(original)
    original_job = dict(key='off_quarter',config={'grid':72})
    report = dict(job=original_job,protocol_sha256='original',evidence_sha256={})
    path = original/'gpu-result.json'; path.write_text(json.dumps(report))
    provenance = dict(folder=str(original.resolve()),report_sha256=s.digest(path),
        job=original_job,protocol_sha256='original')
    job = dict(key='off72',config={'grid':72},reused=True)
    assert s.case_report(tmp_path,job,dict(reuse_provenance=provenance)) == report
    with pytest.raises(ValueError,match='configuration'):
        s.case_report(tmp_path,dict(job,config={'grid':88}),dict(reuse_provenance=provenance))
    with pytest.raises(ValueError,match='provenance'):
        s.case_report(tmp_path,job,dict(reuse_provenance=dict(provenance,protocol_sha256='new')))


def network(time,n=1):
    return dict(time=time,volumes=[1.]*n,delta=[[0.]*n for _ in range(n)],ids=list(range(n)),
        parents=[-1]*n,conductance=[[0.]*n for _ in range(n)],edges=0)


def write_network(folder,rows):
    folder.mkdir(); (folder/'network.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in rows))


def test_dense_temporal_check_reports_count_gaps_without_fabricating_spectra(tmp_path):
    a = [network(i*.015) for i in range(1201)]; b = [dict(r) for r in a]
    b[17] = network(17*.015,2)
    write_network(tmp_path/'a',a);write_network(tmp_path/'b',b)
    result = s.dense_pair(tmp_path,'a','b')
    assert result['unmatched_times'] == [.255] and result['peak']['relative_spectrum'] == 0
    b.pop(18); (tmp_path/'b/network.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in b))
    with pytest.raises(ValueError,match='Incomplete'): s.dense_pair(tmp_path,'a','b')


@pytest.mark.parametrize('failure',['physical','spatial','time','dense','trend'])
def test_scoped_pass_cannot_bypass_a_failed_required_screen(failure):
    paths = [{'passed':True}]; space = [{'passed':False},{'passed':True}]
    time = {'passed':True}; dense = {'peak':{'relative_spectrum':.0099}}; trend=True
    # Earlier coarse failures are kept; all finest-pair and temporal gates still apply.
    assert s.scoped_decision(paths,space,time,dense,trend,.01)
    if failure == 'physical': paths[0]['passed']=False
    if failure == 'spatial': space[-1]['passed']=False
    if failure == 'time': time['passed']=False
    if failure == 'dense': dense['peak']['relative_spectrum']=.0101
    if failure == 'trend': trend=False
    assert not s.scoped_decision(paths,space,time,dense,trend,.01)
