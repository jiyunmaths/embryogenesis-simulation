"""Read-only dashboard views must faithfully expose benchmark artifacts."""

import http.client
import json
import shutil
import threading

import numpy as np
import pytest

from embryo.dashboard import DashboardController, make_server
from embryo.dashboard_benchmark import BenchmarkData
from embryo.model import Config
from embryo.nonlinear import run_benchmark


@pytest.fixture(scope='module')
def benchmark_path(tmp_path_factory):
    path = tmp_path_factory.mktemp('dashboard-benchmark') / 'recorded'
    run_benchmark(path, resolutions=(4, 8, 16), duration=1., dt=.05, sample_interval=.25)
    return path


def test_all_named_runs_and_fixed_scales_match_source(benchmark_path):
    reader = BenchmarkData(benchmark_path)
    meta = reader.metadata()
    assert meta['available']
    assert meta['default_run'] == 'n16'
    assert [r['id'] for r in meta['runs']] == ['n4', 'n8', 'n16', 'temporal_2', 'temporal_4', 'stable_control']
    report = json.loads((benchmark_path / 'analysis.json').read_text())
    assert meta['checks'] == report['checks']
    assert not meta['checks']['persistent_fine_mesh_contrast']  # Preserve failed science checks.
    with np.load(benchmark_path / 'fields.npz', allow_pickle=False) as archive:
        ranges = {s: [min(archive[f'{r["id"]}_{s}'].min() for r in meta['runs']),
                      max(archive[f'{r["id"]}_{s}'].max() for r in meta['runs'])] for s in ('a', 'h')}
        assert meta['ranges'] == ranges
        for run in meta['runs']:
            result = reader.run(run['id'])
            np.testing.assert_array_equal(result['a'], archive[f'{run["id"]}_a'])
            np.testing.assert_array_equal(result['h'], archive[f'{run["id"]}_h'])
            np.testing.assert_array_equal(result['times'], archive[f'{run["id"]}_time'])
            np.testing.assert_array_equal(result['indices'], np.argwhere(archive[f'n{run["n"]}_mask']))
            assert result is reader.run(run['id'])
            json.dumps(result, allow_nan=False)
    with pytest.raises(KeyError):
        reader.run('../../README.md')


def test_missing_data_can_become_available_without_resetting_embryo(tmp_path, benchmark_path):
    path = tmp_path / 'later'
    reader = BenchmarkData(path)
    assert not reader.metadata()['available']
    shutil.copytree(benchmark_path, path)
    assert reader.metadata()['available']


@pytest.mark.parametrize('damage', ['incomplete', 'nonfinite', 'shape', 'times'])
def test_invalid_or_partial_benchmarks_are_explicitly_unavailable(tmp_path, benchmark_path, damage):
    destination = tmp_path / 'bad'
    shutil.copytree(benchmark_path, destination)
    if damage == 'incomplete':
        report = json.loads((destination / 'analysis.json').read_text())
        report.pop('checks')
        (destination / 'analysis.json').write_text(json.dumps(report))
    else:
        with np.load(destination / 'fields.npz', allow_pickle=False) as archive:
            fields = {key: archive[key].copy() for key in archive.files}
        if damage == 'nonfinite':
            fields['n4_a'][0, 0] = np.nan
        elif damage == 'shape':
            fields['n4_a'] = fields['n4_a'][:, :-1]
        else:
            fields['n4_time'][-1] = 0
        np.savez_compressed(destination / 'fields.npz', **fields)
    reader = BenchmarkData(destination)
    assert not reader.metadata()['available']
    with pytest.raises(ValueError, match='No readable'):
        reader.run('n4')


def test_http_benchmark_routes_do_not_change_live_simulation(benchmark_path):
    controller = DashboardController(Config(grid=12, interface_width=.18, max_cells=1, competence_cells=1))
    server = make_server(controller, 0, benchmark_path)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    def request(route):
        conn = http.client.HTTPConnection('127.0.0.1', server.server_address[1], timeout=5)
        conn.request('GET', route)
        response = conn.getresponse()
        raw = response.read()
        result = json.loads(raw) if 'application/json' in response.getheader('Content-Type') else raw
        conn.close()
        return response.status, result

    try:
        before = controller.snapshot()
        assert request('/api/benchmark')[1]['available']
        assert request('/api/benchmark/run?id=temporal_2')[1]['id'] == 'temporal_2'
        assert request('/api/benchmark/run?id=../../README.md')[0] == 404
        assert request('/api/benchmark/run')[0] == 404
        assert request('/benchmark.js')[0] == 200
        assert b'Continuum signaling' in request('/')[1]
        after = controller.snapshot()
        assert before['frames'] == after['frames']
        assert after['state'] == 'ready' and after['step'] == 0
        assert before['generation'] == after['generation']
    finally:
        server.shutdown()
        server.server_close()
        thread.join(5)
        controller.close()
