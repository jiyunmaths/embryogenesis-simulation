"""Read-only visualization adapter for completed nonlinear benchmark artifacts."""

import json
from pathlib import Path
import threading

import numpy as np


class BenchmarkData:
    """Expose named runs from one server-configured directory, never URL paths.

    Loading is lazy and locked. A missing/incomplete benchmark does not prevent
    the live embryo dashboard from starting. Restart after replacing artifacts.
    """

    def __init__(self, directory=None):
        if directory is None:
            candidates = [Path('outputs/nonlinear-bridge-verified'), Path('outputs/nonlinear-bridge')]
            directory = next((p for p in candidates if (p / 'analysis.json').is_file()
                              and (p / 'fields.npz').is_file()), candidates[0])
        self.directory = Path(directory)
        self._lock = threading.RLock()
        self._metadata = None
        self._runs = {}
        self._cached_run = None

    def metadata(self):
        with self._lock:
            if self._metadata is not None:
                return self._metadata
            try:
                report = json.loads((self.directory / 'analysis.json').read_text())
                if report.get('schema') != 1 or not report.get('checks') or report.get('numerical_failure'):
                    raise ValueError('The nonlinear benchmark is incomplete or has an unsupported format')
                runs = {}
                for row in report['spatial_runs']:
                    runs[f'n{row["n"]}'] = {**row, 'kind': 'spatial', 'label': f'Spatial · n={row["n"]} · dt={row["dt"]:g}'}
                finest = report['parameters']['resolutions'][-1]
                for factor, row in zip((2, 4), report['temporal_runs'][1:]):
                    runs[f'temporal_{factor}'] = {**row, 'kind': 'temporal', 'label': f'Time refinement · n={row["n"]} · dt={row["dt"]:g}'}
                row = report['stable_control']
                runs['stable_control'] = {**row, 'kind': 'control', 'label': f'Equal diffusion control · n={row["n"]}'}
                # One fixed color range per chemical across ALL saved runs and
                # times makes the stable control and early amplitudes comparable.
                ranges = {'a': [float('inf'), -float('inf')], 'h': [float('inf'), -float('inf')]}
                with np.load(self.directory / 'fields.npz', allow_pickle=False) as archive:
                    for name, row in runs.items():
                        n = row['n']
                        if type(n) is not int or not 4 <= n <= 32 or n % 2:
                            raise ValueError('Unsupported benchmark mesh')
                        times = archive[f'{name}_time']
                        mask = archive[f'n{n}_mask']
                        if mask.shape != (n, n, n) or mask.dtype != np.bool_:
                            raise ValueError('Invalid benchmark mask')
                        if (times.ndim != 1 or len(times) < 2 or not np.isfinite(times).all()
                                or np.any(np.diff(times) <= 0) or times[0] != 0
                                or not np.isclose(times[-1], row['duration'])):
                            raise ValueError('Invalid saved benchmark times')
                        row['times'] = times.tolist()
                        row['compartments'] = int(mask.sum())
                        for axis in 'xyz':
                            edges = archive[f'n{n}_{axis}_edges']
                            if edges.shape != (n+1,) or not np.isfinite(edges).all() or np.any(np.diff(edges) <= 0):
                                raise ValueError('Invalid benchmark mesh edges')
                        for species in ('a', 'h'):
                            values = archive[f'{name}_{species}']
                            if values.shape != (len(times), int(mask.sum())) or not np.isfinite(values).all() or np.any(values <= 0):
                                raise ValueError('Invalid benchmark concentration fields')
                            ranges[species][0] = min(ranges[species][0], float(values.min()))
                            ranges[species][1] = max(ranges[species][1], float(values.max()))
                self._runs = runs
                self._metadata = {'available': True, 'source': self.directory.name,
                                  'scope': report['scope'], 'parameters': report['parameters'],
                                  'checks': report['checks'], 'acceptance': report['acceptance'],
                                  'spatial_comparisons': report['spatial_comparisons'],
                                  'temporal_comparisons': report['temporal_comparisons'],
                                  'observed_temporal_order': report['observed_temporal_order'],
                                  'default_run': f'n{finest}', 'ranges': ranges,
                                  'runs': [{'id': name, **row} for name, row in runs.items()]}
                # Reject malformed JSON numbers before sending an HTTP response.
                json.dumps(self._metadata, allow_nan=False)
                return self._metadata
            except (OSError, ValueError, KeyError, TypeError) as error:
                self._metadata = None
                return {'available': False, 'source': self.directory.name,
                        'message': 'No readable completed nonlinear benchmark is available.',
                        'detail': str(error),
                        'command': 'python -m embryo.nonlinear --output outputs/nonlinear-bridge',
                        'hint': 'Start the dashboard with --benchmark PATH to read a different completed run.'}

    def run(self, name):
        with self._lock:
            if not self.metadata()['available']:
                raise ValueError('No readable completed nonlinear benchmark is available')
            if name not in self._runs:
                raise KeyError('Unknown benchmark run')
            if self._cached_run is not None and self._cached_run['id'] == name:
                return self._cached_run
            row = self._runs[name]
            n = row['n']
            with np.load(self.directory / 'fields.npz', allow_pickle=False) as archive:
                self._cached_run = {'id': name, 'n': n, 'times': row['times'],
                                    'edges': [archive[f'n{n}_{axis}_edges'].tolist() for axis in 'xyz'],
                                    'indices': np.argwhere(archive[f'n{n}_mask']).tolist(),
                                    'a': archive[f'{name}_a'].tolist(), 'h': archive[f'{name}_h'].tolist()}
            return self._cached_run
