"""Local, live dashboard for the actual 3D simulation.

Run ``python -m embryo.dashboard`` and open the printed loopback URL. The
single simulation worker advances only between Run and Pause; display requests
read immutable snapshots and never advance or perturb the model's RNG streams.
"""

import argparse
from collections import deque
from dataclasses import asdict, fields
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from importlib.resources import files
import json
import math
import os
from pathlib import Path
import threading
import time
from urllib.parse import parse_qs, urlsplit

# Small contact matrices benefit from a single BLAS thread; retain an explicit
# user choice. This takes effect when this module is the command-line entrypoint.
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np

from .graph_analysis import cycle
from .model import Config, Simulation
from .signaling import normalized_graph, stability


MAX_BODY = 64 * 1024
MAX_HISTORY_BYTES = 16 * 1024 * 1024
MAX_HISTORY_FRAMES = 400
CONTROL_TIMEOUT = 5.0


class ConflictError(ValueError):
    """A valid request is incompatible with the current worker state."""


def schema():
    defaults = asdict(Config())
    names = {bool: "boolean", int: "integer", float: "number", str: "string"}
    return {
        "defaults": defaults,
        "types": {field.name: names[field.type] for field in fields(Config)},
        "choices": {"division_orientation": ["shape", "isotropic"]},
        "limits": {"grid": [12, 128], "max_cells": [1, 128],
                   "seed": [0, 2**32 - 1], "cell_grid_voxels": 16_777_216,
                   "dt_times_signal_rate": 100, "history_frames": MAX_HISTORY_FRAMES,
                   "history_json_bytes": MAX_HISTORY_BYTES},
    }


def validated_config(values):
    """Validate JSON primitive types before Config's numerical constraints.

    Local dashboard resource bounds prevent an accidental parameter edit from
    allocating an unbounded dense cell/grid array or a huge signal substep loop.
    The command-line batch simulator remains available for larger experiments.
    """
    if not isinstance(values, dict):
        raise ValueError("config must be a JSON object")
    description = schema()
    unknown = values.keys() - description["defaults"].keys()
    if unknown:
        raise ValueError("unknown configuration fields: " + ", ".join(sorted(unknown)))
    for name, value in values.items():
        kind = description["types"][name]
        if kind == "integer" and type(value) is not int:
            raise ValueError(f"{name} must be an integer")
        if kind == "boolean" and type(value) is not bool:
            raise ValueError(f"{name} must be a boolean")
        if kind == "string" and type(value) is not str:
            raise ValueError(f"{name} must be a string")
        if kind == "number":
            try:
                finite = type(value) in (float, int) and math.isfinite(value)
            except OverflowError:
                finite = False
            if not finite:
                raise ValueError(f"{name} must be a finite number")
    config = Config(**values)
    if not 0 <= config.seed <= 2**32 - 1:
        raise ValueError("seed must be between 0 and 4294967295")
    if config.grid > 128 or config.max_cells > 128:
        raise ValueError("the live dashboard supports grid <= 128 and max_cells <= 128")
    if config.grid**3 * config.max_cells > 16_777_216:
        raise ValueError("grid**3 * max_cells exceeds the dashboard's 16777216-voxel resource limit")
    if config.dt * max(config.signal_beta + config.signal_dh, 1 + config.signal_da) > 100:
        raise ValueError("reduce signal rates or dt: excessive signaling substeps for the live dashboard")
    config.validate()
    return config


def reference_preflight(config, simulation):
    """Reference finite-graph spectra, evaluated before any 3D time step."""
    reports = {}
    for count in (4, 8, 16):
        graphs = {f"cycle_{count}": cycle(count),
                  f"complete_{count}": normalized_graph(np.ones((count, count)) - np.eye(count))}
        for name, graph in graphs.items():
            reports[name] = stability(graph, config.signal_beta, config.signal_da, config.signal_dh)
    return {"parameters": {"beta": config.signal_beta, "da": config.signal_da,
                           "dh": config.signal_dh, "cutoff": config.graph_contact_cutoff},
            "graphs": reports, "initial_graph": simulation.graph_snapshot(),
            "signaling_enabled": config.signaling,
            "interpretation": "Reference cycle/complete graphs are not predicted embryo contacts. "
                              "Frozen-graph growth rates apply near the homogeneous equilibrium; "
                              "normalized eigenvalues are not physical wavenumbers."}


class DashboardController:
    """One worker, serialized controls, and immutable published frames.

    A worker owns its Simulation while running. HTTP readers access only cached
    metadata/frames. Pausing and resuming preserves the very same Simulation and
    random generators; resetting creates a new generation only after validation
    and a safe step boundary. No checkpoint files are written implicitly.
    """

    def __init__(self, config=None):
        self._lock = threading.RLock()
        self._control = threading.Lock()
        self._stop = threading.Event()
        self._worker = None
        self._generation = -1
        self._closed = False
        self._replace(validated_config(asdict(config) if isinstance(config, Config)
                                       else ({} if config is None else config)))

    def _make_frame(self, simulation, revision):
        frame = {"revision": revision, "step": simulation.step_number,
                 "metrics": simulation.metrics(), "cells": simulation.surfaces(max_points=1000),
                 "graph": simulation.graph_snapshot()}
        # Refuse to publish a corrupted/non-finite state as a successful frame.
        size = len(json.dumps(frame, allow_nan=False, separators=(",", ":")))
        return frame, size

    def _replace(self, config):
        simulation = Simulation(config)
        preflight = reference_preflight(config, simulation)
        frame, size = self._make_frame(simulation, 0)
        with self._lock:
            self._simulation = simulation
            self._config = asdict(config)
            self._generation += 1
            self._revision = 0
            self._frames = deque([(frame, size)])
            self._history_bytes = size
            self._preflight = {"generation": self._generation, **preflight}
            self._step = 0
            self._state = "ready"
            self._error = None
            self._elapsed = 0.0
            self._started_at = None
            self._stop.clear()

    def _publish(self, simulation):
        with self._lock:
            if self._frames[-1][0]["step"] == simulation.step_number:
                return
            revision = self._revision + 1
        frame, size = self._make_frame(simulation, revision)
        with self._lock:
            self._revision = revision
            self._frames.append((frame, size))
            self._history_bytes += size
            while len(self._frames) > 1 and (len(self._frames) > MAX_HISTORY_FRAMES or
                                             self._history_bytes > MAX_HISTORY_BYTES):
                self._history_bytes -= self._frames.popleft()[1]

    def snapshot(self, after=-1):
        with self._lock:
            elapsed = self._elapsed
            if self._started_at is not None:
                elapsed += time.monotonic() - self._started_at
            # If a reset took place after the caller's last revision, send its
            # new history. The generation number disambiguates revision reuse.
            if after > self._revision:
                after = -1
            return {"state": self._state, "generation": self._generation,
                    "revision": self._revision, "step": self._step,
                    "total_steps": self._config["steps"], "config": dict(self._config),
                    "error": self._error, "elapsed_seconds": elapsed,
                    "history_start_revision": self._frames[0][0]["revision"],
                    "frames": [frame for frame, _ in self._frames if frame["revision"] > after]}

    def preflight(self):
        with self._lock:
            return self._preflight

    def _work(self):
        simulation = self._simulation
        try:
            while simulation.step_number < simulation.config.steps and not self._stop.is_set():
                simulation.step()
                with self._lock:
                    self._step = simulation.step_number
                if simulation.step_number % simulation.config.save_every == 0:
                    self._publish(simulation)
            self._publish(simulation)
            with self._lock:
                self._state = "completed" if simulation.step_number >= simulation.config.steps else "paused"
        except Exception as error:
            # Preserve the latest valid frame when a failed step has left
            # non-finite arrays. Error status is authoritative over that frame.
            try:
                self._publish(simulation)
            except Exception:
                pass
            with self._lock:
                self._step = simulation.step_number
                self._state = "error"
                self._error = f"{type(error).__name__}: {error}"
        finally:
            with self._lock:
                self._elapsed += time.monotonic() - self._started_at
                self._started_at = None

    def _join_worker(self):
        self._stop.set()
        if self._worker is not None:
            self._worker.join(CONTROL_TIMEOUT)
            if self._worker.is_alive():
                raise ConflictError("The current numerical step is still finishing. Pause is requested; "
                                    "wait for paused status and retry.")

    def run(self, values=None):
        candidate = validated_config(values) if values is not None else None
        with self._control:
            if self._closed:
                raise ConflictError("dashboard is shutting down")
            with self._lock:
                state, current = self._state, self._config
            if state == "running":
                if candidate is not None and asdict(candidate) != current:
                    raise ConflictError("Reset before changing parameters of a running simulation")
                return self.snapshot()
            if state not in ("ready", "paused"):
                raise ConflictError("Reset before starting a completed or failed simulation")
            if state == "paused" and candidate is not None and asdict(candidate) != current:
                raise ConflictError("Reset before changing parameters of a paused simulation")
            # A worker may have published its final state immediately before
            # its finally block; join before starting another elapsed timer.
            self._join_worker()
            if state == "ready" and candidate is not None and asdict(candidate) != current:
                self._replace(candidate)
            self._stop.clear()
            with self._lock:
                self._state = "running"
                self._started_at = time.monotonic()
                self._worker = threading.Thread(target=self._work, name="embryo-simulation", daemon=True)
                self._worker.start()
            return self.snapshot()

    def pause(self):
        with self._control:
            self._join_worker()
            return self.snapshot()

    def reset(self, values):
        candidate = validated_config(values)  # Never interrupt a run for invalid edits.
        with self._control:
            if self._closed:
                raise ConflictError("dashboard is shutting down")
            self._join_worker()
            self._replace(candidate)
            return self.snapshot()

    def close(self):
        with self._control:
            self._closed = True
            self._join_worker()


def make_server(controller, port=8765):
    """Bind IPv4 loopback only; expose a fixed, package-owned static allowlist."""

    class Handler(BaseHTTPRequestHandler):
        server_version = "EmbryoDashboard/1"

        def log_message(self, format, *args):
            # Frequent state polling should not obscure numerical/HTTP errors.
            if len(args) > 1 and str(args[1]) not in ("200", "304"):
                super().log_message(format, *args)

        def _send(self, status, data, content_type="application/json; charset=utf-8"):
            if isinstance(data, (dict, list)):
                data = json.dumps(data, allow_nan=False, separators=(",", ":")).encode()
            self.send_response(status)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", str(len(data)))
            self.send_header("Cache-Control", "no-store")
            self.send_header("X-Content-Type-Options", "nosniff")
            self.send_header("Content-Security-Policy", "default-src 'self'; script-src 'self'; "
                             "style-src 'self'; img-src 'self' data:; connect-src 'self'; "
                             "frame-ancestors 'none'; base-uri 'none'; form-action 'none'")
            self.end_headers()
            self.wfile.write(data)

        def _local_request(self):
            port = self.server.server_address[1]
            allowed = {f"localhost:{port}", f"127.0.0.1:{port}"}
            if port == 80:
                allowed.update(("localhost", "127.0.0.1"))
            hosts = self.headers.get_all("Host", [])
            if len(hosts) != 1 or hosts[0].lower() not in allowed:
                self._send(403, {"error": "Only local dashboard Host headers are accepted"})
                return False
            origin = self.headers.get("Origin")
            if origin is not None and origin.lower() != "http://" + hosts[0].lower():
                self._send(403, {"error": "Cross-origin dashboard requests are not accepted"})
                return False
            if self.headers.get("Sec-Fetch-Site") == "cross-site":
                self._send(403, {"error": "Cross-site dashboard requests are not accepted"})
                return False
            return True

        def do_GET(self):
            if not self._local_request():
                return
            url = urlsplit(self.path)
            if url.path == "/api/schema":
                self._send(200, schema())
            elif url.path == "/api/preflight":
                self._send(200, controller.preflight())
            elif url.path == "/api/state":
                try:
                    after = int(parse_qs(url.query).get("after", ["-1"])[0])
                except ValueError:
                    self._send(400, {"error": "after must be an integer revision"})
                    return
                self._send(200, controller.snapshot(after))
            elif url.path in ("/", "/dashboard.js", "/dashboard.css"):
                filename = "dashboard.html" if url.path == "/" else url.path[1:]
                content_type = {"dashboard.html": "text/html; charset=utf-8",
                                "dashboard.js": "text/javascript; charset=utf-8",
                                "dashboard.css": "text/css; charset=utf-8"}[filename]
                try:
                    self._send(200, files("embryo").joinpath(filename).read_bytes(), content_type)
                except FileNotFoundError:
                    self._send(404, {"error": "Dashboard asset is missing from the package"})
            else:
                self._send(404, {"error": "Unknown dashboard route"})

        def do_POST(self):
            if not self._local_request():
                return
            route = urlsplit(self.path).path
            if route not in ("/api/run", "/api/pause", "/api/reset"):
                self._send(404, {"error": "Unknown dashboard route"})
                return
            if self.headers.get("Transfer-Encoding"):
                self._send(400, {"error": "Transfer-Encoding is not supported"})
                return
            if self.headers.get_content_type() != "application/json":
                self._send(415, {"error": "Use Content-Type: application/json"})
                return
            lengths = self.headers.get_all("Content-Length", [])
            try:
                if len(lengths) != 1:
                    raise ValueError
                size = int(lengths[0])
                if size < 0:
                    raise ValueError
            except ValueError:
                self._send(400, {"error": "A valid Content-Length is required"})
                return
            if size > MAX_BODY:
                self._send(413, {"error": "Request body exceeds 64 KiB"})
                return
            try:
                self.connection.settimeout(5)
                body = json.loads(self.rfile.read(size))
                if not isinstance(body, dict):
                    raise ValueError("request body must be a JSON object")
                allowed = {"config"} if route != "/api/pause" else set()
                if body.keys() - allowed:
                    raise ValueError("unknown request fields")
                if "config" in body and not isinstance(body["config"], dict):
                    raise ValueError("config must be a JSON object")
                if route == "/api/reset" and "config" not in body:
                    raise ValueError("reset requires config")
                if route == "/api/run":
                    result = controller.run(body.get("config"))
                elif route == "/api/pause":
                    result = controller.pause()
                else:
                    result = controller.reset(body["config"])
                self._send(200, result)
            except ConflictError as error:
                self._send(409, {"error": str(error)})
            except (ValueError, TypeError, UnicodeDecodeError, OverflowError) as error:
                self._send(400, {"error": str(error)})
            except TimeoutError:
                self._send(408, {"error": "Timed out reading request body"})
            except Exception as error:
                self._send(500, {"error": f"{type(error).__name__}: {error}"})

    return ThreadingHTTPServer(("127.0.0.1", port), Handler)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--port", type=int, default=8765, help="Loopback HTTP port (default: 8765)")
    parser.add_argument("--config", type=Path, help="Initial model configuration JSON")
    args = parser.parse_args()
    if not 0 <= args.port <= 65535:
        parser.error("port must be between 0 and 65535")
    try:
        values = json.loads(args.config.read_text()) if args.config else {}
        controller = DashboardController(validated_config(values))
        server = make_server(controller, args.port)
    except (ValueError, TypeError, OSError) as error:
        parser.error(str(error))
    print(f"Live embryo dashboard: http://127.0.0.1:{server.server_address[1]}", flush=True)
    print("The server is local to this computer. Ctrl+C stops it. Runs are kept in memory.", flush=True)
    try:
        server.serve_forever(poll_interval=.2)
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
        try:
            controller.close()
        except ConflictError:
            # Worker is a daemon; no files or checkpoint writes can be left half-finished.
            pass


if __name__ == "__main__":
    main()
