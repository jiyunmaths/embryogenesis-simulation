"""The dashboard advances the real model and does not change its trajectory."""

from dataclasses import asdict
import http.client
import json
import threading
import time

import numpy as np
import pytest

from embryo.dashboard import (ConflictError, DashboardController, make_server,
                              schema, validated_config)
from embryo.model import Config, Simulation


def small_config(**changes):
    return Config(grid=12, interface_width=.18, max_cells=1, steps=80,
                  save_every=10, competence_cells=1, fate_noise=.02, **changes)


def until(predicate, timeout=5):
    deadline = time.monotonic() + timeout
    while not predicate():
        if time.monotonic() > deadline:
            raise AssertionError("worker did not reach expected state")
        time.sleep(.002)


@pytest.fixture
def controller():
    control = DashboardController(small_config())
    yield control
    control.close()


def test_schema_and_reference_preflight(controller):
    description = schema()
    assert description["defaults"] == asdict(Config())
    assert description["types"]["signaling"] == "boolean"
    assert description["choices"]["division_orientation"] == ["shape", "isotropic"]
    report = controller.preflight()
    assert report["graphs"]["cycle_4"]["unstable_modes"] == []
    assert len(report["graphs"]["cycle_8"]["unstable_modes"]) == 2
    assert len(report["graphs"]["cycle_16"]["unstable_modes"]) == 4
    assert report["initial_graph"]["ids"] == [0]
    initial = controller.snapshot()
    assert initial["state"] == "ready"
    assert initial["revision"] == initial["step"] == 0
    assert initial["frames"][0]["cells"][0]["id"] == 0
    assert initial["frames"][0]["cells"][0]["points"]


def test_pause_resume_is_identical_to_uninterrupted_model(controller):
    controller.run()
    until(lambda: controller.snapshot()["step"] >= 4)
    paused = controller.pause()
    assert paused["state"] == "paused"
    assert 4 <= paused["step"] < paused["total_steps"]
    assert paused["frames"][-1]["step"] == paused["step"]
    reference = Simulation(small_config())
    for _ in range(paused["step"]):
        reference.step()
    np.testing.assert_array_equal(controller._simulation.phi, reference.phi)
    np.testing.assert_array_equal(controller._simulation.fate, reference.fate)
    assert controller._simulation.fate_rng.bit_generator.state == reference.fate_rng.bit_generator.state
    assert controller.snapshot(paused["revision"])["frames"] == []
    elapsed = paused["elapsed_seconds"]
    assert controller.snapshot()["elapsed_seconds"] == elapsed
    controller.run()
    until(lambda: controller.snapshot()["state"] == "completed")
    controller.pause()  # Also join the worker's final timing block.
    for _ in range(small_config().steps - paused["step"]):
        reference.step()
    np.testing.assert_array_equal(controller._simulation.phi, reference.phi)
    np.testing.assert_array_equal(controller._simulation.fate, reference.fate)
    np.testing.assert_array_equal(controller._simulation.polarity, reference.polarity)
    assert controller._simulation.rng.bit_generator.state == reference.rng.bit_generator.state
    final = controller.snapshot(paused["revision"])
    assert final["frames"][-1]["step"] == 80
    assert all(frame["revision"] > paused["revision"] for frame in final["frames"])
    with pytest.raises(ConflictError, match="Reset"):
        controller.run()


@pytest.mark.parametrize("changes, message", [
    ({"seed": -1}, "seed"), ({"seed": True}, "integer"),
    ({"steps": 1.0}, "integer"), ({"signaling": 1}, "boolean"),
    ({"dt": float("nan")}, "finite"), ({"grid": 1000}, "grid"),
    ({"max_cells": 1000}, "max_cells"), ({"grid": 128, "max_cells": 128}, "voxel"),
    ({"division_orientation": "random"}, "division_orientation"),
    ({"signal_dh": 1e8}, "substeps"), ({"unknown": 1}, "unknown"),
])
def test_invalid_reset_preserves_current_run(controller, changes, message):
    original = controller.snapshot()
    with pytest.raises(ValueError, match=message):
        controller.reset({**original["config"], **changes})
    assert controller.snapshot() == original


def test_reset_isolates_generation_and_run_applies_initial_edits(controller):
    original = controller.snapshot()
    edited = {**original["config"], "seed": 31}
    controller.run(edited)
    until(lambda: controller.snapshot()["step"] >= 3)
    controller.pause()
    state = controller.snapshot()
    assert state["generation"] == original["generation"] + 1
    with pytest.raises(ConflictError, match="Reset"):
        controller.run({**edited, "seed": 42})
    controller.run()
    reset = controller.reset({**edited, "seed": 42})
    assert reset["generation"] == state["generation"] + 1
    assert reset["state"] == "ready"
    assert reset["step"] == reset["revision"] == 0
    assert len(reset["frames"]) == 1
    assert reset["config"]["seed"] == 42
    assert controller.preflight()["generation"] == reset["generation"]
    assert not controller._worker.is_alive()
    assert controller.snapshot(after=100)["frames"] == reset["frames"]


def test_failed_step_retains_a_valid_frame_and_reports_error(controller, monkeypatch):
    def fail():
        raise FloatingPointError("deliberate numerical failure")
    monkeypatch.setattr(controller._simulation, "step", fail)
    controller.run()
    until(lambda: controller.snapshot()["state"] == "error")
    state = controller.pause()
    assert "deliberate numerical failure" in state["error"]
    assert state["frames"][-1]["step"] == 0
    with pytest.raises(ConflictError):
        controller.run()
    assert controller.reset(state["config"])["state"] == "ready"


def test_slow_step_pause_is_bounded_and_invalid_reset_does_not_interrupt(controller, monkeypatch):
    entered, release = threading.Event(), threading.Event()
    original_step = controller._simulation.step

    def slow_step():
        entered.set()
        if not release.wait(2):
            raise RuntimeError("test step was not released")
        original_step()

    monkeypatch.setattr(controller._simulation, "step", slow_step)
    monkeypatch.setattr("embryo.dashboard.CONTROL_TIMEOUT", .01)
    try:
        controller.run()
        assert entered.wait(2)
        generation = controller.snapshot()["generation"]
        with pytest.raises(ValueError):
            controller.reset({"grid": 1})
        assert not controller._stop.is_set()
        assert controller.snapshot()["state"] == "running"
        with pytest.raises(ConflictError, match="still finishing"):
            controller.pause()
        assert controller._stop.is_set()
        assert controller.snapshot()["generation"] == generation
        release.set()
        until(lambda: controller.snapshot()["state"] == "paused")
        assert controller.snapshot()["step"] == 1
    finally:
        release.set()


def test_history_is_bounded_and_reports_first_retained_revision(controller, monkeypatch):
    monkeypatch.setattr("embryo.dashboard.MAX_HISTORY_FRAMES", 3)
    controller.run()
    until(lambda: controller.snapshot()["state"] == "completed")
    state = controller.snapshot()
    assert len(state["frames"]) == 3
    assert state["history_start_revision"] == state["frames"][0]["revision"] > 0
    assert state["frames"][-1]["step"] == 80


@pytest.fixture
def http_dashboard(controller):
    server = make_server(controller, 0)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    yield server.server_address[1]
    server.shutdown()
    server.server_close()
    thread.join(5)


def request(port, method, route, body=None, headers=None):
    connection = http.client.HTTPConnection("127.0.0.1", port, timeout=5)
    raw = json.dumps(body) if body is not None else None
    outgoing = {"Content-Type": "application/json", **(headers or {})}
    connection.request(method, route, raw, outgoing)
    response = connection.getresponse()
    data = response.read()
    status, content_type = response.status, response.getheader("Content-Type")
    connection.close()
    return status, json.loads(data) if "application/json" in content_type else data


def test_http_routes_controls_and_validation(http_dashboard):
    port = http_dashboard
    assert request(port, "GET", "/api/schema")[1]["defaults"]["grid"] == 40
    status, initial = request(port, "GET", "/api/state?after=-1")
    assert status == 200 and initial["state"] == "ready"
    assert request(port, "GET", "/api/state?after=invalid")[0] == 400
    assert request(port, "GET", "/api/preflight")[1]["initial_graph"]["ids"] == [0]
    status, page = request(port, "GET", "/")
    assert status == 200 and b"<html" in page
    assert request(port, "GET", "/dashboard.css")[0] == 200
    assert request(port, "POST", "/api/run", {})[0] == 200
    status, paused = request(port, "POST", "/api/pause", {})
    assert status == 200 and paused["state"] in ("paused", "completed")
    assert request(port, "POST", "/api/reset", {"config": {"grid": 1}})[0] == 400
    assert request(port, "POST", "/api/reset", {})[0] == 400
    assert request(port, "POST", "/api/run", {"config": None})[0] == 400
    assert request(port, "POST", "/api/pause", {"unexpected": 1})[0] == 400
    status, reset = request(port, "POST", "/api/reset", {"config": initial["config"]})
    assert status == 200 and reset["state"] == "ready" and reset["step"] == 0
    assert request(port, "GET", "/../../pyproject.toml")[0] == 404
    assert request(port, "GET", "/README.md")[0] == 404
    assert request(port, "POST", "/unknown", {})[0] == 404


def test_http_rejects_foreign_origins_and_oversized_or_malformed_bodies(http_dashboard):
    port = http_dashboard
    assert request(port, "POST", "/api/run", {}, {"Origin": "https://example.com"})[0] == 403
    assert request(port, "POST", "/api/run", {}, {"Origin": "null"})[0] == 403
    assert request(port, "POST", "/api/run", {}, {"Host": f"example.com:{port}"})[0] == 403
    assert request(port, "POST", "/api/run", {}, {"Sec-Fetch-Site": "cross-site"})[0] == 403
    assert request(port, "POST", "/api/run", {}, {"Content-Type": "text/plain"})[0] == 415
    assert request(port, "POST", "/api/run", {"config": "x" * 65536})[0] == 413
    assert request(port, "POST", "/api/run", [], {})[0] == 400
    assert request(port, "POST", "/api/run", {}, {"Origin": f"http://127.0.0.1:{port}"})[0] == 200


def test_config_non_objects_are_rejected():
    with pytest.raises(ValueError, match="JSON object"):
        validated_config([])
