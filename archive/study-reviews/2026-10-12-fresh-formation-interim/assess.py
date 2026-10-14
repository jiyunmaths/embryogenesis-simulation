"""Read-only interim calculations; imports no simulation modules and uses no GPU."""
import datetime as dt
import hashlib
import json
from pathlib import Path
from zoneinfo import ZoneInfo

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
STUDY = ROOT / "output/studies/fresh-zygote-formation"
REVIEW = Path(__file__).resolve().parent


def capture(path):
    raw = path.read_bytes()
    evidence[str(path.relative_to(ROOT))] = hashlib.sha256(raw).hexdigest()
    return json.loads(raw)


def formation(times, contrast, threshold, duration, window):
    above = contrast > threshold
    first = float(times[np.flatnonzero(above)[0]]) if above.any() else None
    sustained = None
    for index in np.flatnonzero(above):
        end = times[index] + window
        if times[-1] < end - 1e-9:
            break
        if np.all(above[(times >= times[index] - 1e-9) & (times <= end + 1e-9)]):
            sustained = float(times[index])
            break
    complete = bool(times[-1] >= duration - 1e-9)
    late = times >= duration - window - 1e-9
    transitions = []
    previous = False
    for time, state in zip(times, above):
        if bool(state) != previous:
            transitions.append({"elapsed": float(time), "direction": "up" if state else "down"})
        previous = bool(state)
    return {
        "threshold": threshold,
        "first_crossing_elapsed": first,
        "first_sustained_onset_elapsed": sustained,
        "sustained_formation": sustained is not None,
        "late_persistent": bool(np.all(above[late])) if complete else None,
        "late_min_log_sd": float(np.min(contrast[late])) if complete else None,
        "final_log_sd": float(contrast[-1]),
        "sampled_transitions": transitions,
    }


evidence = {}
observed = dt.datetime.now(dt.timezone.utc)
protocol = capture(STUDY / "protocol.json")
status = capture(STUDY / "status.json")
protocol_hash = evidence[str((STUDY / "protocol.json").relative_to(ROOT))]
assert protocol_hash == status["protocol_sha256"]
source_failures = [name for name, digest in protocol["source_sha256"].items()
                   if hashlib.sha256(Path(name).read_bytes()).hexdigest() != digest]
assert not source_failures, source_failures
launch = capture(ROOT / "archive/study-launches/2026-10-12-fresh-zygote-formation/process-launch.json")
completed = []
prefix_seconds = []
active = None
recompute_errors = {"log_contrast_abs": 0.0, "spectrum_relative": 0.0,
                    "growth_abs": 0.0, "conductance_symmetry_abs": 0.0}
for job in protocol["jobs"]:
    directory = STUDY / job["key"]
    prefix_path = directory / "actual-start-prefix.json"
    if prefix_path.exists():
        prefix = capture(prefix_path)
        assert prefix["passed"] and prefix["exact_gpu_restart"]
        prefix_seconds.append(prefix["seconds"])
    result_path = directory / "gpu-result.json"
    is_active = job["key"] == status["active"]
    if not result_path.exists() and not is_active:
        continue
    result = capture(result_path) if result_path.exists() else None
    progress = capture(directory / "gpu-progress.json") if is_active else None
    history = capture(directory / "gpu-history.json")
    times = np.array([row["elapsed"] for row in history])
    contrast = np.array([np.std(np.log(np.asarray(row["chemistry"])[0])) for row in history])
    assert np.all(np.isfinite(contrast))
    contrast_error = float(np.max(np.abs(contrast - [row["log_activator_sd"] for row in history])))
    recompute_errors["log_contrast_abs"] = max(recompute_errors["log_contrast_abs"], contrast_error)
    conductances, mode_counts, growth_values = [], [], []
    for row in history:
        volumes = np.asarray(row["volumes"])
        delta = np.asarray(row["delta"])
        conductance = volumes[:, None] * delta
        recompute_errors["conductance_symmetry_abs"] = max(recompute_errors["conductance_symmetry_abs"],
                                                            float(np.max(np.abs(conductance - conductance.T))))
        conductances.append(float(np.triu(conductance, 1).sum()))
        root_volume = np.sqrt(volumes)
        similar = -root_volume[:, None] * delta / root_volume[None, :]
        eigenvalues = np.linalg.eigvalsh((similar + similar.T) / 2)
        recorded = np.asarray(row["spectrum"])
        error = np.linalg.norm(eigenvalues - recorded) / max(np.linalg.norm(recorded), 1e-30)
        recompute_errors["spectrum_relative"] = max(recompute_errors["spectrum_relative"], float(error))
        lam = eigenvalues[eigenvalues > 1e-10]
        beta = job["config"]["signal_beta"]
        a = 1 - job["config"]["signal_da"] * lam
        d = -beta - job["config"]["signal_dh"] * lam
        growth = np.real(((a + d) + np.sqrt(((a - d)**2 - 8 * beta).astype(complex))) / 2)
        maximum = float(growth.max())
        count = int(np.sum(growth > 1e-10))
        assert count == row["frozen_uniform_unstable_modes"]
        recompute_errors["growth_abs"] = max(recompute_errors["growth_abs"], abs(maximum - row["frozen_uniform_max_growth"]))
        mode_counts.append(count)
        growth_values.append(maximum)
    outcomes = [formation(times, contrast, threshold, protocol["duration"], protocol["sustained_window"])
                for threshold in protocol["thresholds"]]
    if result:
        assert result["passed"] and result["protocol_sha256"] == protocol_hash
        assert result["exact_restart_steps"] == {"12": 8, "36": 8}
        assert result["lineage_cycle_rng_unchanged"]
        for independent, recorded in zip(outcomes, result["formation"]):
            for key in ["first_crossing_elapsed", "first_sustained_onset_elapsed", "sustained_formation",
                        "late_persistent", "sampled_transitions"]:
                assert independent[key] == recorded[key], (job["key"], key)
        assert len(history) == 401
    summary = {"case": job["key"], "history": job["seed"], "branch": job["branch"],
               "level": job["level"], "records": len(history), "last_elapsed": float(times[-1]),
               "last_physical_time": history[-1]["time"], "formation": outcomes,
               "final_axis_ratio": history[-1]["axis_ratio"],
               "conductance_change_fraction": conductances[-1] / conductances[0] - 1,
               "frozen_uniform_unstable_modes": {"start": mode_counts[0], "end": mode_counts[-1],
                                                 "min": min(mode_counts), "max": max(mode_counts)},
               "frozen_uniform_growth": {"start": growth_values[0], "end": growth_values[-1],
                                         "min": min(growth_values), "max": max(growth_values)},
               "boundary_max": max(row["boundary_occupancy"] for row in history),
               "passed": result["passed"] if result else None,
               "seconds": result["seconds"] if result else progress["seconds"],
               "audit": result["audit"] if result else progress["audit"]}
    if result:
        completed.append(summary)
    else:
        active = summary

pairs = []
pair_seconds = []
for path in sorted(STUDY.glob("time-h*.json")):
    pair = capture(path)
    assert pair["passed"] and not pair["failed_gates"]
    pairs.append({key: pair[key] for key in ["seed", "branch", "passed", "homologous_errors",
                  "field_max", "dense_peak", "growth_abs_max", "topology_mismatch_times", "formation_checks"]})
    last_result = max((STUDY / case / "gpu-result.json").stat().st_mtime for case in pair["cases"])
    pair_seconds.append(path.stat().st_mtime - last_result)

means = {level["level"]: float(np.mean([case["seconds"] for case in completed
                                      if case["level"] == level["level"]])) for level in protocol["levels"]}
pending = [job for job in protocol["jobs"] if job["key"] not in [case["case"] for case in completed]
           and job["key"] != status["active"]]
pending_by_clock = {level: sum(job["level"] == level for job in pending) for level in means}
active_remaining = max(0.0, means[active["level"]] - active["seconds"]) if active else 0.0
prefix_mean = float(np.mean(prefix_seconds))
pair_mean = float(np.mean(pair_seconds))
remaining = (active_remaining + sum(pending_by_clock[level] * mean for level, mean in means.items())
             + len(pending) * prefix_mean + (len(protocol["time_pairs"]) - len(pairs)) * pair_mean
             + 120.0)
start = dt.datetime.fromisoformat(launch["utc"])
finish = observed + dt.timedelta(seconds=remaining)
report = {"observed_utc": observed.isoformat(), "observed_local": observed.astimezone(ZoneInfo("America/New_York")).isoformat(),
          "study": str(STUDY), "protocol_sha256": protocol_hash, "status_snapshot": status,
          "scope": "Read-only interim recomputation from recorded chemistry, volumes and transport matrices; no full checkpoint/dense-graph audit, GPU use, solver imports, or scientific acceptance.",
          "source_pins_checked": len(protocol["source_sha256"]), "source_pin_failures": source_failures,
          "independent_histories": 3, "new_independent_histories": 0,
          "completed": completed, "active": active, "completed_pairs": pairs,
          "independent_recomputation_errors": recompute_errors,
          "eta": {"production_mean_seconds_by_clock": means, "prefix_mean_seconds": prefix_mean,
                  "pair_post_result_mean_seconds": pair_mean, "pending_by_clock_excluding_active": pending_by_clock,
                  "active_remaining_seconds": active_remaining, "final_integrity_allowance_seconds": 120,
                  "remaining_hours": remaining / 3600, "total_estimated_hours": ((observed - start).total_seconds() + remaining) / 3600,
                  "central_finish_utc": finish.isoformat(), "central_finish_local": finish.astimezone(ZoneInfo("America/New_York")).isoformat(),
                  "planning_range_remaining_hours": [6, 7],
                  "condition": "All remaining registered gates pass, comparable GPU availability and I/O, and no stop; final independent review is additional time."},
          "acceptance": {"numerical_accepted": False, "scientific_accepted": False,
                         "reason": "Cohort incomplete; final independent saved-array audit required."},
          "evidence_sha256": evidence}
filename = observed.strftime("%Y-%m-%d-fresh-formation-interim-%H%M%S")
destination = ROOT / "docs/lab_notebook" / (filename + ".json")
destination.write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps({"report": str(destination.relative_to(ROOT)), "observed_local": report["observed_local"],
                  "completed": len(completed), "pairs": len(pairs), "active": active,
                  "eta": report["eta"], "recomputation": recompute_errors,
                  "history1101_fine": [case for case in completed if case["history"] == 1101 and case["level"] == "fine"]}, indent=2))
