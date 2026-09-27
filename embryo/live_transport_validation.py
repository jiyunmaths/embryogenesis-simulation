"""Calibrated-contact refinement on a fixed 3D slab domain (unit cross-section).

The slabs are manufactured numerical compartments, not biological cleavage.
This isolates the overlap-area closure and spectral scaling from cell mechanics.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from .transport import contact_transport, transport_graph
from .signaling import stability, normalized_graph
from .graph_analysis import verify_mode


def slab_graph(count, samples_per_width=6):
    width = 1 / (16 * count)
    points = count * 16 * samples_per_width
    dx = 1 / points
    x = (np.arange(points) + .5) * dx
    boundaries = np.arange(1, count) / count
    transitions = .5 * (1 + np.tanh((x[None] - boundaries[:, None]) / (np.sqrt(2) * width)))
    phi = np.concatenate([np.ones((1, points)), transitions], axis=0) - np.concatenate(
        [transitions, np.zeros((1, points))], axis=0)
    occupancy = phi**2 * (3 - 2 * phi)
    volumes = occupancy.sum(axis=1) * dx
    centers = np.zeros((count, 3))
    centers[:, 0] = occupancy @ x * dx / volumes
    shell = phi * (1 - phi)
    contacts = shell @ shell.T * dx  # unit transverse area
    np.fill_diagonal(contacts, 0)
    return transport_graph(contact_transport(contacts, volumes, centers, width, .02)), contacts


def run():
    rows = []
    for count in (8, 16, 32, 64):
        graph, contacts = slab_graph(count)
        report = stability(graph, da=.02, dh=.4)
        rows.append({"compartments": count,
                     "first_eigenvalue": float(graph.eigenvalues[1]),
                     "relative_first_mode_error": float(abs(graph.eigenvalues[1] / np.pi**2 - 1)),
                     "unstable_modes": report["unstable_modes"],
                     "physical_lambda_band": report["continuous_lambda_band"],
                     "mass_residual": float(np.max(abs(graph.masses @ graph.delta))),
                     "constant_residual": float(np.max(abs(graph.delta @ np.ones(count)))),
                     "random_walk_first_eigenvalue": float(normalized_graph(contacts).eigenvalues[1])})
    graph, _ = slab_graph(16)
    verification = verify_mode(graph, da=.02, dh=.4)
    errors = np.array([row["relative_first_mode_error"] for row in rows])
    orders = np.log2(errors[:-1] / errors[1:])
    checks = {"first_mode_error_decreases": bool(np.all(np.diff(errors) < 0)),
              "last_order_above_1_8": bool(orders[-1] > 1.8),
              "finest_first_mode_error_below_0_001": bool(errors[-1] < .001),
              "resolved_unstable_modes_match_continuum": all(row["unstable_modes"] == [1] for row in rows[1:]),
              "linear_growth_matches_prediction": bool(verification["absolute_error"] < 2e-5),
              "mass_and_constants_preserved": all(max(row["mass_residual"], row["constant_residual"]) < 1e-9 for row in rows)}
    return {"description": __doc__, "coefficients": {"beta": 2., "da": .02, "dh": .4},
            "continuum_first_eigenvalue": float(np.pi**2), "continuum_unstable_modes": [1],
            "refinement": rows, "orders": orders.tolist(), "linear_mode_verification": verification,
            "checks": checks, "passed": all(checks.values()),
            "limitations": "Flat complementary interfaces and orthogonal slab centers only. Does not validate arbitrary deformed contacts or the full coupled embryo continuum limit."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    report = run()
    (args.output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))
    if not report['passed']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
