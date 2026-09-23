"""Finite-graph spectral preflight, without a 3D simulation."""

import argparse
import json
from pathlib import Path
import numpy as np

from .signaling import normalized_graph, stability, mode_growth, mode_transfer, integrate, gm_jacobian


def cycle(count):
    w = np.zeros((count, count))
    for i in range(count):
        w[i, (i + 1) % count] = w[(i + 1) % count, i] = 1
    return normalized_graph(w)


def verify_mode(graph, beta=2., da=1., dh=20.):
    rates = mode_growth(graph.eigenvalues, beta, da, dh)
    eligible = np.flatnonzero(graph.eigenvalues > 1e-10)
    if not len(eligible):
        return None
    mode = eligible[np.argmax(rates[eligible])]
    block = gm_jacobian(beta) - graph.eigenvalues[mode] * np.diag([da, dh])
    values, vectors = np.linalg.eig(block)
    index = np.argmax(values.real)
    if abs(values[index].imag) > 1e-10:
        return {"mode": int(mode), "verification": "oscillatory mode; no scalar exponential fit"}
    chemical = vectors[:, index].real
    scale = np.sqrt(np.where(graph.degree > 0, graph.degree, 1))
    spatial = graph.eigenvectors[:, mode] / scale
    a, h = 1 + 1e-6 * spatial * chemical[0], 1 + 1e-6 * spatial * chemical[1]
    ts, amplitude = [], []
    for k in range(501):
        ts.append(k * .01)
        amplitude.append(float(graph.eigenvectors[:, mode] @ (scale * (a - 1))))
        if k < 500:
            a, h = integrate(a, h, graph, .01, beta, da, dh)
    fitted = float(np.polyfit(ts, np.log(np.abs(amplitude)), 1)[0])
    return {"mode": int(mode), "predicted": float(rates[mode]), "fitted": fitted,
            "absolute_error": abs(fitted - rates[mode])}


def refinement_experiment(beta=2., da=1., dh=20.):
    """Illustrative 4→8→16 cycle refinement with inherited, unreset signals."""
    rng = np.random.default_rng(7)
    a = 1 + .001 * rng.normal(size=4)
    h = np.ones(4)
    graph = cycle(4)
    frames = []
    for step in range(1201):
        if step in (200, 700):  # Time 10 and 35; no new signal perturbation.
            a, h = np.repeat(a, 2), np.repeat(h, 2)
            graph = cycle(len(a))
        if step % 20 == 0:
            frames.append({"time": step * .05, "cells": len(a), "activator": a.tolist(),
                           "inhibitor": h.tolist(), "stability": stability(graph, beta, da, dh)})
        if step < 1200:
            a, h = integrate(a, h, graph, .05, beta, da, dh)
    return frames


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="New output directory")
    parser.add_argument("--checkpoint", type=Path, help="Read geometry only, including historical checkpoints")
    parser.add_argument("--contacts", type=Path, help="User-supplied symmetric contact matrix as .npy")
    parser.add_argument("--beta", type=float, default=2.)
    parser.add_argument("--da", type=float, default=1.)
    parser.add_argument("--dh", type=float, default=20.)
    parser.add_argument("--cutoff", type=float, default=.02)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("output directory already exists")
    graphs = {}
    for n in (4, 8, 16):
        graphs[f"cycle_{n}"] = cycle(n)
        graphs[f"complete_{n}"] = normalized_graph(np.ones((n, n)) - np.eye(n))
    if args.checkpoint:
        with np.load(args.checkpoint, allow_pickle=False) as data:
            field = data["phi"]
            config = json.loads(str(data["metadata"]))["config"]
        shell = (field * (1 - field)).reshape(len(field), -1)
        contacts = (shell @ shell.T) * (2 * config["extent"] / config["grid"])**3
        graphs["checkpoint_contacts"] = normalized_graph(contacts, args.cutoff)
    if args.contacts:
        graphs["supplied_contacts"] = normalized_graph(np.load(args.contacts, allow_pickle=False), args.cutoff)
    reports = {}
    for name, graph in graphs.items():
        report = stability(graph, args.beta, args.da, args.dh)
        report["linear_mode_verification"] = verify_mode(graph, args.beta, args.da, args.dh)
        reports[name] = report
        print(f'{name}: gap={report["spectral_gap"]}, unstable={report["unstable_modes"]}, '
              f'verification={report["linear_mode_verification"]}', flush=True)
    transfers = {}
    for n in (4, 8):
        # An explicit illustrative refinement, not a claim about real cell contacts.
        p = np.repeat(np.eye(n), 2, axis=0)
        transfers[f"cycle_{n}_to_{2*n}"] = mode_transfer(cycle(n), cycle(2*n), p)
    args.output.mkdir(parents=True)
    refinement = refinement_experiment(args.beta, args.da, args.dh)
    (args.output / "refinement_dynamics.json").write_text(json.dumps(refinement, indent=2) + "\n")
    (args.output / "analysis.json").write_text(json.dumps({
        "parameters": {"beta": args.beta, "da": args.da, "dh": args.dh, "cutoff": args.cutoff},
        "graphs": reports, "refinement_transfers": transfers,
        "interpretation": "Frozen-graph linear stability; normalized eigenvalues are not physical wavenumbers."
    }, indent=2) + "\n")
    for name, graph in graphs.items():
        np.savez_compressed(args.output / f"{name}.npz", contacts=graph.weights,
                            delta=graph.delta, eigenvalues=graph.eigenvalues, eigenvectors=graph.eigenvectors)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(9, 5), constrained_layout=True)
    values = np.linspace(0, 2, 500)
    ax.plot(values, mode_growth(values, args.beta, args.da, args.dh), color="black", label="2×2 modal growth envelope")
    for name, graph in graphs.items():
        ax.scatter(graph.eigenvalues, mode_growth(graph.eigenvalues, args.beta, args.da, args.dh),
                   s=42, alpha=.75, label=name)
    ax.axhline(0, color="gray", linewidth=1)
    ax.set(xlabel="Discrete normalized Laplacian eigenvalue λ", ylabel="Largest real growth rate",
           title="Only graph-supported eigenvalues can become unstable")
    ax.legend(fontsize=8, ncol=2)
    fig.savefig(args.output / "spectrum.png", dpi=150)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    for ax, (name, item) in zip(axes, transfers.items()):
        picture = ax.imshow(item["energy_fraction_new_by_old"], origin="lower", vmin=0, vmax=1, aspect="auto")
        ax.set_xticks(range(len(item["old_eigenvalues"])), [f'{v:.3f}' for v in item["old_eigenvalues"]])
        ax.set_yticks(range(len(item["new_eigenvalues"])), [f'{v:.3f}' for v in item["new_eigenvalues"]])
        ax.set(xlabel="Old eigenvalue group", ylabel="New eigenvalue group", title=name)
    fig.colorbar(picture, ax=axes, label="Inherited modal energy fraction")
    fig.savefig(args.output / "mode_transfer.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
