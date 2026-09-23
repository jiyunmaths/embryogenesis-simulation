"""Gierer–Meinhardt kinetics and exact finite normalized-graph stability."""

from dataclasses import dataclass
import numpy as np


@dataclass
class Graph:
    weights: np.ndarray
    degree: np.ndarray
    delta: np.ndarray
    symmetric: np.ndarray
    eigenvalues: np.ndarray
    eigenvectors: np.ndarray


def normalized_graph(contacts, relative_cutoff=0.02):
    """Random-walk Delta, similar to -L_sym; isolated vertices have Delta=0.

    Signals are nonconserved per-cell activities. This operator preserves constants
    and degree-weighted sums on a fixed graph, not volume-weighted molecular mass.
    """
    w = np.asarray(contacts, dtype=float).copy()
    if w.ndim != 2 or w.shape[0] != w.shape[1] or not len(w):
        raise ValueError("contacts must be a nonempty square matrix")
    if not np.all(np.isfinite(w)) or np.any(w < 0) or not np.allclose(w, w.T):
        raise ValueError("contacts must be finite, nonnegative and symmetric")
    if not 0 <= relative_cutoff < 1:
        raise ValueError("relative_cutoff must be in [0,1)")
    np.fill_diagonal(w, 0)
    w[w < max(1e-12, relative_cutoff * w.max())] = 0
    degree = w.sum(axis=1)
    active = degree > 0
    inverse = np.divide(1.0, degree, out=np.zeros_like(degree), where=active)
    delta = inverse[:, None] * w - np.diag(active.astype(float))
    symmetric = np.diag(active.astype(float)) - np.sqrt(inverse[:, None] * inverse[None, :]) * w
    eigenvalues, eigenvectors = np.linalg.eigh(symmetric)
    eigenvalues[np.abs(eigenvalues) < 1e-12] = 0
    return Graph(w, degree, delta, symmetric, eigenvalues, eigenvectors)


def gm_reaction(a, h, beta=2.0):
    """Positive homogeneous equilibrium (a,h)=(1,1); beta>1 gives local stability."""
    return a * a / h - a, beta * (a * a - h)


def gm_jacobian(beta=2.0):
    return np.array([[1.0, -1.0], [2 * beta, -beta]])


def mode_growth(eigenvalues, beta=2.0, da=1.0, dh=20.0):
    j = gm_jacobian(beta)
    return np.array([np.max(np.linalg.eigvals(j - value * np.diag([da, dh])).real)
                     for value in np.atleast_1d(eigenvalues)])


def stability(graph, beta=2.0, da=1.0, dh=20.0):
    if min(beta, da, dh) <= 0 or not np.isfinite([beta, da, dh]).all():
        raise ValueError("beta, da and dh must be positive and finite")
    growth = mode_growth(graph.eigenvalues, beta, da, dh)
    local = float(mode_growth([0], beta, da, dh)[0])
    nonzero = graph.eigenvalues > 1e-10
    # det(J-lambda D)=beta+(beta*Da-Dh)*lambda+Da*Dh*lambda^2.
    roots = np.roots([da * dh, beta * da - dh, beta])
    band = sorted(float(x.real) for x in roots if abs(x.imag) < 1e-10 and x.real > 0)
    band = band if len(band) == 2 and local < 0 else None
    unstable = nonzero & (growth > 1e-10)
    return {
        "cells": len(graph.degree), "edges": int(np.count_nonzero(np.triu(graph.weights, 1))),
        "components": int(np.count_nonzero(~nonzero)),
        "smallest_positive_eigenvalue": float(graph.eigenvalues[nonzero].min()) if np.any(nonzero) else None,
        "spectral_gap": float(graph.eigenvalues[1]) if len(graph.degree) > 1 else None,
        "eigenvalues": graph.eigenvalues.tolist(), "growth_rates": growth.tolist(),
        "local_growth_rate": local, "local_stable": local < -1e-10,
        "continuous_lambda_band": band,
        "unstable_modes": np.flatnonzero(unstable).tolist(),
        "diffusion_driven_instability": bool(local < -1e-10 and np.any(unstable)),
        "maximum_spatial_growth": float(growth[nonzero].max()) if np.any(nonzero) else None,
    }


def integrate(a, h, graph, dt, beta=2.0, da=1.0, dh=20.0):
    """SSP-RK2 with positive Euler stages for this production/loss system."""
    if dt <= 0 or min(beta, da, dh) <= 0 or not np.isfinite([dt, beta, da, dh]).all():
        raise ValueError("dt and signaling parameters must be positive and finite")
    a, h = np.asarray(a, dtype=float).copy(), np.asarray(h, dtype=float).copy()
    if a.shape != graph.degree.shape or h.shape != a.shape or np.any(a < 0) or np.any(h <= 0):
        raise ValueError("signals must match graph size, with a>=0 and h>0")
    substeps = max(1, int(np.ceil(dt * max(1 + da, beta + dh) / 0.2)))
    step = dt / substeps

    def rhs(x, y):
        f, g = gm_reaction(x, y, beta)
        return f + da * (graph.delta @ x), g + dh * (graph.delta @ y)

    for _ in range(substeps):
        fa, fh = rhs(a, h)
        first_a, first_h = a + step * fa, h + step * fh
        fa, fh = rhs(first_a, first_h)
        a, h = .5 * a + .5 * (first_a + step * fa), .5 * h + .5 * (first_h + step * fh)
        if not np.isfinite(a).all() or not np.isfinite(h).all() or np.any(a < 0) or np.any(h <= 0):
            raise FloatingPointError("GM signals lost positivity or became nonfinite; reduce dt/check kinetics")
    return a, h


def eigenvalue_groups(values, tolerance=1e-8):
    groups = []
    for index, value in enumerate(values):
        if not groups or abs(value - values[groups[-1][0]]) > tolerance:
            groups.append([index])
        else:
            groups[-1].append(index)
    return groups


def mode_transfer(before, after, prolongation):
    """Basis-invariant energy transfer between degenerate spectral subspaces.

    P copies parent activities to daughters. Transform random-walk right modes
    through P into the new degree-weighted symmetric eigenbasis. Summing squares
    over old and new degenerate groups avoids arbitrary eigenvector signs/rotations.
    """
    p = np.asarray(prolongation, dtype=float)
    if p.shape != (len(after.degree), len(before.degree)) or not np.allclose(p.sum(axis=1), 1):
        raise ValueError("prolongation must have new-by-old shape and preserve constants")
    old_scale = np.sqrt(np.where(before.degree > 0, before.degree, 1))
    new_scale = np.sqrt(np.where(after.degree > 0, after.degree, 1))
    transfer = after.eigenvectors.T @ (new_scale[:, None] * p) @ (before.eigenvectors / old_scale[:, None])
    old_groups, new_groups = eigenvalue_groups(before.eigenvalues), eigenvalue_groups(after.eigenvalues)
    energy = np.array([[np.sum(transfer[np.ix_(new, old)]**2) for old in old_groups] for new in new_groups])
    energy /= np.maximum(energy.sum(axis=0, keepdims=True), 1e-30)
    return {"old_eigenvalues": [float(before.eigenvalues[g[0]]) for g in old_groups],
            "new_eigenvalues": [float(after.eigenvalues[g[0]]) for g in new_groups],
            "old_multiplicities": [len(g) for g in old_groups],
            "new_multiplicities": [len(g) for g in new_groups],
            "energy_fraction_new_by_old": energy.tolist()}


def cleavage_prolongation(count, mother):
    return np.eye(count)[[i for i in range(count) if i != mother] + [mother, mother]]
