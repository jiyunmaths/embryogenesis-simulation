"""Apical–basal orientation from free cortex, with polar cortical mechanics."""

import numpy as np


def exposure_cue(phi, dx):
    """Interface-weighted outward normal on cortex unoccupied by other cells."""
    h = phi * phi * (3 - 2 * phi)
    occupied = h.sum(axis=0)
    cues = []
    for i, field in enumerate(phi):
        gradient = np.stack(np.gradient(field, dx))
        length = np.sqrt(np.sum(gradient**2, axis=0))
        normal = -gradient / np.maximum(length, 1e-10)
        shell = field * (1 - field)
        free = 1 - np.clip(2 * (occupied - h[i]), 0, 1)
        cue = np.sum(normal * (shell * free)[None], axis=(1, 2, 3)) / max(float(shell.sum()), 1e-12)
        # Avoid amplifying float32 roundoff on exactly symmetric free cells.
        if np.linalg.norm(cue) < 1e-6:
            cue[:] = 0
        cues.append(cue)
    return np.asarray(cues)


def evolve(polarity, cue, activator, delta, dt, rate=1., alignment=.25, decay=.5):
    drive = 2 * activator / (1 + activator)
    norm2 = np.sum(polarity**2, axis=1)
    drift = rate * drive[:, None] * cue + alignment * (delta @ polarity)
    drift -= (decay + norm2)[:, None] * polarity
    result = polarity + dt * drift
    result /= np.maximum(1, np.linalg.norm(result, axis=1))[:, None]
    return result


def tension_field(relative, p, base, contrast, width):
    direction = relative / np.sqrt(np.sum(relative**2, axis=0) + width**2)[None]
    return base * (1 - contrast * np.einsum("d,dijk->ijk", p, direction))


def flux_divergence(field, coefficient, dx):
    """Conservative face-flux div(gamma grad(phi)), with no boundary flux.

    Includes the spatial coefficient gradient; gamma*laplace(phi) alone is wrong.
    """
    result = np.zeros_like(field, dtype=float)
    for axis in range(3):
        left, right = [slice(None)] * 3, [slice(None)] * 3
        left[axis], right[axis] = slice(None, -1), slice(1, None)
        left, right = tuple(left), tuple(right)
        flux = .5 * (coefficient[left] + coefficient[right]) * (field[right] - field[left]) / dx**2
        result[left] += flux
        result[right] -= flux
    return result
