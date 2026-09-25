"""Conservative concentration transport between fixed 3D compartments.

Conductances contain geometry (interface area / center separation), not a
diffusion coefficient. With compartment volumes M = diag(V), the positive
stiffness matrix K = diag(G 1) - G gives Delta = -M^-1 K. This is a physical
finite-volume discretization, separate from the normalized activity exchange
used by the existing deformable-cell simulation.
"""

from dataclasses import dataclass
import math

import numpy as np
from scipy import sparse

from .signaling import gm_reaction


@dataclass
class Transport:
    volumes: np.ndarray
    conductance: sparse.csr_matrix
    delta: sparse.csr_matrix
    symmetric: sparse.csr_matrix
    centers: np.ndarray | None = None
    shape: tuple[int, int, int] | None = None


def conservative_transport(conductance, volumes):
    """Construct conservative exchange from symmetric, nonnegative G.

    G_ij has dimensions length^(d-2) for geometric bulk diffusion in d
    dimensions; volumes have dimensions length^d. Multiplication by a physical
    diffusivity then gives concentration change per time. Diagonal conductances
    must be zero. Isolated compartments are supported. Tiny relative asymmetry
    (at most 1e-12 of the largest conductance) is averaged to enforce exact
    reciprocal fluxes, while larger asymmetry is rejected.

    No dense N-by-N matrices or eigenvectors are created for sparse input.
    ``symmetric`` is M^-1/2 K M^-1/2, a positive-semidefinite matrix with the
    same eigenvalues as -Delta; physical right modes are q / sqrt(V).
    """
    if np.iscomplexobj(conductance) or np.iscomplexobj(volumes):
        raise ValueError("conductances and volumes must be real")
    if not sparse.issparse(conductance) and np.asarray(conductance).ndim != 2:
        raise ValueError("conductance must be a nonempty square matrix")
    try:
        matrix = sparse.csr_matrix(conductance, dtype=float, copy=True)
        volume = np.asarray(volumes, dtype=float).copy()
    except (TypeError, ValueError) as error:
        raise ValueError("conductances must form a matrix and volumes a vector") from error
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1] or matrix.shape[0] == 0:
        raise ValueError("conductance must be a nonempty square matrix")
    count = matrix.shape[0]
    if volume.shape != (count,) or not np.all(np.isfinite(volume)) or np.any(volume <= 0):
        raise ValueError("volumes must be a matching vector of positive finite values")
    matrix.sum_duplicates()
    matrix.eliminate_zeros()
    if not np.all(np.isfinite(matrix.data)) or np.any(matrix.data < 0):
        raise ValueError("conductances must be finite and nonnegative")
    if np.any(matrix.diagonal() != 0):
        raise ValueError("conductance must have a zero diagonal")
    difference = matrix - matrix.T
    maximum = float(np.max(matrix.data)) if matrix.nnz else 0.0
    if difference.nnz and np.max(np.abs(difference.data)) > 1e-12 * maximum:
        raise ValueError("conductance must be symmetric")
    matrix = (matrix * .5 + matrix.T * .5).tocsr()
    degree = np.asarray(matrix.sum(axis=1)).ravel()
    stiffness = sparse.diags(degree, format="csr") - matrix
    with np.errstate(over="raise", divide="raise", invalid="raise"):
        try:
            inverse = 1.0 / volume
            inverse_sqrt = 1.0 / np.sqrt(volume)
            delta = (-sparse.diags(inverse) @ stiffness).tocsr()
            symmetric = (sparse.diags(inverse_sqrt) @ stiffness @ sparse.diags(inverse_sqrt)).tocsr()
        except FloatingPointError as error:
            raise ValueError("volumes and conductances produce nonfinite transport rates") from error
    if not np.all(np.isfinite(delta.data)) or not np.all(np.isfinite(symmetric.data)):
        raise ValueError("volumes and conductances produce nonfinite transport rates")
    return Transport(volume, matrix, delta, symmetric)


def cartesian_transport(n=None, length=1.0, edges=None):
    """Finite-volume diffusion on an orthogonal Cartesian box, no exterior flux.

    For a uniform cube [0, length]^3, ``n`` is a positive integer or a tuple
    (nx, ny, nz). Alternatively, provide ``edges=(x_edges,y_edges,z_edges)`` and
    leave n=None; strictly increasing axis edges may be nonuniform and determine
    the box dimensions. ``length`` is used only for the uniform-cube option.

    Flatten arrays in C order after meshgrid(indexing="ij"): the z index varies
    fastest. Volumes are products of the local axis widths. Each interior-face
    conductance is its true area divided by the distance between the adjacent
    compartment centers. No conductance is added through exterior faces.
    """
    if edges is None:
        if isinstance(length, (bool, np.bool_)) or not np.isscalar(length):
            raise ValueError("length must be a positive finite number")
        try:
            length = float(length)
        except (ValueError, TypeError, OverflowError) as error:
            raise ValueError("length must be a positive finite number") from error
        if not np.isfinite(length) or length <= 0:
            raise ValueError("length must be a positive finite number")
        counts = (n, n, n) if isinstance(n, (int, np.integer)) else n
        if not isinstance(counts, (tuple, list)) or len(counts) != 3 or any(
            isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value < 1
            for value in counts
        ):
            raise ValueError("n must be a positive integer or three positive integers")
        axis_edges = tuple(np.linspace(0.0, length, int(value) + 1) for value in counts)
    else:
        if n is not None:
            raise ValueError("provide n or explicit edges, not both")
        if not isinstance(edges, (tuple, list)) or len(edges) != 3:
            raise ValueError("edges must contain three one-dimensional axis edge arrays")
        if any(np.iscomplexobj(axis) for axis in edges):
            raise ValueError("axis edges must be real")
        axis_edges = tuple(np.asarray(axis, dtype=float) for axis in edges)
        if any(axis.ndim != 1 or len(axis) < 2 or not np.all(np.isfinite(axis)) or
               np.any(np.diff(axis) <= 0) for axis in axis_edges):
            raise ValueError("axis edges must be finite, strictly increasing vectors of length >= 2")
    widths = tuple(np.diff(axis) for axis in axis_edges)
    axes = tuple(axis[:-1] + .5 * width for axis, width in zip(axis_edges, widths))
    shape = tuple(len(axis) for axis in axes)
    volume = (widths[0][:, None, None] * widths[1][None, :, None] * widths[2][None, None, :]).ravel()
    centers = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(-1, 3)
    indices = np.arange(math.prod(shape)).reshape(shape)
    rows, columns, values = [], [], []
    for axis in range(3):
        if shape[axis] == 1:
            continue
        left, right = [slice(None)] * 3, [slice(None)] * 3
        left[axis], right[axis] = slice(None, -1), slice(1, None)
        face_shape = list(shape)
        face_shape[axis] -= 1
        area = 1.0
        for transverse in range(3):
            if transverse != axis:
                area_shape = [1, 1, 1]
                area_shape[transverse] = shape[transverse]
                area = area * widths[transverse].reshape(area_shape)
        distance_shape = [1, 1, 1]
        distance_shape[axis] = shape[axis] - 1
        # Adjacent center separation, without subtracting large coordinates.
        distance = (.5 * (widths[axis][:-1] + widths[axis][1:])).reshape(distance_shape)
        conductances = np.broadcast_to(area / distance, face_shape).ravel()
        low, high = indices[tuple(left)].ravel(), indices[tuple(right)].ravel()
        rows.extend((low, high))
        columns.extend((high, low))
        values.extend((conductances, conductances))
    if rows:
        matrix = sparse.coo_matrix((np.concatenate(values), (np.concatenate(rows), np.concatenate(columns))),
                                   shape=(len(volume), len(volume))).tocsr()
    else:
        matrix = sparse.csr_matrix((1, 1), dtype=float)
    transport = conservative_transport(matrix, volume)
    transport.centers, transport.shape = centers, shape
    return transport


def integrate_gm(a, h, transport, dt, beta=2.0, da=.02, dh=.4):
    """Advance GM concentrations with conservative transport and SSP-RK2.

    Production terms are nonnegative. The maximum loss coefficient is
    max(1 + da*max_exit_rate, beta + dh*max_exit_rate), where exit rates are
    -diag(Delta). Internal substeps keep dt_sub*loss <= .2, making the forward
    Euler stages and their SSP convex combination positivity-preserving. This
    restriction grows with inverse compartment spacing squared; rates are not
    normalized by graph degree. Zero diffusivities are allowed as controls.

    Reaction changes molecular amounts; only the exchange contribution conserves
    volume-weighted amount. Input arrays are copied once, not changed in place.
    """
    if any(isinstance(value, (bool, np.bool_)) or not np.isscalar(value)
           for value in (dt, beta, da, dh)):
        raise ValueError("dt, beta and diffusion coefficients must be real numbers")
    try:
        dt, beta, da, dh = map(float, (dt, beta, da, dh))
    except (ValueError, TypeError, OverflowError) as error:
        raise ValueError("dt, beta and diffusion coefficients must be real numbers") from error
    if not np.isfinite([dt, beta, da, dh]).all() or dt <= 0 or beta <= 0 or min(da, dh) < 0:
        raise ValueError("dt and beta must be positive finite; diffusion coefficients must be nonnegative finite")
    if np.iscomplexobj(a) or np.iscomplexobj(h):
        raise ValueError("signals must be real")
    a, h = np.asarray(a, dtype=float).copy(), np.asarray(h, dtype=float).copy()
    if (a.shape != transport.volumes.shape or h.shape != a.shape or
        not np.isfinite(a).all() or not np.isfinite(h).all() or np.any(a < 0) or np.any(h <= 0)):
        raise ValueError("signals must match transport size, be finite, and satisfy a>=0 and h>0")
    exit_rate = float(np.max(-transport.delta.diagonal()))
    loss = max(1.0 + da * exit_rate, beta + dh * exit_rate)
    required = dt * loss / .2
    if not np.isfinite(required):
        raise ValueError("dt and transport rates require a nonfinite substep count")
    substeps = max(1, math.ceil(required))
    step = dt / substeps

    def rhs(x, y):
        fa, fh = gm_reaction(x, y, beta)
        return fa + da * (transport.delta @ x), fh + dh * (transport.delta @ y)

    with np.errstate(over="raise", divide="raise", invalid="raise"):
        for _ in range(substeps):
            fa, fh = rhs(a, h)
            first_a, first_h = a + step * fa, h + step * fh
            fa, fh = rhs(first_a, first_h)
            a *= .5
            h *= .5
            a += .5 * (first_a + step * fa)
            h += .5 * (first_h + step * fh)
            if np.any(a < 0) or np.any(h <= 0):
                raise FloatingPointError("GM concentrations lost positivity; reduce dt/check kinetics")
    if not np.isfinite(a).all() or not np.isfinite(h).all():
        raise FloatingPointError("GM concentrations became nonfinite; reduce dt/check kinetics")
    return a, h
