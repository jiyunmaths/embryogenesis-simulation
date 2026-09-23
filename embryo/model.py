"""Deformable cells, discrete activator–inhibitor signaling, and polar mechanics.

All quantities are dimensionless. See docs/model.md for assumptions and energy.
"""

from dataclasses import asdict, dataclass
import numpy as np
from scipy.ndimage import laplace
from .signaling import normalized_graph, stability, integrate, mode_transfer, cleavage_prolongation
from .polarity import exposure_cue, evolve, tension_field, flux_divergence


@dataclass
class Config:
    seed: int = 7
    grid: int = 40
    extent: float = 1.6
    dt: float = 0.015
    steps: int = 1000
    max_cells: int = 16
    division_interval: float = 2.0
    cycle_jitter: float = 0.12
    division_orientation: str = "shape"
    axis_degeneracy: float = 0.03
    cytokinesis_duration: float = 0.9
    ring_strength: float = 8.0
    neck_threshold: float = 0.15
    division_overlap_tolerance: float = 0.002
    interface_width: float = 0.085
    surface_tension: float = 1.0
    volume_stiffness: float = 12.0
    repulsion: float = 3.0
    adhesion: float = 4.0
    fate_adhesion: float = 0.35
    fate_tension: float = 0.25
    fate_rate: float = 0.8
    neighbor_inhibition: float = 0.0
    exposure_bias: float = 0.0
    fate_noise: float = 0.0
    partition_noise: float = 0.0
    signaling: bool = True
    signal_beta: float = 2.0
    signal_da: float = 1.0
    signal_dh: float = 20.0
    signal_partition_noise: float = 0.001
    signal_fate_gain: float = 1.0
    graph_contact_cutoff: float = 0.02
    polarity_enabled: bool = True
    polarity_rate: float = 1.0
    polarity_alignment: float = 0.25
    polarity_decay: float = 0.5
    polarity_tension: float = 0.35
    competence_cells: int = 4
    fate_threshold: float = 0.55
    save_every: int = 20
    feedback: bool = True
    differentiation: bool = True

    def validate(self):
        for name in ("grid", "steps", "max_cells", "competence_cells", "save_every"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 1:
                raise ValueError(f"{name} must be a positive integer")
        if self.grid < 12:
            raise ValueError("grid must be at least 12")
        for name in ("extent", "dt", "division_interval", "interface_width",
                     "surface_tension", "volume_stiffness", "fate_threshold",
                     "cytokinesis_duration", "ring_strength", "signal_beta", "signal_da", "signal_dh"):
            if not np.isfinite(getattr(self, name)) or getattr(self, name) <= 0:
                raise ValueError(f"{name} must be finite and positive")
        for name in ("repulsion", "adhesion", "fate_rate", "neighbor_inhibition",
                     "exposure_bias", "fate_noise", "partition_noise", "signal_partition_noise",
                     "signal_fate_gain", "polarity_rate", "polarity_alignment", "polarity_decay"):
            if not np.isfinite(getattr(self, name)) or getattr(self, name) < 0:
                raise ValueError(f"{name} must be finite and nonnegative")
        for name in ("cycle_jitter", "fate_adhesion", "fate_tension", "graph_contact_cutoff", "polarity_tension"):
            if not 0 <= getattr(self, name) < 1:
                raise ValueError(f"{name} must be in [0, 1)")
        if self.division_orientation not in ("shape", "isotropic"):
            raise ValueError("division_orientation must be 'shape' or 'isotropic'")
        for name in ("axis_degeneracy", "neck_threshold", "division_overlap_tolerance"):
            if not 0 < getattr(self, name) < 1:
                raise ValueError(f"{name} must be in (0, 1)")
        if self.cytokinesis_duration < 4 * self.dt:
            raise ValueError("cytokinesis_duration must span at least four time steps")
        if self.dt * self.ring_strength > 0.2:
            raise ValueError("reduce dt: dt * ring_strength must be <= 0.2")
        dx = 2 * self.extent / self.grid
        if self.interface_width < 0.6 * dx:
            raise ValueError("interface_width must be >= 0.6 grid spacings")
        # Conservative explicit diffusion restriction; reactions are checked by convergence tests.
        if self.dt * (self.polarity_rate * 2 + self.polarity_alignment * 2 + self.polarity_decay + 3) > .5:
            raise ValueError("reduce dt for polarity dynamics")
        if self.dt * self.surface_tension * (1 + self.fate_tension) * (1 + self.polarity_tension) * self.interface_width**2 / dx**2 > 1 / 6:
            raise ValueError("dt is too large for explicit 3D diffusion")


def occupancy(phi):
    return phi * phi * (3 - 2 * phi)


class Simulation:
    def __init__(self, config=None):
        self.config = config or Config()
        self.config.validate()
        c = self.config
        mechanical_seed, fate_seed, signal_seed = np.random.SeedSequence(c.seed).spawn(3)
        self.rng = np.random.default_rng(mechanical_seed)
        self.fate_rng = np.random.default_rng(fate_seed)
        self.signal_rng = np.random.default_rng(signal_seed)
        self.dx = 2 * c.extent / c.grid
        axis = (np.arange(c.grid) + 0.5) * self.dx - c.extent
        self.xyz = np.stack(np.meshgrid(axis, axis, axis, indexing="ij"))
        radius = np.sqrt(np.sum(self.xyz**2, axis=0))
        initial = 0.5 * (1 - np.tanh((radius - 0.8) / (np.sqrt(2) * c.interface_width)))
        self.phi = initial[None].astype(np.float32)
        self.target = self.volumes().copy()
        self.initial_volume = float(self.target.sum())
        self.fate = np.array([0.0])
        self.activator = np.ones(1)
        self.inhibitor = np.ones(1)
        self.polarity = np.zeros((1, 3))
        self.graph_events = []
        self.ids = np.array([0], dtype=int)
        self.parents = np.array([-1], dtype=int)
        self.next_id = 1
        self.time = 0.0
        self.step_number = 0
        self.due = np.array([self._cycle()])
        self.lineage = [{"id": 0, "parent": -1, "birth": 0.0, "division": None}]
        self.divisions = {}
        self.volume_projection_max = 0.0
        self.clipped_fraction = 0.0
        self.last_contact = np.zeros((1, 1))
        self.last_exposure = np.ones(1)

    def _cycle(self):
        c = self.config
        return c.division_interval * self.rng.uniform(1 - c.cycle_jitter, 1 + c.cycle_jitter)

    def volumes(self):
        return occupancy(self.phi).sum(axis=(1, 2, 3), dtype=np.float64) * self.dx**3

    def centers(self):
        h = occupancy(self.phi)
        return np.einsum("nijk,dijk->nd", h, self.xyz) * self.dx**3 / self.volumes()[:, None]

    def division_direction(self, index):
        """Spindle axis; the cleavage plane is perpendicular to this direction."""
        sample = self.rng.normal(size=3)
        if self.config.division_orientation == "isotropic":
            return sample / np.linalg.norm(sample)
        h = occupancy(self.phi[index]).ravel()
        offset = (self.xyz - self.centers()[index, :, None, None, None]).reshape(3, -1)
        covariance = (offset * h) @ offset.T / h.sum()
        eigenvalues, eigenvectors = np.linalg.eigh(covariance)
        # Project onto the nearly maximal eigenspace. A sphere has no unique
        # long axis; an oblate cell has a plane of equally valid spindle axes.
        basis = eigenvectors[:, eigenvalues >= eigenvalues[-1] * (1 - self.config.axis_degeneracy)]
        axis = basis @ (basis.T @ sample)
        return axis / np.linalg.norm(axis)

    def _division_coordinates(self, index, event):
        relative = self.xyz - self.centers()[index, :, None, None, None]
        axial = np.einsum("d,dijk->ijk", np.asarray(event["axis"]), relative)
        radial = np.sqrt(np.maximum(np.sum(relative**2, axis=0) - axial**2, 0))
        return axial - event["offset"], radial

    def _project_volume(self, field, volume):
        """Interface-local scalar pressure correction for a dividing mother.

        This enforces incompressibility numerically, not a fluid pressure solve.
        The monotone bounded update leaves remote empty grid points empty.
        """
        weight = 6 * field * (1 - field)
        shift = 0.0
        for _ in range(20):
            candidate = np.clip(field + shift * weight, 0, 1)
            error = volume - occupancy(candidate).sum(dtype=np.float64) * self.dx**3
            if abs(error) <= 1e-7 * volume:
                self.volume_projection_max = max(self.volume_projection_max, abs(shift))
                return candidate
            slope = np.sum(6 * candidate * (1 - candidate) * weight, dtype=np.float64) * self.dx**3
            if slope <= 1e-12:
                break
            shift += error / slope
        raise FloatingPointError("cytokinesis volume projection failed; reduce dt or refine the grid")

    def contacts(self):
        shell = self.phi * (1 - self.phi)
        flat = shell.reshape(len(shell), -1)
        contacts = (flat @ flat.T).astype(float) * self.dx**3
        np.fill_diagonal(contacts, 0)
        others = occupancy(self.phi).sum(axis=0)[None] - occupancy(self.phi)
        # An interface-weighted exposure proxy, not a measured geometric surface fraction.
        blocked = np.clip(2 * others, 0, 1)
        exposure = 1 - (shell * blocked).sum(axis=(1, 2, 3)) / np.maximum(shell.sum(axis=(1, 2, 3)), 1e-12)
        return contacts, exposure

    def update_fate(self, contacts, exposure):
        c = self.config
        if not c.differentiation or len(self.fate) < c.competence_cells:
            return
        neighbor = contacts @ self.fate / np.maximum(contacts.sum(axis=1), 1e-12)
        # Local inhibition is a hypothesis to explore, not a species-specific pathway.
        drift = c.fate_rate * (self.fate - self.fate**3
                 - c.neighbor_inhibition * neighbor + c.exposure_bias * (exposure - 0.5)
                 + (c.signal_fate_gain * (self.activator - 1) if c.signaling else 0))
        self.fate += c.dt * drift + c.fate_noise * np.sqrt(c.dt) * self.fate_rng.normal(size=len(self.fate))

    def graph_snapshot(self, contacts=None):
        if contacts is None:
            contacts = self.contacts()[0]
        c = self.config
        graph = normalized_graph(contacts, c.graph_contact_cutoff)
        return {"time": self.time, "ids": self.ids.tolist(), "weights": graph.weights.tolist(),
                "activator": self.activator.tolist(), "inhibitor": self.inhibitor.tolist(),
                "polarity": self.polarity.tolist(),
                **stability(graph, c.signal_beta, c.signal_da, c.signal_dh)}

    def mechanical_step(self):
        c = self.config
        phi = self.phi
        shell2 = (phi * (1 - phi))**2
        derivative = 2 * phi * (1 - phi) * (1 - 2 * phi)
        fate = np.tanh(self.fate)
        adhesion = np.full((len(phi), len(phi)), c.adhesion)
        tension = np.full(len(phi), c.surface_tension)
        if c.feedback:
            adhesion *= 1 + c.fate_adhesion * fate[:, None] * fate[None, :]
            tension *= 1 + c.fate_tension * fate
        np.fill_diagonal(adhesion, 0)
        attract = (adhesion @ shell2.reshape(len(phi), -1)).reshape(phi.shape)
        exclude = np.sum(phi**2, axis=0)[None] - phi**2
        volume_force = c.volume_stiffness * (self.target - self.volumes()) / self.target
        updated = np.empty_like(phi)
        self.volume_projection_max = 0.0
        centers = self.centers()
        for i in range(len(phi)):
            diffusion = laplace(phi[i], mode="nearest") / self.dx**2
            force = tension[i] * (c.interface_width**2 * diffusion - derivative[i])
            if c.feedback and c.polarity_enabled and np.linalg.norm(self.polarity[i]) > 1e-10:
                gamma = tension_field(self.xyz - centers[i, :, None, None, None], self.polarity[i],
                                      tension[i], c.polarity_tension, c.interface_width)
                force = c.interface_width**2 * flux_divergence(phi[i], gamma, self.dx) - gamma * derivative[i]
            force += volume_force[i] * 6 * phi[i] * (1 - phi[i])
            force -= c.repulsion * phi[i] * exclude[i]
            force += derivative[i] * attract[i]
            event = self.divisions.get(int(self.ids[i]))
            if event is not None:
                axial, radial = self._division_coordinates(i, event)
                progress = np.clip((self.time - event["start"]) / c.cytokinesis_duration, 0, 1)
                ramp = progress**2 * (3 - 2 * progress)
                ring_radius = event["radius"] * (1 - ramp)
                band = np.exp(-0.5 * (axial / c.interface_width)**2)
                outside = 0.5 * (1 + np.tanh((radial - ring_radius) / c.interface_width))
                # A moving annular potential acts at the diffuse cell interface;
                # it progressively excludes cytoplasm from the advancing furrow.
                force -= c.ring_strength * ramp * band * outside * 6 * phi[i] * (1 - phi[i])
            updated[i] = phi[i] + c.dt * force
        self.clipped_fraction = float(np.mean((updated < 0) | (updated > 1)))
        self.phi = np.clip(updated, 0, 1)
        for i, cell_id in enumerate(self.ids):
            event = self.divisions.get(int(cell_id))
            if event is not None:
                self.phi[i] = self._project_volume(self.phi[i], event["volume"])

    def divide(self, index, direction=None):
        """Initiate cytokinesis without changing occupancy, cell count, or contacts.

        Each active mother reserves one future cell slot. Explicit directions are
        available for controlled experiments; the default follows cell shape.
        """
        c = self.config
        parent = int(self.ids[index])
        if parent in self.divisions or len(self.phi) + len(self.divisions) >= c.max_cells:
            return False
        direction = self.division_direction(index) if direction is None else np.asarray(direction, dtype=float)
        if direction.shape != (3,):
            raise ValueError("division direction must have three components")
        norm = np.linalg.norm(direction)
        if not np.isfinite(norm) or norm == 0:
            raise ValueError("division direction must be finite and nonzero")
        direction = direction / norm
        center = self.centers()[index]
        relative = self.xyz - center[:, None, None, None]
        projection = np.einsum("d,dijk->ijk", direction, relative)
        mother = occupancy(self.phi[index])
        lo, hi = float(projection.min()), float(projection.max())
        for _ in range(35):
            offset = (lo + hi) / 2
            weight = 0.5 * (1 + np.tanh((projection - offset) / c.interface_width))
            if np.sum(mother * weight) > 0.5 * np.sum(mother):
                lo = offset
            else:
                hi = offset
        radial = np.sqrt(np.maximum(np.sum(relative**2, axis=0) - projection**2, 0))
        radius = float(np.max(radial[self.phi[index] >= 0.5]))
        self.divisions[parent] = {
            "start": self.time, "axis": direction.tolist(), "offset": offset,
            "radius": radius, "volume": float(self.volumes()[index]),
        }
        self.due[index] = np.inf
        self.lineage[parent]["division_start"] = self.time
        self.lineage[parent]["division_axis"] = direction.tolist()
        return True

    def _cleavage_fields(self, index, event):
        axial, _ = self._division_coordinates(index, event)
        mother = occupancy(self.phi[index])
        weight = 0.5 * (1 + np.tanh(axial / self.config.interface_width))
        first_h = mother * weight
        second_h = mother - first_h
        inverse = lambda h: 0.5 - np.sin(np.arcsin(np.clip(1 - 2 * h, -1, 1)) / 3)
        daughters = np.stack([inverse(first_h), inverse(second_h)]).astype(np.float32)
        # Test the entire plane band, not a single center voxel. This is a
        # resolved neck criterion; coarse grids may require refinement to pass.
        neck_band = np.abs(axial) <= max(self.config.interface_width * 0.5, self.dx * 0.5)
        neck = float(self.phi[index][neck_band].max())
        overlap = float(np.sum(daughters[0]**2 * daughters[1]**2, dtype=np.float64)
                        * self.dx**3 / event["volume"])
        fraction = float(first_h.sum(dtype=np.float64) / mother.sum(dtype=np.float64))
        return daughters, neck, overlap, fraction

    def _finish_divisions(self):
        c = self.config
        for parent, event in list(self.divisions.items()):
            if self.time - event["start"] < c.cytokinesis_duration:
                continue
            index = int(np.flatnonzero(self.ids == parent)[0])
            daughters, neck, overlap, fraction = self._cleavage_fields(index, event)
            if neck > c.neck_threshold or overlap > c.division_overlap_tolerance or not 0.1 < fraction < 0.9:
                # Do not force an unresolved cell to split just because a timer expired.
                continue
            self._complete_division(index, daughters, fraction, neck, overlap)

    def _complete_division(self, index, daughters, fraction, neck, overlap):
        c = self.config
        old_graph = normalized_graph(self.contacts()[0], c.graph_contact_cutoff)
        old_ids = self.ids.tolist()
        old_a, old_h = self.activator.copy(), self.inhibitor.copy()
        prolongation = cleavage_prolongation(len(self.phi), index)
        parent = int(self.ids[index])
        first_id, second_id = self.next_id, self.next_id + 1
        self.next_id += 2
        self.phi = np.concatenate([self.phi[:index], self.phi[index + 1:], daughters])
        target = self.target[index]
        # Preserve the mother target and the same relative volume error in both
        # daughters, avoiding an imposed 50:50 pressure jump in asymmetric lobes.
        self.target = np.r_[np.delete(self.target, index), target * fraction, target * (1 - fraction)]
        fate = self.fate[index]
        perturbation = self.fate_rng.normal(scale=c.partition_noise) if c.differentiation else 0.0
        self.fate = np.r_[np.delete(self.fate, index), fate + 2 * perturbation * (1 - fraction),
                         fate - 2 * perturbation * fraction]
        for name in ("activator", "inhibitor"):
            values = getattr(self, name)
            parent_value = values[index]
            draw = self.signal_rng.normal(scale=c.signal_partition_noise) if c.signaling else 0.0
            variation = parent_value * np.clip(draw, -.1, .1)
            inherited = [parent_value + 2 * variation * (1 - fraction), parent_value - 2 * variation * fraction]
            setattr(self, name, np.r_[np.delete(values, index), inherited])
        self.polarity = prolongation @ self.polarity
        self.ids = np.r_[np.delete(self.ids, index), first_id, second_id]
        self.parents = np.r_[np.delete(self.parents, index), parent, parent]
        self.due = np.r_[np.delete(self.due, index), self.time + self._cycle(), self.time + self._cycle()]
        self.lineage[parent].update(division=self.time, neck_at_abscission=neck,
                                    overlap_at_abscission=overlap)
        self.lineage.extend({"id": child, "parent": parent, "birth": self.time, "division": None}
                            for child in (first_id, second_id))
        del self.divisions[parent]
        new_graph = normalized_graph(self.contacts()[0], c.graph_contact_cutoff)
        self.graph_events.append({"time": self.time, "parent": parent,
                                  "ids_before": old_ids, "ids_after": self.ids.tolist(),
                                  "weights_before": old_graph.weights.tolist(),
                                  "weights_after": new_graph.weights.tolist(),
                                  "activator_before": old_a.tolist(), "inhibitor_before": old_h.tolist(),
                                  "activator_after": self.activator.tolist(), "inhibitor_after": self.inhibitor.tolist(),
                                  "activator_partition_jump": (self.activator - prolongation @ old_a).tolist(),
                                  "inhibitor_partition_jump": (self.inhibitor - prolongation @ old_h).tolist(),
                                  "before": stability(old_graph, c.signal_beta, c.signal_da, c.signal_dh),
                                  "after": stability(new_graph, c.signal_beta, c.signal_da, c.signal_dh),
                                  "prolongation": prolongation.tolist(),
                                  "mode_transfer": mode_transfer(old_graph, new_graph, prolongation)})

    def step(self):
        contact, exposure = self.contacts()
        c = self.config
        graph = normalized_graph(contact, c.graph_contact_cutoff)
        if c.signaling:
            self.activator, self.inhibitor = integrate(self.activator, self.inhibitor, graph, c.dt,
                                                       c.signal_beta, c.signal_da, c.signal_dh)
        if c.polarity_enabled:
            cue = exposure_cue(self.phi, self.dx)
            self.polarity = evolve(self.polarity, cue, self.activator, graph.delta, c.dt,
                                   c.polarity_rate, c.polarity_alignment, c.polarity_decay)
        self.update_fate(contact, exposure)
        self.mechanical_step()
        self.step_number += 1
        self.time = self.step_number * self.config.dt
        self._finish_divisions()
        while len(self.phi) + len(self.divisions) < self.config.max_cells and np.any(self.due <= self.time):
            self.divide(int(np.argmin(self.due)))
        if not np.all(np.isfinite(self.phi)) or not np.all(np.isfinite(self.fate)):
            raise FloatingPointError("non-finite state; reduce time step")

    def metrics(self):
        h = occupancy(self.phi)
        # A capped diffuse union, independent of arbitrary cell-center covariance.
        union = np.minimum(h.sum(axis=0), 1)
        weight = union.sum()
        centroid = np.einsum("ijk,dijk->d", union, self.xyz) / weight
        offset = (self.xyz - centroid[:, None, None, None]).reshape(3, -1)
        covariance = (offset * union.ravel()) @ offset.T / weight
        eigen = np.linalg.eigvalsh(covariance)
        contact, exposure = self.contacts()
        self.last_contact, self.last_exposure = contact, exposure
        volumes = self.volumes()
        fate_a = self.fate > self.config.fate_threshold
        fate_b = self.fate < -self.config.fate_threshold
        centers = self.centers()
        separation = None
        if np.any(fate_a) and np.any(fate_b):
            ca = np.average(centers[fate_a], axis=0, weights=volumes[fate_a])
            cb = np.average(centers[fate_b], axis=0, weights=volumes[fate_b])
            separation = float(np.linalg.norm(ca - cb) / np.sqrt(np.trace(covariance)))
        edge = np.zeros_like(union, dtype=bool)
        edge[[0, -1], :, :] = True
        edge[:, [0, -1], :] = True
        edge[:, :, [0, -1]] = True
        return {
            "time": self.time, "cells": len(self.phi),
            "dividing_cells": len(self.divisions),
            "overdue_divisions": sum(self.time - e["start"] > 3 * self.config.cytokinesis_duration
                                     for e in self.divisions.values()),
            "volume_projection_max": self.volume_projection_max,
            "fate_a": int(fate_a.sum()), "fate_b": int(fate_b.sum()),
            "uncommitted": int((~(fate_a | fate_b)).sum()),
            "total_volume": float(volumes.sum()), "target_volume": float(self.target.sum()),
            "relative_volume_error": float((volumes.sum() - self.target.sum()) / self.target.sum()),
            "max_cell_volume_error": float(np.max(np.abs(volumes / self.target - 1))),
            "axis_ratio": float(np.sqrt(eigen[-1] / max(eigen[0], 1e-12))),
            "asphericity": float(1.5 * np.sum((eigen - eigen.mean())**2) / np.sum(eigen)**2),
            "fate_separation": separation,
            "mean_exposure": float(exposure.mean()),
            "activator_std": float(self.activator.std()),
            "inhibitor_std": float(self.inhibitor.std()),
            "mean_polarity": float(np.linalg.norm(self.polarity, axis=1).mean()),
            "unstable_graph_modes": len(stability(normalized_graph(contact, self.config.graph_contact_cutoff),
                                                   self.config.signal_beta, self.config.signal_da,
                                                   self.config.signal_dh)["unstable_modes"]),
            "boundary_occupancy": float(union[edge].max()),
            "clipped_fraction": self.clipped_fraction,
            "min_radius_grid_cells": float(np.min((3 * volumes / (4 * np.pi))**(1 / 3)) / self.dx),
        }

    def surfaces(self, max_points=1400):
        """Interpolate phi=0.5 edge crossings for visualization (not physical measurement)."""
        cells = []
        volumes = self.volumes()
        centers = self.centers()
        for i, field in enumerate(self.phi):
            points = []
            for axis in range(3):
                left, right = [slice(None)] * 3, [slice(None)] * 3
                left[axis], right[axis] = slice(None, -1), slice(1, None)
                a, b = field[tuple(left)], field[tuple(right)]
                indices = np.argwhere((a >= 0.5) != (b >= 0.5))
                if not len(indices):
                    continue
                ix = tuple(indices.T)
                coords = (indices.astype(float) + 0.5) * self.dx - self.config.extent
                coords[:, axis] += (0.5 - a[ix]) / (b[ix] - a[ix]) * self.dx
                points.append(coords)
            points = np.concatenate(points) if points else np.empty((0, 3))
            if len(points) > max_points:
                points = points[np.linspace(0, len(points) - 1, max_points).astype(int)]
            cells.append({"id": int(self.ids[i]), "parent": int(self.parents[i]),
                          "fate": float(self.fate[i]), "volume": float(volumes[i]),
                          "dividing": int(self.ids[i]) in self.divisions,
                          "activator": float(self.activator[i]), "inhibitor": float(self.inhibitor[i]),
                          "polarity": self.polarity[i].tolist(), "center": centers[i].tolist(),
                          "points": np.round(points, 4).tolist()})
        return cells

    def checkpoint(self, path):
        """Save the full state, including RNG state, for exact continuation."""
        import json
        np.savez_compressed(path, phi=self.phi, target=self.target, fate=self.fate,
                            activator=self.activator, inhibitor=self.inhibitor, polarity=self.polarity,
                            ids=self.ids, parents=self.parents, due=self.due,
                            metadata=json.dumps({"config": asdict(self.config), "time": self.time,
                                "step_number": self.step_number, "next_id": self.next_id,
                                "initial_volume": self.initial_volume, "lineage": self.lineage,
                                "rng": self.rng.bit_generator.state,
                                "fate_rng": self.fate_rng.bit_generator.state,
                                "signal_rng": self.signal_rng.bit_generator.state,
                                "schema_version": 3, "divisions": self.divisions, "graph_events": self.graph_events,
                                "volume_projection_max": self.volume_projection_max,
                                "clipped_fraction": self.clipped_fraction}))

    @classmethod
    def restore(cls, path):
        import json
        with np.load(path, allow_pickle=False) as data:
            meta = json.loads(str(data["metadata"]))
            if meta.get("schema_version") != 3:
                raise ValueError("checkpoint uses a previous model without graph signaling/polarity; start a new run")
            sim = cls(Config(**meta["config"]))
            for key in ("phi", "target", "fate", "ids", "parents", "due", "activator", "inhibitor", "polarity"):
                setattr(sim, key, data[key].copy())
            for key in ("time", "step_number", "next_id", "initial_volume", "lineage"):
                setattr(sim, key, meta[key])
            sim.rng.bit_generator.state = meta["rng"]
            sim.fate_rng.bit_generator.state = meta["fate_rng"]
            sim.signal_rng.bit_generator.state = meta["signal_rng"]
            sim.graph_events = meta["graph_events"]
            sim.divisions = {int(k): v for k, v in meta.get("divisions", {}).items()}
            sim.volume_projection_max = meta.get("volume_projection_max", 0.0)
            sim.clipped_fraction = meta.get("clipped_fraction", 0.0)
        return sim
