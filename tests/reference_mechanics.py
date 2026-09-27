"""Frozen pre-optimization kernel for trajectory-equivalence regression tests.

Kept separately from production code: changing physics requires an explicit
review of this oracle, not silently updating it with the optimized kernel.
"""
import numpy as np
from scipy.ndimage import laplace
from embryo.model import Simulation, occupancy
from embryo.polarity import tension_field, flux_divergence


class ReferenceSimulation(Simulation):
    def centers(self):
        h = occupancy(self.phi)
        return np.einsum("nijk,dijk->nd", h, self.xyz) * self.dx**3 / self.volumes()[:, None]


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


    def mechanical_step(self):
        c = self.config
        phi = self.phi
        shell2 = (phi * (1 - phi))**2
        derivative = 2 * phi * (1 - phi) * (1 - 2 * phi)
        fate = np.tanh(self.fate)
        # Whole-valued parameters (e.g. JSON 4 instead of 4.0) still need
        # floating-point storage for the in-place mechanical feedback below.
        adhesion = np.full((len(phi), len(phi)), c.adhesion, dtype=float)
        tension = np.full(len(phi), c.surface_tension, dtype=float)
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
