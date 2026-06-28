"""LBM/IBM demo: coalescence and rebound of two water droplets on a superhydrophobic wall.

This is a compact research/teaching prototype rather than a quantitatively validated
solver.  It combines:

* D2Q9 lattice-Boltzmann hydrodynamics with a single-relaxation-time BGK collision.
* A conservative diffuse-interface order parameter ``phi`` to represent two droplets
  that merge on a non-wetting wall.
* A simple immersed-boundary (IBM) hydrophilic particle model: a Lagrangian circular
  particle (50 micrometre diameter) is advected by interpolated fluid velocity, while
  a short-range attraction pulls the particle toward the water phase during merger.
* A phenomenological wall-repulsion impulse that mimics superhydrophobic rebound.

Default physical scales are two 200 micrometre droplets and one 50 micrometre
particle.  The mapping to lattice units is intentionally explicit so the script can
be used as a starting point for more detailed calibration.

Example
-------
    python lbm_ibm_coalescence.py --steps 2500 --output frames --save-every 100

The script writes PNG frames when matplotlib is installed and always writes a CSV
trajectory for the hydrophilic particle.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np

try:
    import matplotlib.pyplot as plt
except ImportError:  # plotting is optional; the numerical demo still runs
    plt = None


# D2Q9 discrete velocities and weights.
C = np.array(
    [[0, 0], [1, 0], [0, 1], [-1, 0], [0, -1], [1, 1], [-1, 1], [-1, -1], [1, -1]],
    dtype=np.int8,
)
W = np.array([4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 36, 1 / 36, 1 / 36, 1 / 36])
CS2 = 1.0 / 3.0


@dataclass
class PhysicalScales:
    droplet_diameter_um: float = 200.0
    particle_diameter_um: float = 50.0
    water_density: float = 997.0
    water_kinematic_viscosity: float = 1.0e-6
    surface_tension: float = 0.072
    lattice_nodes_per_droplet: int = 60
    lattice_dt_s: float = 2.0e-7

    @property
    def dx_m(self) -> float:
        return self.droplet_diameter_um * 1.0e-6 / self.lattice_nodes_per_droplet

    @property
    def particle_radius_lu(self) -> float:
        return 0.5 * self.particle_diameter_um / self.droplet_diameter_um * self.lattice_nodes_per_droplet

    @property
    def viscosity_lu(self) -> float:
        return self.water_kinematic_viscosity * self.lattice_dt_s / (self.dx_m**2)


@dataclass
class SimulationConfig:
    nx: int = 220
    ny: int = 140
    steps: int = 2500
    save_every: int = 100
    tau: float = 0.82
    phase_tau: float = 1.0
    interface_width: float = 4.0
    surface_tension_force: float = 0.014
    wall_repulsion: float = 0.022
    wall_decay: float = 8.0
    particle_phase_attraction: float = 0.018
    output: Path = Path("frames")
    seed: int = 7


@dataclass
class Particle:
    position: np.ndarray
    velocity: np.ndarray
    radius: float


def equilibrium(rho: np.ndarray, ux: np.ndarray, uy: np.ndarray) -> np.ndarray:
    """Return D2Q9 equilibrium populations."""
    u2 = ux * ux + uy * uy
    feq = np.empty((9, *rho.shape), dtype=np.float64)
    for i, (cx, cy) in enumerate(C):
        cu = cx * ux + cy * uy
        feq[i] = W[i] * rho * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * u2)
    return feq


def scalar_equilibrium(phi: np.ndarray, ux: np.ndarray, uy: np.ndarray) -> np.ndarray:
    geq = np.empty((9, *phi.shape), dtype=np.float64)
    for i, (cx, cy) in enumerate(C):
        geq[i] = W[i] * phi * (1.0 + 3.0 * (cx * ux + cy * uy))
    return geq


def stream(pop: np.ndarray) -> np.ndarray:
    out = np.empty_like(pop)
    for i, (cx, cy) in enumerate(C):
        out[i] = np.roll(np.roll(pop[i], cx, axis=1), cy, axis=0)
    return out


def initialise(scales: PhysicalScales, cfg: SimulationConfig) -> tuple[np.ndarray, np.ndarray, Particle]:
    y, x = np.mgrid[0 : cfg.ny, 0 : cfg.nx]
    radius = 0.5 * scales.lattice_nodes_per_droplet
    centres = [(cfg.nx * 0.43, radius + 10.0), (cfg.nx * 0.57, radius + 10.0)]
    signed = np.full((cfg.ny, cfg.nx), -1.0e9)
    for cx, cy in centres:
        distance = np.sqrt((x - cx) ** 2 + (y - cy) ** 2)
        signed = np.maximum(signed, radius - distance)
    phi = np.tanh(signed / cfg.interface_width)

    rho = 1.0 + 0.08 * np.maximum(phi, 0.0)
    ux = np.zeros_like(rho)
    uy = np.zeros_like(rho)
    feq = equilibrium(rho, ux, uy)
    geq = scalar_equilibrium(phi, ux, uy)
    particle = Particle(
        position=np.array([cfg.nx * 0.50, radius + 18.0], dtype=np.float64),
        velocity=np.array([0.0, 0.0], dtype=np.float64),
        radius=scales.particle_radius_lu,
    )
    return feq.copy(), geq.copy(), particle


def gradients(a: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return 0.5 * (np.roll(a, -1, axis=1) - np.roll(a, 1, axis=1)), 0.5 * (
        np.roll(a, -1, axis=0) - np.roll(a, 1, axis=0)
    )


def compute_forces(phi: np.ndarray, particle: Particle, cfg: SimulationConfig) -> tuple[np.ndarray, np.ndarray]:
    gx, gy = gradients(phi)
    curvature_proxy = np.roll(phi, 1, 0) + np.roll(phi, -1, 0) + np.roll(phi, 1, 1) + np.roll(phi, -1, 1) - 4.0 * phi
    fx = -cfg.surface_tension_force * curvature_proxy * gx
    fy = -cfg.surface_tension_force * curvature_proxy * gy

    y = np.arange(cfg.ny, dtype=np.float64)[:, None]
    water_near_wall = np.maximum(phi, 0.0)
    fy += cfg.wall_repulsion * np.exp(-y / cfg.wall_decay) * water_near_wall
    return fx, fy


def bilinear(field: np.ndarray, pos: np.ndarray) -> float:
    x, y = pos
    x0 = int(np.floor(x)) % field.shape[1]
    y0 = int(np.floor(y)) % field.shape[0]
    x1 = (x0 + 1) % field.shape[1]
    y1 = min(y0 + 1, field.shape[0] - 1)
    tx = x - np.floor(x)
    ty = y - np.floor(y)
    return float((1 - tx) * (1 - ty) * field[y0, x0] + tx * (1 - ty) * field[y0, x1] + (1 - tx) * ty * field[y1, x0] + tx * ty * field[y1, x1])


def update_particle(particle: Particle, ux: np.ndarray, uy: np.ndarray, phi: np.ndarray, cfg: SimulationConfig) -> None:
    fluid_u = np.array([bilinear(ux, particle.position), bilinear(uy, particle.position)])
    gx, gy = gradients(phi)
    phase_pull = cfg.particle_phase_attraction * np.array([bilinear(gx, particle.position), bilinear(gy, particle.position)])
    particle.velocity = 0.82 * particle.velocity + 0.18 * fluid_u + phase_pull
    particle.position += particle.velocity
    particle.position[0] = np.clip(particle.position[0], particle.radius + 1.0, cfg.nx - particle.radius - 2.0)
    if particle.position[1] < particle.radius + 1.0:
        particle.position[1] = particle.radius + 1.0
        particle.velocity[1] = abs(particle.velocity[1]) * 0.55
    if particle.position[1] > cfg.ny - particle.radius - 2.0:
        particle.position[1] = cfg.ny - particle.radius - 2.0
        particle.velocity[1] *= -0.2


def apply_particle_feedback(fx: np.ndarray, fy: np.ndarray, particle: Particle, ux: np.ndarray, uy: np.ndarray) -> None:
    y, x = np.mgrid[0 : fx.shape[0], 0 : fx.shape[1]]
    r2 = (x - particle.position[0]) ** 2 + (y - particle.position[1]) ** 2
    kernel = np.exp(-r2 / (2.0 * (0.65 * particle.radius) ** 2))
    kernel /= kernel.sum() + 1.0e-30
    desired = particle.velocity
    fx += 0.08 * kernel * (desired[0] - ux)
    fy += 0.08 * kernel * (desired[1] - uy)


def save_frame(step: int, phi: np.ndarray, particle: Particle, cfg: SimulationConfig) -> None:
    if plt is None:
        return
    cfg.output.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.5, 4.2), dpi=140)
    ax.imshow(phi, origin="lower", cmap="coolwarm", vmin=-1, vmax=1)
    ax.add_patch(plt.Circle(particle.position, particle.radius, color="gold", ec="black", lw=0.8))
    ax.axhline(0, color="black", lw=2)
    ax.set_title(f"LBM/IBM coalescence rebound, step {step}")
    ax.set_xlim(0, cfg.nx - 1)
    ax.set_ylim(0, cfg.ny - 1)
    ax.set_xlabel("x [lattice units]")
    ax.set_ylabel("y [lattice units]")
    fig.tight_layout()
    fig.savefig(cfg.output / f"frame_{step:06d}.png")
    plt.close(fig)


def run(scales: PhysicalScales, cfg: SimulationConfig) -> None:
    np.random.seed(cfg.seed)
    f, g, particle = initialise(scales, cfg)
    trajectory: list[tuple[int, float, float, float, float]] = []

    for step in range(cfg.steps + 1):
        rho = np.maximum(f.sum(axis=0), 1.0e-9)
        phi = np.clip(g.sum(axis=0), -1.0, 1.0)
        ux = (f * C[:, 0, None, None]).sum(axis=0) / rho
        uy = (f * C[:, 1, None, None]).sum(axis=0) / rho

        fx, fy = compute_forces(phi, particle, cfg)
        apply_particle_feedback(fx, fy, particle, ux, uy)
        ux_eff = ux + 0.5 * fx / rho
        uy_eff = uy + 0.5 * fy / rho

        f += -(f - equilibrium(rho, ux_eff, uy_eff)) / cfg.tau
        for i, (cx, cy) in enumerate(C):
            f[i] += W[i] * (3.0 * (cx * fx + cy * fy))
        g += -(g - scalar_equilibrium(phi, ux_eff, uy_eff)) / cfg.phase_tau

        f = stream(f)
        g = stream(g)

        # Bounce-back / neutral wetting boundary at bottom and open-ish top boundary.
        f[[2, 5, 6], 0, :] = f[[4, 7, 8], 0, :]
        g[:, 0, :] = scalar_equilibrium(-0.75 * np.ones(cfg.nx), np.zeros(cfg.nx), np.zeros(cfg.nx))
        f[:, -1, :] = f[:, -2, :]
        g[:, -1, :] = g[:, -2, :]

        if step % 2 == 0:
            update_particle(particle, ux_eff, uy_eff, phi, cfg)
        if step % cfg.save_every == 0:
            save_frame(step, phi, particle, cfg)
            trajectory.append((step, particle.position[0], particle.position[1], particle.velocity[0], particle.velocity[1]))

    cfg.output.mkdir(parents=True, exist_ok=True)
    with (cfg.output / "particle_trajectory.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["step", "x_lu", "y_lu", "vx_lu_per_step", "vy_lu_per_step"])
        writer.writerows(trajectory)

    print(f"dx = {scales.dx_m:.3e} m, dt = {scales.lattice_dt_s:.3e} s")
    print(f"water viscosity in lattice units = {scales.viscosity_lu:.4f}; configured tau = {cfg.tau:.3f}")
    print(f"particle radius = {particle.radius:.2f} lu; output written to {cfg.output}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--steps", type=int, default=SimulationConfig.steps)
    parser.add_argument("--save-every", type=int, default=SimulationConfig.save_every)
    parser.add_argument("--output", type=Path, default=SimulationConfig.output)
    parser.add_argument("--nx", type=int, default=SimulationConfig.nx)
    parser.add_argument("--ny", type=int, default=SimulationConfig.ny)
    parser.add_argument("--tau", type=float, default=SimulationConfig.tau)
    parser.add_argument("--droplet-diameter-um", type=float, default=PhysicalScales.droplet_diameter_um)
    parser.add_argument("--particle-diameter-um", type=float, default=PhysicalScales.particle_diameter_um)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    scales = PhysicalScales(
        droplet_diameter_um=args.droplet_diameter_um,
        particle_diameter_um=args.particle_diameter_um,
    )
    cfg = SimulationConfig(nx=args.nx, ny=args.ny, steps=args.steps, save_every=args.save_every, output=args.output, tau=args.tau)
    run(scales, cfg)


if __name__ == "__main__":
    main()
