# LBM/IBM droplet coalescence rebound demo

This repository contains a Python prototype for simulating two 200 μm water droplets
coalescing and rebounding from a superhydrophobic surface while entraining a 50 μm
hydrophilic particle.

## Run

```bash
pip install -r requirements.txt
python lbm_ibm_coalescence.py --steps 2500 --output frames --save-every 100
```

Outputs:

- `frames/frame_*.png`: phase-field snapshots with the immersed particle overlay.
- `frames/particle_trajectory.csv`: particle position and velocity history.

## Model notes

The script is a compact teaching/research starting point, not a fully validated
quantitative solver. It uses a D2Q9 BGK lattice-Boltzmann fluid, a conservative
diffuse-interface scalar for the water/air morphology, a phenomenological
superhydrophobic wall-repulsion force, and a light immersed-boundary coupling for
the hydrophilic particle.

Important physical defaults are exposed in `PhysicalScales`:

- droplet diameter: 200 μm
- particle diameter: 50 μm
- water density: 997 kg/m³
- water kinematic viscosity: 1.0e-6 m²/s
- water surface tension: 0.072 N/m

For publishable predictions, calibrate the lattice-unit force coefficients against
contact angle, Weber number, Ohnesorge number, and restitution data for the exact
surface and particle chemistry.
