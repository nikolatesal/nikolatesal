# OpenLB droplet-particle coalescence jump case

This repository contains an OpenLB-based 3-D case for simulating two equal droplets
coalescing on a superhydrophobic substrate and transporting a spherical hydrophilic
solid particle during the coalescence-induced jump.

The implementation is intended to be copied or built against an OpenLB source tree.
It uses a two-lattice free-energy multiphase model for the liquid/gas interface, a
wetting wall for the superhydrophobic lower substrate, and a lightweight
Lagrangian spherical-particle model that samples the OpenLB velocity/phase field to
advance the particle state and exports the particle trajectory for ParaView.

## Repository layout

- `src/droplet_particle_jump3d.cpp` — OpenLB case implementation.
- `configs/default.ini` — physically grouped lattice-unit parameters.
- `scripts/run_case.sh` — builds and runs the case against an OpenLB tree.
- `scripts/postprocess.py` — derives jump height, take-off time, particle lag, and
  liquid-particle center-of-mass diagnostics from CSV output.
- `docs/model.md` — model equations, assumptions, and extension points.

## Quick start

```bash
export OPENLB_DIR=/path/to/openlb
./scripts/run_case.sh configs/default.ini
python3 scripts/postprocess.py tmp/droplet_particle_jump3d/diagnostics.csv
```

The run writes VTK/VTM files under `tmp/droplet_particle_jump3d/` and a CSV file
named `diagnostics.csv` with droplet and particle metrics.

## Notes

OpenLB releases move APIs over time. The case follows the OpenLB free-energy
examples (`multiComponent/contactAngle3d`) and keeps all case-specific coupling in
this repository so that adapting to a newer OpenLB version is localized to
`src/droplet_particle_jump3d.cpp`.
