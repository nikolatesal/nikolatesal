# Model notes

## Fluid model

The case uses OpenLB's two-lattice free-energy multiphase formulation. Lattice 1
stores total density and momentum with forced BGK dynamics; lattice 2 stores the
order parameter `phi` with free-energy BGK dynamics. The coupling order is:

1. chemical-potential coupling from the order parameter,
2. fluid collide-and-stream,
3. force coupling back to the hydrodynamic lattice.

The initial condition is the union of two smooth spherical droplets resting on the
lower substrate. A narrow smooth bridge is added between the droplets to avoid a
purely grid-dependent first contact time.

## Wetting

The lower wall is assigned a superhydrophobic wetting parameter equivalent to a
large apparent contact angle. The particle model stores its own hydrophilic
contact angle and converts it to an adhesion coefficient that attracts the sphere
toward liquid-rich cells.

## Particle model

The spherical solid particle is represented by a Lagrangian rigid sphere with
position, velocity, radius, effective mass, and contact-angle-derived adhesion.
Each time step:

1. local fluid velocity and order parameter are interpolated at the particle
   center,
2. Stokes-like drag relaxes the particle toward the local liquid velocity,
3. a capillary adhesion term pulls the particle toward the liquid center of mass
   while the sampled order parameter indicates liquid contact,
4. a soft collision prevents penetration into the substrate,
5. particle position is advanced by explicit Euler integration.

This coupling is intentionally localized. For production immersed-boundary work,
replace the Lagrangian force closure with OpenLB's resolved-particle module while
keeping the same initialization, diagnostics, and configuration structure.

## Diagnostics

The CSV output records:

- simulation step and lattice time,
- liquid center of mass,
- liquid volume proxy,
- particle position and velocity,
- particle-liquid lag distance,
- a take-off flag indicating the merged droplet no longer intersects the wall
  sampling band.
