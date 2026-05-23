# navieRstokes: 2D incompressible Navier-Stokes equation solver

The navieRstokes package provides a comprehensive finite difference
method solver for the 2D incompressible Navier-Stokes equations using
Chorin's fractional-step projection method. The package implements
efficient Rcpp-based computational kernels for optimal performance.

## Note

This package is the original creation of the author in all conceptual,
theoretical, and design aspects. Implementation was assisted by
Anthropic's Claude Code v.2 (Opus 4.5) to streamline package
development. All original ideas, creativity, and scientific
contributions belong to the author, who maintains full responsibility
for the package's correctness and reliability. While all code has been
thoroughly reviewed and tested, users are encouraged to report any
issues through the package's issue tracker.

## Main functions

- [`simulate_navier_stokes`](https://robustecologies.github.io/navieRstokes/reference/simulate_navier_stokes.md):
  Main solver function

- [`plot.navieRstokes`](https://robustecologies.github.io/navieRstokes/reference/plot.navieRstokes.md):
  S3 plot method for visualization

- [`flow`](https://robustecologies.github.io/navieRstokes/reference/flow.md):
  Generate animations from simulation results

- [`compute_vorticity`](https://robustecologies.github.io/navieRstokes/reference/compute_vorticity.md):
  Compute vorticity from velocity field

## Initial conditions

- [`vortex_ic`](https://robustecologies.github.io/navieRstokes/reference/vortex_ic.md)
  (alias:
  [`taylor_green_ic`](https://robustecologies.github.io/navieRstokes/reference/vortex_ic.md)):
  Taylor-Green vortex

- [`shear_layer_ic`](https://robustecologies.github.io/navieRstokes/reference/shear_layer_ic.md):
  Shear layer (Kelvin-Helmholtz instability)

- [`rotating_ic`](https://robustecologies.github.io/navieRstokes/reference/rotating_ic.md):
  Solid-body rotation

- [`two_vortex_ic`](https://robustecologies.github.io/navieRstokes/reference/two_vortex_ic.md):
  Two counter-rotating vortices

- [`single_vortex_ic`](https://robustecologies.github.io/navieRstokes/reference/single_vortex_ic.md):
  Single Gaussian vortex

## Forcing functions

- [`pressure_gradient_force`](https://robustecologies.github.io/navieRstokes/reference/pressure_gradient_force.md):
  Constant pressure gradient (channel flow)

- [`oscillatory_force`](https://robustecologies.github.io/navieRstokes/reference/oscillatory_force.md):
  Time-periodic forcing

- [`localized_vortex_force`](https://robustecologies.github.io/navieRstokes/reference/localized_vortex_force.md):
  Spatially localized rotational forcing

- [`taylor_green_force`](https://robustecologies.github.io/navieRstokes/reference/taylor_green_force.md):
  Forcing to maintain Taylor-Green vortex

## Features

- Complete 2D Navier-Stokes solver with Chorin's projection method

- Multiple boundary condition types (Dirichlet, Neumann, periodic)

- Iterative pressure solvers (Jacobi, SOR)

- External forcing support

- Comprehensive diagnostics and CFL monitoring

- Rcpp-optimized computational kernels

- Interactive Shiny dashboard for simulations

## See also

Useful links:

- <https://github.com/robustecologies/navieRstokes>

- <https://robustecologies.github.io/navieRstokes>

- Report bugs at
  <https://github.com/robustecologies/navieRstokes/issues>

## Author

**Maintainer**: Pablo Almaraz <pablo.almaraz@csic.es>
([ORCID](https://orcid.org/0000-0003-1416-2695))
