# Package index

## Package overview

Top-level documentation entry.

- [`navieRstokes-package`](https://robustecologies.github.io/navieRstokes/reference/navieRstokes-package.md)
  [`navieRstokes`](https://robustecologies.github.io/navieRstokes/reference/navieRstokes-package.md)
  : navieRstokes: 2D incompressible Navier-Stokes equation solver

## Simulation core

Single entry point for a 2D incompressible Navier-Stokes integration on
a uniform staggered grid.

- [`simulate_navier_stokes()`](https://robustecologies.github.io/navieRstokes/reference/simulate_navier_stokes.md)
  : Simulate 2D incompressible Navier-Stokes equations with the finite
  difference method

## Initial conditions

Analytical velocity fields for benchmarks and reproducible experiments.

- [`vortex_ic()`](https://robustecologies.github.io/navieRstokes/reference/vortex_ic.md)
  [`taylor_green_ic()`](https://robustecologies.github.io/navieRstokes/reference/vortex_ic.md)
  : Taylor-Green vortex initial condition
- [`shear_layer_ic()`](https://robustecologies.github.io/navieRstokes/reference/shear_layer_ic.md)
  : Shear layer initial condition (Kelvin-Helmholtz instability)
- [`rotating_ic()`](https://robustecologies.github.io/navieRstokes/reference/rotating_ic.md)
  : Rotating flow initial condition
- [`single_vortex_ic()`](https://robustecologies.github.io/navieRstokes/reference/single_vortex_ic.md)
  : Single Gaussian vortex initial condition
- [`two_vortex_ic()`](https://robustecologies.github.io/navieRstokes/reference/two_vortex_ic.md)
  : Two counter-rotating vortices initial condition

## Forcing functions

Body-force constructors injected into the momentum equation.

- [`pressure_gradient_force()`](https://robustecologies.github.io/navieRstokes/reference/pressure_gradient_force.md)
  : Constant pressure gradient forcing (Poiseuille/channel flow)
- [`oscillatory_force()`](https://robustecologies.github.io/navieRstokes/reference/oscillatory_force.md)
  : Time-periodic forcing
- [`taylor_green_force()`](https://robustecologies.github.io/navieRstokes/reference/taylor_green_force.md)
  : Taylor-Green vortex forcing function
- [`localized_vortex_force()`](https://robustecologies.github.io/navieRstokes/reference/localized_vortex_force.md)
  : Spatially localized forcing (vortex generation)

## Diagnostics

Post-processing operators on the velocity and pressure fields.

- [`compute_vorticity()`](https://robustecologies.github.io/navieRstokes/reference/compute_vorticity.md)
  : Compute vorticity from velocity field

## Visualization

ggplot2-based field plots and animation helpers.

- [`flow()`](https://robustecologies.github.io/navieRstokes/reference/flow.md)
  : Create animation from navieRstokes simulation results
- [`plot(`*`<navieRstokes>`*`)`](https://robustecologies.github.io/navieRstokes/reference/plot.navieRstokes.md)
  : Plot method for navieRstokes objects

## Interactive dashboard

Shiny dashboard for exploratory parameter sweeps.

- [`shiny_navieRstokes()`](https://robustecologies.github.io/navieRstokes/reference/shiny_navieRstokes.md)
  : Launch interactive Shiny dashboard for Navier-Stokes simulations

## S3 methods

Print, summary, plot and flow methods for the navieRstokes class.

- [`print(`*`<navieRstokes>`*`)`](https://robustecologies.github.io/navieRstokes/reference/print.navieRstokes.md)
  : Print method for navieRstokes objects
- [`print(`*`<summary.navieRstokes>`*`)`](https://robustecologies.github.io/navieRstokes/reference/print.summary.navieRstokes.md)
  : Print method for summary.navieRstokes objects
- [`summary(`*`<navieRstokes>`*`)`](https://robustecologies.github.io/navieRstokes/reference/summary.navieRstokes.md)
  : Summary method for navieRstokes objects
