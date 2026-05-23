# Launch interactive Shiny dashboard for Navier-Stokes simulations

Launches an interactive Shiny app using `shinydashboard` that provides
access to the theoretical background of the Navier-Stokes equations and
allows interactive simulation of various fluid dynamics scenarios
including lid-driven cavity, vortex decay, shear layer instability, and
more.

## Usage

``` r
shiny_navieRstokes()
```

## Value

No return value. Launches the Shiny application in the default web
browser.

## Details

The dashboard includes:

- **Theory tab**: Mathematical formulation of the Navier-Stokes
  equations and Chorin's projection method

- **Simulation tabs**: Interactive scenarios including:

  - Lid-driven cavity flow

  - Taylor-Green vortex decay

  - Kelvin-Helmholtz instability

  - Solid body rotation

  - Poiseuille channel flow

- **Interactive controls**: Real-time parameter adjustment for grid
  resolution, viscosity, Reynolds number, time stepping, etc.

- **Visualizations**: Velocity fields, vorticity contours, and pressure
  distributions with customizable color schemes

- **Diagnostics**: Mass conservation errors, CFL numbers, and solver
  convergence metrics

## Note

Requires the `shiny` and `shinydashboard` packages. Install with:
`install.packages(c("shiny", "shinydashboard"))`

## See also

[`simulate_navier_stokes`](https://robustecologies.github.io/navieRstokes/reference/simulate_navier_stokes.md),
[`plot.navieRstokes`](https://robustecologies.github.io/navieRstokes/reference/plot.navieRstokes.md),
[`flow.navieRstokes`](https://robustecologies.github.io/navieRstokes/reference/flow.md),
[`vortex_ic`](https://robustecologies.github.io/navieRstokes/reference/vortex_ic.md),
[`pressure_gradient_force`](https://robustecologies.github.io/navieRstokes/reference/pressure_gradient_force.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Launch the interactive dashboard
shiny_navieRstokes()
} # }
```
