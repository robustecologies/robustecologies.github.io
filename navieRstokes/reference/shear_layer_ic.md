# Shear layer initial condition (Kelvin-Helmholtz instability)

Creates a shear layer with small perturbations that trigger the
Kelvin-Helmholtz instability, leading to vortex roll-up.

## Usage

``` r
shear_layer_ic(x, y, U0 = 1, delta = 0.05, epsilon = 0.05, lx = 1, ly = 1)
```

## Arguments

- x:

  Numeric. x-coordinate

- y:

  Numeric. y-coordinate

- U0:

  Numeric. Velocity amplitude (default = 1.0)

- delta:

  Numeric. Shear layer thickness (default = 0.05)

- epsilon:

  Numeric. Perturbation amplitude (default = 0.05)

- lx:

  Numeric. Domain length in x-direction (default = 1.0)

- ly:

  Numeric. Domain length in y-direction (default = 1.0)

## Value

List with components u and v

## See also

[`vortex_ic`](https://robustecologies.github.io/navieRstokes/reference/vortex_ic.md),
[`rotating_ic`](https://robustecologies.github.io/navieRstokes/reference/rotating_ic.md),
[`single_vortex_ic`](https://robustecologies.github.io/navieRstokes/reference/single_vortex_ic.md),
[`two_vortex_ic`](https://robustecologies.github.io/navieRstokes/reference/two_vortex_ic.md),
[`simulate_navier_stokes`](https://robustecologies.github.io/navieRstokes/reference/simulate_navier_stokes.md)

## Examples

``` r
if (FALSE) { # \dontrun{
result <- simulate_navier_stokes(
  nx = 128, ny = 128, lx = 1.0, ly = 1.0,
  dt = 0.001, nt = 2000, nu = 0.0001,
  initial_condition = shear_layer_ic,
  boundary_condition = list(type = "periodic")
)
} # }
```
