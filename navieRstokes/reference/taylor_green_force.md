# Taylor-Green vortex forcing function

Generates the forcing term for the Taylor-Green vortex flow, which
produces an exact solution to the Navier-Stokes equations. This is a
steady forcing that counteracts viscous decay, maintaining the vortex
structure over time.

## Usage

``` r
taylor_green_force(x, y, t, A = 1, k = 2 * pi, nu = 0.01, lx = 1, ly = 1)
```

## Arguments

- x:

  Numeric. x-coordinate(s) in the domain

- y:

  Numeric. y-coordinate(s) in the domain

- t:

  Numeric. Current simulation time

- A:

  Numeric. Amplitude of the vortex (default: 1.0)

- k:

  Numeric. Wave number (deprecated; use `lx`, `ly` instead). Retained
  for backward compatibility. If provided, overrides domain-derived wave
  numbers with \\k_x = k_y = k\\. Default: \\2\pi\\.

- nu:

  Numeric. Kinematic viscosity (default: 0.01)

- lx:

  Numeric. Domain length in x-direction (default: 1.0). Used to compute
  \\k_x = 2\pi/l_x\\.

- ly:

  Numeric. Domain length in y-direction (default: 1.0). Used to compute
  \\k_y = 2\pi/l_y\\.

## Value

A list with components:

- fx:

  Numeric. Forcing in x-direction

- fy:

  Numeric. Forcing in y-direction

## Details

This forcing function is designed to match the initial condition from
[`vortex_ic`](https://robustecologies.github.io/navieRstokes/reference/vortex_ic.md)
(also aliased as
[`taylor_green_ic`](https://robustecologies.github.io/navieRstokes/reference/vortex_ic.md)):
\$\$u(x,y,t) = A\sin(k_x x)\cos(k_y y) e^{-\nu(k_x^2+k_y^2) t}\$\$
\$\$v(x,y,t) = -A\cos(k_x x)\sin(k_y y) e^{-\nu(k_x^2+k_y^2) t}\$\$
where \\k_x = 2\pi/l_x\\ and \\k_y = 2\pi/l_y\\.

The forcing compensates for viscous decay: \$\$f_x = \nu(k_x^2 + k_y^2)
u\$\$ \$\$f_y = \nu(k_x^2 + k_y^2) v\$\$

Note: The sign convention here differs from some textbooks (which use
\\u = -\cos(kx)\sin(ky)\\), but matches the package's
[`vortex_ic()`](https://robustecologies.github.io/navieRstokes/reference/vortex_ic.md).
Both are valid divergence-free fields with opposite rotation sense.

## References

Taylor, G. I., & Green, A. E. (1937). Mechanism of the production of
small eddies from large ones. *Proceedings of the Royal Society of
London A*, 158(895), 499-521. DOI:
[doi:10.1098/rspa.1937.0036](https://doi.org/10.1098/rspa.1937.0036)

## See also

[`vortex_ic`](https://robustecologies.github.io/navieRstokes/reference/vortex_ic.md),
[`taylor_green_ic`](https://robustecologies.github.io/navieRstokes/reference/vortex_ic.md),
[`pressure_gradient_force`](https://robustecologies.github.io/navieRstokes/reference/pressure_gradient_force.md),
[`oscillatory_force`](https://robustecologies.github.io/navieRstokes/reference/oscillatory_force.md),
[`localized_vortex_force`](https://robustecologies.github.io/navieRstokes/reference/localized_vortex_force.md),
[`simulate_navier_stokes`](https://robustecologies.github.io/navieRstokes/reference/simulate_navier_stokes.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Taylor-Green vortex with forcing to maintain structure
result <- simulate_navier_stokes(
  nx = 64, ny = 64,
  lx = 1.0, ly = 1.0,
  dt = 0.001, nt = 1000,
  nu = 0.01,
  initial_condition = function(x, y) vortex_ic(x, y, A = 1.0),
  forcing_function = taylor_green_force,
  boundary_condition = list(type = "periodic")
)

# Custom parameters
custom_force <- function(x, y, t) {
  taylor_green_force(x, y, t, A = 0.5, k = 2 * pi, nu = 0.01)
}
par(mfrow = c(2, 2))
plot(result, time_index = 1, plot_type = "vorticity")
plot(result, time_index = 30, plot_type = "vorticity")
plot(result, time_index = 70, plot_type = "vorticity")
plot(result, time_index = 150, plot_type = "vorticity")
} # }
```
