# Taylor-Green vortex initial condition

Creates a Taylor-Green vortex, a classical test case for Navier-Stokes
solvers. The flow is periodic and has an exact analytical solution for
the decay of kinetic energy.

## Usage

``` r
vortex_ic(x, y, A = 1, lx = 1, ly = 1)

taylor_green_ic(x, y, A = 1, lx = 1, ly = 1)
```

## Arguments

- x:

  Numeric. x-coordinate

- y:

  Numeric. y-coordinate

- A:

  Numeric. Amplitude of the vortex (default = 1.0)

- lx:

  Numeric. Domain length in x-direction (default = 1.0)

- ly:

  Numeric. Domain length in y-direction (default = 1.0)

## Value

List with components u and v

## Details

The velocity field is: \$\$u(x,y) = A \sin(2\pi x / l_x) \cos(2\pi y /
l_y)\$\$ \$\$v(x,y) = -A \cos(2\pi x / l_x) \sin(2\pi y / l_y)\$\$

This is divergence-free by construction. The kinetic energy decays as:
\$\$E(t) = E_0 \exp(-2\nu (k_x^2 + k_y^2) t)\$\$ where \\k_x =
2\pi/l_x\\ and \\k_y = 2\pi/l_y\\. When \\l_x = l_y\\, this simplifies
to \\E(t) = E_0 \exp(-4\nu k^2 t)\\ with \\k = 2\pi/l_x\\.

Note: The sign convention used here differs from some textbooks (which
use \\u = -\cos(kx)\sin(ky)\\), but both are valid divergence-free
fields. The package's
[`taylor_green_force`](https://robustecologies.github.io/navieRstokes/reference/taylor_green_force.md)
matches this convention.

## See also

`taylor_green_ic` (alias),
[`single_vortex_ic`](https://robustecologies.github.io/navieRstokes/reference/single_vortex_ic.md),
[`two_vortex_ic`](https://robustecologies.github.io/navieRstokes/reference/two_vortex_ic.md),
[`shear_layer_ic`](https://robustecologies.github.io/navieRstokes/reference/shear_layer_ic.md),
[`rotating_ic`](https://robustecologies.github.io/navieRstokes/reference/rotating_ic.md),
[`taylor_green_force`](https://robustecologies.github.io/navieRstokes/reference/taylor_green_force.md),
[`simulate_navier_stokes`](https://robustecologies.github.io/navieRstokes/reference/simulate_navier_stokes.md)

## Examples

``` r
if (FALSE) { # \dontrun{
result <- simulate_navier_stokes(
  nx = 64, ny = 64, lx = 1.0, ly = 1.0,
  dt = 0.001, nt = 500, nu = 0.01,
  initial_condition = vortex_ic,
  boundary_condition = list(type = "periodic")
)
} # }
```
