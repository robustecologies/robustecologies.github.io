# Rotating flow initial condition

Creates a solid-body rotation flow pattern centered at domain center.

## Usage

``` r
rotating_ic(x, y, omega = 2 * pi, x0 = NULL, y0 = NULL, lx = 1, ly = 1)
```

## Arguments

- x:

  Numeric. x-coordinate

- y:

  Numeric. y-coordinate

- omega:

  Numeric. Angular velocity (default = 2\*pi)

- x0:

  Numeric. Center x-coordinate (default = NULL, uses lx/2)

- y0:

  Numeric. Center y-coordinate (default = NULL, uses ly/2)

- lx:

  Numeric. Domain length in x-direction (default = 1.0)

- ly:

  Numeric. Domain length in y-direction (default = 1.0)

## Value

List with components u and v

## Details

The velocity field corresponds to solid-body rotation about the point
\\(x_0, y_0)\\ with angular velocity \\\omega\\. The velocity components
are: \$\$u = -\omega(y - y_0)\$\$ \$\$v = \omega(x - x_0)\$\$ This field
is exactly divergence-free (\\\nabla \cdot \mathbf{u} = 0\\) and has
uniform vorticity \\\zeta = 2\omega\\ throughout the domain. With
no-slip boundary conditions, the rotation decays due to viscous friction
at the walls, providing a classical spin-down problem. The
characteristic decay time scales as \\L^2/\nu\\ where \\L\\ is the
domain size.

## References

Kundu, P. K., Cohen, I. M., & Dowling, D. R. (2015). *Fluid mechanics*
(6th ed.). Academic Press. ISBN: 978-0-12-405935-1.

Batchelor, G. K. (1967). *An introduction to fluid dynamics*. Cambridge
University Press. ISBN: 978-0-521-66396-0. DOI:
[doi:10.1017/CBO9780511800955](https://doi.org/10.1017/CBO9780511800955)

## See also

[`vortex_ic`](https://robustecologies.github.io/navieRstokes/reference/vortex_ic.md),
[`single_vortex_ic`](https://robustecologies.github.io/navieRstokes/reference/single_vortex_ic.md),
[`two_vortex_ic`](https://robustecologies.github.io/navieRstokes/reference/two_vortex_ic.md),
[`shear_layer_ic`](https://robustecologies.github.io/navieRstokes/reference/shear_layer_ic.md),
[`simulate_navier_stokes`](https://robustecologies.github.io/navieRstokes/reference/simulate_navier_stokes.md)

## Examples

``` r
if (FALSE) { # \dontrun{
result <- simulate_navier_stokes(
  nx = 64, ny = 64, lx = 1.0, ly = 1.0,
  dt = 0.0005, nt = 10000, nu = 0.01,
  initial_condition = rotating_ic,
  boundary_condition = list(
    type = "dirichlet",
    values = list(
      u_left = 0, u_right = 0, u_top = 0, u_bottom = 0,
      v_left = 0, v_right = 0, v_top = 0, v_bottom = 0
    )
  )
)
} # }
```
