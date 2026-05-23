# Constant pressure gradient forcing (Poiseuille/channel flow)

Applies a constant pressure gradient in the x-direction, driving flow
through a channel. Useful for simulating Poiseuille flow with no-slip
walls at top and bottom.

## Usage

``` r
pressure_gradient_force(x, y, t, dp_dx = 0.1, rho = 1)
```

## Arguments

- x:

  Numeric. x-coordinate (not used, for interface compatibility)

- y:

  Numeric. y-coordinate (not used, for interface compatibility)

- t:

  Numeric. Time (not used, for interface compatibility)

- dp_dx:

  Numeric. Pressure gradient in x-direction (default = 0.1)

- rho:

  Numeric. Fluid density (default = 1.0)

## Value

List with components fx and fy

## Details

The forcing term represents a constant body force equivalent to a
pressure gradient driving the flow. The x-component of forcing is
computed as: \$\$f_x = -\frac{1}{\rho}\frac{dp}{dx}\$\$ with \\f_y =
0\\. For plane Poiseuille flow between parallel plates at \\y = 0\\ and
\\y = h\\ with no-slip boundary conditions, the steady-state analytical
solution is the parabolic velocity profile: \$\$u(y) =
-\frac{1}{2\mu}\frac{dp}{dx}y(h - y)\$\$ where \\\mu = \rho\nu\\ is the
dynamic viscosity. The maximum velocity occurs at the channel centerline
\\y = h/2\\.

**Sign convention:** With the default `dp_dx = 0.1` (positive), the
pressure increases in the +x direction, so the resulting force \\f_x =
-dp\\dx / \rho = -0.1\\ drives flow in the **negative** x-direction. To
drive flow in the positive x-direction, pass a negative value for
`dp_dx` (e.g., `dp_dx = -0.1`).

## References

Tritton, D. J. (1988). *Physical fluid dynamics* (2nd ed.). Oxford
University Press. ISBN: 978-0-19-854493-7.

Batchelor, G. K. (1967). *An introduction to fluid dynamics*. Cambridge
University Press. ISBN: 978-0-521-66396-0. DOI:
[doi:10.1017/CBO9780511800955](https://doi.org/10.1017/CBO9780511800955)

## See also

[`oscillatory_force`](https://robustecologies.github.io/navieRstokes/reference/oscillatory_force.md),
[`taylor_green_force`](https://robustecologies.github.io/navieRstokes/reference/taylor_green_force.md),
[`localized_vortex_force`](https://robustecologies.github.io/navieRstokes/reference/localized_vortex_force.md),
[`simulate_navier_stokes`](https://robustecologies.github.io/navieRstokes/reference/simulate_navier_stokes.md)

## Examples

``` r
if (FALSE) { # \dontrun{
result <- simulate_navier_stokes(
  nx = 32, ny = 64, lx = 1.0, ly = 1.0,
  dt = 0.0005, nt = 500, nu = 0.1,
  initial_condition = function(x, y) list(u = 0, v = 0),
  forcing_function = pressure_gradient_force,
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
