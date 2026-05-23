# Time-periodic forcing

Applies a sinusoidally varying force in x-direction, useful for studying
oscillatory flows.

## Usage

``` r
oscillatory_force(x, y, t, A = 0.5, omega = 2 * pi)
```

## Arguments

- x:

  Numeric. x-coordinate (not used)

- y:

  Numeric. y-coordinate (not used)

- t:

  Numeric. Time

- A:

  Numeric. Force amplitude (default = 0.5)

- omega:

  Numeric. Angular frequency (default = 2\*pi)

## Value

List with components fx and fy

## Details

The forcing term oscillates sinusoidally in time with the form:
\$\$f_x(t) = A \sin(\omega t)\$\$ with \\f_y = 0\\. This forcing is
spatially uniform and drives oscillatory Stokes flow when the Reynolds
number is small. For an infinite domain with periodic boundaries, the
analytical solution for the velocity field is: \$\$u(t) =
\frac{A}{\omega}\left(1 - \cos(\omega t)\right)\$\$ starting from rest.
The Womersley number \\Wo = L\sqrt{\omega/\nu}\\ characterizes the ratio
of oscillatory inertia to viscous forces.

## References

Batchelor, G. K. (1967). *An introduction to fluid dynamics*. Cambridge
University Press. ISBN: 978-0-521-66396-0. DOI:
[doi:10.1017/CBO9780511800955](https://doi.org/10.1017/CBO9780511800955)

Pozrikidis, C. (2016). *Fluid dynamics: Theory, computation, and
numerical simulation* (3rd ed.). Springer. ISBN: 978-1-4899-7990-2. DOI:
[doi:10.1007/978-1-4899-7991-9](https://doi.org/10.1007/978-1-4899-7991-9)

## See also

[`pressure_gradient_force`](https://robustecologies.github.io/navieRstokes/reference/pressure_gradient_force.md),
[`taylor_green_force`](https://robustecologies.github.io/navieRstokes/reference/taylor_green_force.md),
[`localized_vortex_force`](https://robustecologies.github.io/navieRstokes/reference/localized_vortex_force.md),
[`simulate_navier_stokes`](https://robustecologies.github.io/navieRstokes/reference/simulate_navier_stokes.md)

## Examples

``` r
if (FALSE) { # \dontrun{
result <- simulate_navier_stokes(
  nx = 64, ny = 64, lx = 1.0, ly = 1.0,
  dt = 0.001, nt = 1000, nu = 0.01,
  initial_condition = function(x, y) list(u = 0, v = 0),
  forcing_function = oscillatory_force,
  boundary_condition = list(type = "periodic")
)
} # }
```
