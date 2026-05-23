# Single Gaussian vortex initial condition

Creates a single Gaussian vortex for studying vortex decay.

## Usage

``` r
single_vortex_ic(
  x,
  y,
  Gamma = 1,
  sigma = 0.15,
  x0 = NULL,
  y0 = NULL,
  lx = 1,
  ly = 1
)
```

## Arguments

- x:

  Numeric. x-coordinate

- y:

  Numeric. y-coordinate

- Gamma:

  Numeric. Circulation strength (default = 1.0)

- sigma:

  Numeric. Vortex core radius (default = 0.15)

- x0:

  Numeric. Vortex center x (default = NULL, uses lx/2)

- y0:

  Numeric. Vortex center y (default = NULL, uses ly/2)

- lx:

  Numeric. Domain length in x-direction (default = 1.0)

- ly:

  Numeric. Domain length in y-direction (default = 1.0)

## Value

List with components u and v

## Details

This initial condition creates a single Lamb-Oseen (Gaussian) vortex
centered at \\(x_0, y_0)\\. The azimuthal velocity profile is:
\$\$u\_\theta(r) =
\frac{\Gamma}{2\pi\sigma}\frac{r}{\sigma}\exp\left(-\frac{r^2}{2\sigma^2}\right)\$\$
where \\r = \sqrt{(x-x_0)^2 + (y-y_0)^2}\\ is the radial distance from
the vortex center. The velocity field is converted to Cartesian
components: \$\$u = -u\_\theta \frac{(y - y_0)}{r}, \quad v = u\_\theta
\frac{(x - x_0)}{r}\$\$ The Gaussian vortex has several desirable
properties: it is smooth at the center (avoiding the \\1/r\\ singularity
of a point vortex), is exactly divergence-free, and has an analytical
solution for viscous decay. Under viscous diffusion, the core radius
grows as \\\sigma(t) = \sqrt{\sigma_0^2 + 4\nu t}\\ while the peak
vorticity decays.

## References

Lamb, H. (1932). *Hydrodynamics* (6th ed.). Cambridge University Press.
ISBN: 978-0-521-45868-9. DOI:
[doi:10.1017/CBO9780511624612](https://doi.org/10.1017/CBO9780511624612)

Saffman, P. G. (1992). *Vortex dynamics*. Cambridge University Press.
ISBN: 978-0-521-42058-7. DOI:
[doi:10.1017/CBO9780511624063](https://doi.org/10.1017/CBO9780511624063)

## See also

[`two_vortex_ic`](https://robustecologies.github.io/navieRstokes/reference/two_vortex_ic.md),
[`vortex_ic`](https://robustecologies.github.io/navieRstokes/reference/vortex_ic.md),
[`shear_layer_ic`](https://robustecologies.github.io/navieRstokes/reference/shear_layer_ic.md),
[`rotating_ic`](https://robustecologies.github.io/navieRstokes/reference/rotating_ic.md),
[`simulate_navier_stokes`](https://robustecologies.github.io/navieRstokes/reference/simulate_navier_stokes.md)

## Examples

``` r
if (FALSE) { # \dontrun{
result <- simulate_navier_stokes(
  nx = 128, ny = 128, lx = 1.0, ly = 1.0,
  dt = 0.001, nt = 2000, nu = 0.001,
  initial_condition = single_vortex_ic,
  boundary_condition = list(type = "periodic")
)
} # }
```
