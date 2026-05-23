# Two counter-rotating vortices initial condition

Creates two Gaussian vortices with opposite circulation that interact
and merge over time.

## Usage

``` r
two_vortex_ic(x, y, Gamma = 1, sigma = 0.1, lx = 1, ly = 1)
```

## Arguments

- x:

  Numeric. x-coordinate

- y:

  Numeric. y-coordinate

- Gamma:

  Numeric. Circulation strength (default = 1.0)

- sigma:

  Numeric. Vortex core radius (default = 0.1)

- lx:

  Numeric. Domain length in x-direction (default = 1.0)

- ly:

  Numeric. Domain length in y-direction (default = 1.0)

## Value

List with components u and v

## Details

This initial condition creates a pair of Lamb-Oseen vortices with
opposite circulation signs, positioned at \\(0.3 l_x, 0.5 l_y)\\ and
\\(0.7 l_x, 0.5 l_y)\\. Each vortex has the smooth Gaussian velocity
profile: \$\$u\_\theta(r) =
\frac{\Gamma}{2\pi\sigma}\frac{r}{\sigma}\exp\left(-\frac{r^2}{2\sigma^2}\right)\$\$
where \\r\\ is the distance from the vortex center and \\\Gamma\\ is the
circulation. The velocity field is converted to Cartesian components:
\$\$u = -u\_\theta \frac{(y - y_c)}{r}, \quad v = u\_\theta \frac{(x -
x_c)}{r}\$\$ This profile is smooth at the vortex center (no
singularity) and decays exponentially with distance. The two vortices
induce mutual advection, causing them to orbit around their common
centroid while viscosity causes core expansion and eventual merger.

## References

Lamb, H. (1932). *Hydrodynamics* (6th ed.). Cambridge University Press.
ISBN: 978-0-521-45868-9. DOI:
[doi:10.1017/CBO9780511624612](https://doi.org/10.1017/CBO9780511624612)

Saffman, P. G. (1992). *Vortex dynamics*. Cambridge University Press.
ISBN: 978-0-521-42058-7. DOI:
[doi:10.1017/CBO9780511624063](https://doi.org/10.1017/CBO9780511624063)

## See also

[`single_vortex_ic`](https://robustecologies.github.io/navieRstokes/reference/single_vortex_ic.md),
[`vortex_ic`](https://robustecologies.github.io/navieRstokes/reference/vortex_ic.md),
[`shear_layer_ic`](https://robustecologies.github.io/navieRstokes/reference/shear_layer_ic.md),
[`rotating_ic`](https://robustecologies.github.io/navieRstokes/reference/rotating_ic.md),
[`simulate_navier_stokes`](https://robustecologies.github.io/navieRstokes/reference/simulate_navier_stokes.md)

## Examples

``` r
if (FALSE) { # \dontrun{
result <- simulate_navier_stokes(
  nx = 128, ny = 128, lx = 1.0, ly = 1.0,
  dt = 0.001, nt = 1000, nu = 0.001,
  initial_condition = two_vortex_ic,
  boundary_condition = list(type = "periodic")
)
} # }
```
