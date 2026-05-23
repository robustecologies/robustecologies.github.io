# Spatially localized forcing (vortex generation)

Applies a localized rotational force to generate vorticity at a specific
location in the domain.

## Usage

``` r
localized_vortex_force(x, y, t, x0 = 0.5, y0 = 0.5, sigma = 0.1, A = 1)
```

## Arguments

- x:

  Numeric. x-coordinate

- y:

  Numeric. y-coordinate

- t:

  Numeric. Time (not used)

- x0:

  Numeric. Center x-coordinate (default = 0.5)

- y0:

  Numeric. Center y-coordinate (default = 0.5)

- sigma:

  Numeric. Forcing radius (default = 0.1)

- A:

  Numeric. Force amplitude (default = 1.0)

## Value

List with components fx and fy

## Details

The forcing is constructed as a rotational field with Gaussian spatial
localization. Defining \\r^2 = (x - x_0)^2 + (y - y_0)^2\\, the forcing
components are: \$\$f_x = -A(y - y_0)
\exp\left(-\frac{r^2}{\sigma^2}\right)\$\$ \$\$f_y = A(x - x_0)
\exp\left(-\frac{r^2}{\sigma^2}\right)\$\$ This forcing injects positive
(counter-clockwise) vorticity into the flow, with the curl of the
forcing field given by: \$\$\nabla \times \mathbf{f} = 2A\left(1 -
\frac{r^2}{\sigma^2}\right)\exp\left(-\frac{r^2}{\sigma^2}\right)\$\$
The parameter \\\sigma\\ controls the spatial extent of the forcing
region, while \\A\\ controls the strength of vorticity injection. This
type of forcing is useful for studying vortex dynamics, wake formation,
and turbulence generation.

## References

Saffman, P. G. (1992). *Vortex dynamics*. Cambridge University Press.
ISBN: 978-0-521-42058-7. DOI:
[doi:10.1017/CBO9780511624063](https://doi.org/10.1017/CBO9780511624063)

Kundu, P. K., Cohen, I. M., & Dowling, D. R. (2015). *Fluid mechanics*
(6th ed.). Academic Press. ISBN: 978-0-12-405935-1.

## See also

[`pressure_gradient_force`](https://robustecologies.github.io/navieRstokes/reference/pressure_gradient_force.md),
[`oscillatory_force`](https://robustecologies.github.io/navieRstokes/reference/oscillatory_force.md),
[`taylor_green_force`](https://robustecologies.github.io/navieRstokes/reference/taylor_green_force.md),
[`simulate_navier_stokes`](https://robustecologies.github.io/navieRstokes/reference/simulate_navier_stokes.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Note: localized_vortex_force adds circulation and may cause
# numerical instabilities. Use with caution and conservative parameters.
result <- simulate_navier_stokes(
  nx = 64, ny = 64, lx = 1.0, ly = 1.0,
  dt = 0.001, nt = 1000, nu = 0.01,
  initial_condition = function(x, y) list(u = 0, v = 0),
  forcing_function = localized_vortex_force,
  boundary_condition = list(type = "periodic")
)
} # }
```
