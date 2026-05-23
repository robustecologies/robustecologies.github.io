# Compute vorticity from velocity field

Computes the vorticity field from velocity components using second-order
central differences: \\\omega = \frac{\partial v}{\partial x} -
\frac{\partial u}{\partial y}\\

This is an R wrapper for the optimized C++ implementation.

## Usage

``` r
compute_vorticity(u, v, dx, dy)
```

## Arguments

- u:

  Numeric matrix. Velocity component in x-direction

- v:

  Numeric matrix. Velocity component in y-direction

- dx:

  Numeric. Grid spacing in x-direction

- dy:

  Numeric. Grid spacing in y-direction

## Value

Numeric matrix containing the vorticity field

## Details

The vorticity is a measure of local rotation in the fluid. Positive
values indicate counter-clockwise rotation, negative values indicate
clockwise rotation.

The calculation uses second-order central differences: \$\$\omega\_{i,j}
= \frac{v\_{i+1,j} - v\_{i-1,j}}{2\Delta x} - \frac{u\_{i,j+1} -
u\_{i,j-1}}{2\Delta y}\$\$

## References

Kundu, P. K., Cohen, I. M., & Dowling, D. R. (2015). *Fluid mechanics*
(6th ed.). Academic Press. ISBN: 978-0-12-405935-1.

Ferziger, J. H., Peric, M., & Street, R. L. (2020). *Computational
methods for fluid dynamics* (4th ed.). Springer. ISBN:
978-3-319-99691-2. DOI:
[doi:10.1007/978-3-319-99693-6](https://doi.org/10.1007/978-3-319-99693-6)

## See also

[`simulate_navier_stokes`](https://robustecologies.github.io/navieRstokes/reference/simulate_navier_stokes.md),
[`plot.navieRstokes`](https://robustecologies.github.io/navieRstokes/reference/plot.navieRstokes.md),
[`flow.navieRstokes`](https://robustecologies.github.io/navieRstokes/reference/flow.md),
[`summary.navieRstokes`](https://robustecologies.github.io/navieRstokes/reference/summary.navieRstokes.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# After running a simulation
result <- simulate_navier_stokes(
  nx = 64, ny = 64,
  lx = 1.0, ly = 1.0,
  dt = 0.001, nt = 1000,
  nu = 0.01,
  initial_condition = function(x, y) list(u = 0, v = 0),
  boundary_condition = list(
    type = "dirichlet",
    values = list(
      u_left = 0, u_right = 0, u_top = 1, u_bottom = 0,
      v_left = 0, v_right = 0, v_top = 0, v_bottom = 0
    )
  )
)

# Compute vorticity at final time
omega <- compute_vorticity(
  result$u[, , dim(result$u)[3]],
  result$v[, , dim(result$v)[3]],
  result$parameters$dx,
  result$parameters$dy
)

# Visualize
image(result$x, result$y, omega,
  col = hcl.colors(50, "RdBu", rev = TRUE),
  xlab = "x", ylab = "y",
  main = "Vorticity field"
)
contour(result$x, result$y, omega, add = TRUE)
} # }
```
