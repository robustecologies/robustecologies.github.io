# Compute 3D basins of attraction

Computes basins of attraction on a 3D grid for higher-dimensional
systems. Uses C++/OpenMP parallelization for high performance with
interrupt support.

## Usage

``` r
compute_basins_3d(
  cpp_dynamics,
  params,
  attractors,
  dim = 3,
  type = c("ode", "map"),
  x_range,
  y_range,
  z_range,
  resolution = 50,
  t_max = 100,
  dt = 0.01,
  n_cores = NULL,
  z_var = 3,
  fixed_vars = NULL,
  verbose = TRUE
)
```

## Arguments

- cpp_dynamics:

  Character string containing C++ code for the dynamics. Same format as
  [`compile_basin_function`](https://robustecologies.github.io/wadaR/reference/compile_basin_function.md).

- params:

  Named list of parameters.

- attractors:

  List of attractor specifications.

- dim:

  Integer. State space dimension (must be \>= 3).

- type:

  Character. System type: "ode" or "map".

- x_range, y_range, z_range:

  Numeric vectors of length 2. Ranges for each axis.

- resolution:

  Integer or vector. Grid resolution(s).

- t_max:

  Numeric. Maximum integration time.

- dt:

  Numeric. Integration time step.

- n_cores:

  Integer. Number of CPU cores to use for parallel computation. Default
  is `NULL`, which uses `parallel::detectCores(logical = FALSE) - 1`.

- z_var:

  Integer. Which state variable maps to z-axis (default 3).

- fixed_vars:

  Named list. Fixed values for dimensions \> 3.

- verbose:

  Logical. Show progress.

## Value

A list of class `basin_result_3d` containing:

- basin_array:

  3D array of basin assignments.

- x_grid,y_grid,z_grid:

  Grid coordinates.

- n_attractors:

  Number of attractors.

- attractors:

  Attractor specifications.

## Details

For systems with dimension \> 3, a 3D projection is computed by:

1.  Selecting 3 state variables for x, y, z axes

2.  Fixing remaining variables at specified values

The computation uses OpenMP parallelization across the 3D grid points.

## References

- Nusse, H. E., & Yorke, J. A. (1998). *Dynamics: Numerical
  explorations*. Springer-Verlag. ISBN: 978-0-387-98267-3

- Kennedy, J., & Yorke, J. A. (1991). Basins of Wada. *Physica D*,
  51(1-3), 213-225.
  [doi:10.1016/0167-2789(91)90234-Z](https://doi.org/10.1016/0167-2789%2891%2990234-Z)

## See also

[`plot_3d_basins`](https://robustecologies.github.io/wadaR/reference/plot_3d_basins.md),
[`slice_3d_basins`](https://robustecologies.github.io/wadaR/reference/slice_3d_basins.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# ============================================================ #
# Example: Lorenz system 3D basins                             #
# ============================================================ #
# The Lorenz system: dx/dt = sigma*(y-x), dy/dt = x*(rho-z)-y,
#                    dz/dt = x*y - beta*z

basins_3d <- compute_basins_3d(
    cpp_dynamics = '
        deriv[0] = sigma * (state[1] - state[0]);
        deriv[1] = state[0] * (rho - state[2]) - state[1];
        deriv[2] = state[0] * state[1] - beta_l * state[2];
    ',
    params = list(sigma = 10, rho = 28, beta_l = 8/3),
    dim = 3,
    attractors = list(
        attractor_point(c(sqrt(72), sqrt(72), 27), 3, "Right wing"),
        attractor_point(c(-sqrt(72), -sqrt(72), 27), 3, "Left wing"),
        attractor_exit(0, "Escape")
    ),
    x_range = c(-20, 20),
    y_range = c(-30, 30),
    z_range = c(0, 50),
    resolution = 50,
    t_max = 50
)

# Interactive 3D plot
plot_3d_basins(basins_3d)

# 2D slice at z = 27
slice <- slice_3d_basins(basins_3d, z = 27)
plot(slice)
} # }
```
