# Extract 2D slice from 3D basins

Extracts a 2D cross-section from 3D basin data at a specified coordinate
along the x, y, or z axis, returning a standard wada_basins object.

## Usage

``` r
slice_3d_basins(x, x_slice = NULL, y_slice = NULL, z_slice = NULL)
```

## Arguments

- x:

  A basin_result_3d object.

- x_slice, y_slice, z_slice:

  Numeric. Slice coordinate (specify exactly one).

## Value

A basin_result object (2D) that can be used with other wadaR functions.

## Examples

``` r
if (FALSE) { # \dontrun{
# First compute 3D basins for the Lorenz system
basins_3d <- compute_basins_3d(
    cpp_dynamics = '
        deriv[0] = sigma * (state[1] - state[0]);
        deriv[1] = state[0] * (rho - state[2]) - state[1];
        deriv[2] = state[0] * state[1] - beta_l * state[2];
    ',
    params = list(sigma = 10, rho = 28, beta_l = 8/3),
    dim = 3,
    attractors = list(
        attractor_point(c(sqrt(72), sqrt(72), 27), 5, "Right"),
        attractor_point(c(-sqrt(72), -sqrt(72), 27), 5, "Left"),
        attractor_exit(0, "Escape")
    ),
    x_range = c(-20, 20),
    y_range = c(-30, 30),
    z_range = c(0, 50),
    resolution = 30,
    t_max = 30
)

# Extract z = 27 slice (Lorenz system at attractor height)
slice_z27 <- slice_3d_basins(basins_3d, z_slice = 27)
plot(slice_z27)

# Extract x = 0 slice
slice_x0 <- slice_3d_basins(basins_3d, x_slice = 0)
plot(slice_x0)

# Compute entropy of slice
entropy <- basin_entropy(slice_x0)
print(entropy)
} # }
```
