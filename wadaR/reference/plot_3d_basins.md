# Plot 3D basins using plotly

Creates an interactive 3D visualization of basins of attraction using
plotly, supporting isosurface, scatter, and volume rendering modes.

## Usage

``` r
plot_3d_basins(
  x,
  type = c("isosurface", "scatter", "volume"),
  opacity = 0.6,
  colors = NULL,
  show_attractors = TRUE,
  subsample = 0.2,
  ...
)
```

## Arguments

- x:

  A basin_result_3d object.

- type:

  Character. Visualization type:

  - "isosurface": 3D isosurfaces showing basin boundaries (default)

  - "scatter": 3D scatter plot of basin points

  - "volume": Volume rendering (requires high memory)

- opacity:

  Numeric. Transparency (0-1).

- colors:

  Character vector. Colors for each basin.

- show_attractors:

  Logical. Mark attractor positions.

- subsample:

  Numeric. Fraction of points to plot for scatter (0-1).

- ...:

  Additional arguments passed to plotly.

## Value

A plotly object.

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

# Interactive 3D visualization (isosurface)
plot_3d_basins(basins_3d)

# Scatter plot with subsampling
plot_3d_basins(basins_3d, type = "scatter", subsample = 0.1)

# Custom colors
plot_3d_basins(basins_3d, colors = c("#E74C3C", "#3498DB", "#2ECC71"))
} # }
```
