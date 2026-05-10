# Create animated rotation of 3D basins

Creates an animated visualization of 3D basins, either as an interactive
HTML plotly rotation or as a GIF of 2D slices sweeping through the
z-axis.

## Usage

``` r
animate_3d_rotation(
  x,
  filename,
  n_frames = 36,
  fps = 10,
  width = 800,
  height = 600,
  ...
)
```

## Arguments

- x:

  A basin_result_3d object.

- filename:

  Character. Output filename (.html for interactive, .gif for animated
  slices).

- n_frames:

  Integer. Number of animation frames (default 36).

- fps:

  Integer. Frames per second (default 10).

- width, height:

  Integer. Animation dimensions in pixels.

- ...:

  Additional arguments.

## Value

Invisibly returns the animation object.

## Details

For HTML output (.html extension), creates an interactive plotly plot
with camera rotation controls. For GIF output (.gif extension), creates
an animated sequence of 2D slices through the z-axis using gganimate.

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

# Create HTML animation (interactive rotation)
animate_3d_rotation(basins_3d, "rotation.html", n_frames = 36)

# Create GIF animation (z-slice sweep)
animate_3d_rotation(basins_3d, "slices.gif", n_frames = 30)
} # }
```
