# Create animated 3D attractor with playback controls

Generates an interactive plotly animation showing the attractor
evolution with play/pause controls. Shows both trajectory lines and
points.

## Usage

``` r
create_animated_attractor(
  sim,
  n_frames = 50,
  trail_length = 100,
  min_separation = 0.1,
  smooth_window = 30,
  point_size = 6,
  line_width = 2,
  fps = 10,
  colorscale = NULL,
  camera = list(eye = list(x = 2, y = 0.5, z = 1)),
  axis_margin = 1.1,
  verbose = TRUE
)
```

## Arguments

- sim:

  Object of class "vonkarman_sim".

- n_frames:

  Number of animation frames.

- trail_length:

  Number of points to show as trailing history.

- min_separation:

  Minimum separation for peak extraction.

- smooth_window:

  Smoothing window size.

- point_size:

  Size of the leading point marker.

- line_width:

  Width of trajectory line.

- fps:

  Frames per second for animation.

- colorscale:

  Plotly colorscale for the trajectory gradient (default "Blues").

- camera:

  List specifying camera position with eye coordinates.

- axis_margin:

  Multiplier for axis range padding (default 1.1).

- verbose:

  Logical; if TRUE, print messages when parameters are adjusted.

## Value

A plotly object with animation controls.

## See also

[`create_animated_attractor_accumulate()`](https://robustecologies.github.io/tuRbulence/reference/create_animated_attractor_accumulate.md),
[`create_attractor_3d()`](https://robustecologies.github.io/tuRbulence/reference/create_attractor_3d.md),
[`create_trajectory_3d()`](https://robustecologies.github.io/tuRbulence/reference/create_trajectory_3d.md),
[`save_attractor_html()`](https://robustecologies.github.io/tuRbulence/reference/save_attractor_html.md),
[`export_animation_frames()`](https://robustecologies.github.io/tuRbulence/reference/export_animation_frames.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate chaotic attractor
sim <- vonkarman_sim(mu = 0.3, n_steps = 100000, seed = 42)

# Create animation with trailing trajectory
anim <- create_animated_attractor(sim, n_frames = 50, trail_length = 100)
anim
} # }
```
