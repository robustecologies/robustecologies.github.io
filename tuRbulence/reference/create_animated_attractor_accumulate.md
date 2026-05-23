# Create animated attractor showing progressive formation

Shows the attractor being built up point by point with connecting lines.
This clearly shows the structure forming over time.

## Usage

``` r
create_animated_attractor_accumulate(
  sim,
  n_frames = 50,
  min_separation = 0.1,
  smooth_window = 30,
  point_size = 3,
  line_width = 1.5,
  fps = 10,
  colorscale = "Plasma",
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

- min_separation:

  Minimum separation for peak extraction.

- smooth_window:

  Smoothing window size.

- point_size:

  Size of scatter points.

- line_width:

  Width of connecting lines.

- fps:

  Frames per second for animation.

- colorscale:

  Plotly colorscale for points (default "Plasma").

- camera:

  List specifying camera position with eye coordinates.

- axis_margin:

  Multiplier for axis range padding (default 1.1).

- verbose:

  Logical; if TRUE, print messages when parameters are adjusted.

## Value

A plotly object with animation controls.

## See also

[`create_animated_attractor()`](https://robustecologies.github.io/tuRbulence/reference/create_animated_attractor.md),
[`create_attractor_3d()`](https://robustecologies.github.io/tuRbulence/reference/create_attractor_3d.md),
[`create_trajectory_3d()`](https://robustecologies.github.io/tuRbulence/reference/create_trajectory_3d.md),
[`save_attractor_html()`](https://robustecologies.github.io/tuRbulence/reference/save_attractor_html.md),
[`export_animation_frames()`](https://robustecologies.github.io/tuRbulence/reference/export_animation_frames.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate chaotic attractor
sim <- vonkarman_sim(mu = 0.3, n_steps = 100000, seed = 42)

# Create accumulating animation
anim <- create_animated_attractor_accumulate(sim, n_frames = 50)
anim
} # }
```
