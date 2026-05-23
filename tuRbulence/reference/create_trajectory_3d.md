# Create static 3D phase space trajectory

Shows the full (x, y, z) phase space as a point cloud colored by time.
Since von Kármán is inherently stochastic, points are used rather than
lines for cleaner visualization.

## Usage

``` r
create_trajectory_3d(
  sim,
  n_points = 10000,
  color_by = c("time", "velocity", "z"),
  point_size = 1.5
)
```

## Arguments

- sim:

  Object of class "vonkarman_sim".

- n_points:

  Maximum points to display.

- color_by:

  Variable to map to color: "time", "velocity", or "z".

- point_size:

  Size of scatter points.

## Value

A plotly object.

## See also

[`create_attractor_3d()`](https://robustecologies.github.io/tuRbulence/reference/create_attractor_3d.md),
[`create_animated_attractor()`](https://robustecologies.github.io/tuRbulence/reference/create_animated_attractor.md),
[`create_animated_attractor_accumulate()`](https://robustecologies.github.io/tuRbulence/reference/create_animated_attractor_accumulate.md),
[`save_attractor_html()`](https://robustecologies.github.io/tuRbulence/reference/save_attractor_html.md),
[`export_animation_frames()`](https://robustecologies.github.io/tuRbulence/reference/export_animation_frames.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Create phase space trajectory
sim <- vonkarman_sim(mu = 0.3, n_steps = 50000, seed = 42)
p <- create_trajectory_3d(sim, n_points = 5000, color_by = "time")
p
} # }
```
