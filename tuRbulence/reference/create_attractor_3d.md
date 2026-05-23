# Create static 3D embedded attractor

Shows the complete embedded attractor structure using points colored by
sequence to visualize the flow direction. Since von Kármán is inherently
stochastic, points are used rather than lines for cleaner visualization.

## Usage

``` r
create_attractor_3d(
  sim,
  min_separation = 0.05,
  smooth_window = 0,
  point_size = NULL,
  colorscale = "Plasma"
)
```

## Arguments

- sim:

  Object of class "vonkarman_sim".

- min_separation:

  Minimum time separation between extracted peaks. Smaller values yield
  more points. Default 0.05.

- smooth_window:

  Smoothing window size. Default 0 (no smoothing).

- point_size:

  Size of scatter points. Auto-adjusted based on point count if not
  specified.

- colorscale:

  Plotly colorscale name ("Viridis", "Plasma", "Inferno", etc).

## Value

A plotly object.

## Details

This function uses peak embedding as described in the original paper.
For a direct phase space view with more points, use
[`create_trajectory_3d`](https://robustecologies.github.io/tuRbulence/reference/create_trajectory_3d.md).

## See also

[`create_trajectory_3d()`](https://robustecologies.github.io/tuRbulence/reference/create_trajectory_3d.md),
[`create_animated_attractor()`](https://robustecologies.github.io/tuRbulence/reference/create_animated_attractor.md),
[`create_animated_attractor_accumulate()`](https://robustecologies.github.io/tuRbulence/reference/create_animated_attractor_accumulate.md),
[`save_attractor_html()`](https://robustecologies.github.io/tuRbulence/reference/save_attractor_html.md),
[`export_animation_frames()`](https://robustecologies.github.io/tuRbulence/reference/export_animation_frames.md)

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- vonkarman_sim(mu = 0.3, n_steps = 100000, seed = 42)

# Peak-embedded attractor (Poincaré section)
p <- create_attractor_3d(sim, min_separation = 0.05)
p

# For denser visualization, use direct phase space
p2 <- create_trajectory_3d(sim, n_points = 10000)
p2
} # }
```
