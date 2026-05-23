# Save plotly attractor as interactive HTML

Save plotly attractor as interactive HTML

## Usage

``` r
save_attractor_html(
  p,
  filename,
  selfcontained = TRUE,
  title = "Von Karman Attractor"
)
```

## Arguments

- p:

  Plotly object from create_animated_attractor or other plot functions.

- filename:

  Output HTML filename.

- selfcontained:

  If TRUE, embed all dependencies in single file.

- title:

  HTML page title.

## See also

[`create_attractor_3d()`](https://robustecologies.github.io/tuRbulence/reference/create_attractor_3d.md),
[`create_trajectory_3d()`](https://robustecologies.github.io/tuRbulence/reference/create_trajectory_3d.md),
[`create_animated_attractor()`](https://robustecologies.github.io/tuRbulence/reference/create_animated_attractor.md),
[`create_animated_attractor_accumulate()`](https://robustecologies.github.io/tuRbulence/reference/create_animated_attractor_accumulate.md),
[`export_animation_frames()`](https://robustecologies.github.io/tuRbulence/reference/export_animation_frames.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Create and save animated attractor
sim <- vonkarman_sim(mu = 0.3, n_steps = 50000)
anim <- create_animated_attractor(sim, n_frames = 30)
save_attractor_html(anim, "attractor_animation.html")
} # }
```
