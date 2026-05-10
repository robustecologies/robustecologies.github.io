# Plot 3D attractor trajectory as line

Create an interactive 3D line plot showing the trajectory through phase
space. Unlike
[`plot.takens_embedding`](https://robustecologies.github.io/chamaeleon/reference/plot.takens_embedding.md)
which shows points, this function draws a continuous line connecting
successive states, making the temporal evolution of the trajectory
visible.

## Usage

``` r
plot_trajectory_3d(
  x,
  color_by = NULL,
  color_label = "Time",
  line_width = 1,
  axes = 1:3,
  main = NULL,
  ...
)
```

## Arguments

- x:

  Object of class `takens_embedding` or numeric matrix with at least 3
  columns representing the embedded trajectory.

- color_by:

  Optional numeric vector for coloring the trajectory. If NULL, defaults
  to time index.

- color_label:

  Character. Label for the colorbar (default "Time").

- line_width:

  Numeric. Width of the trajectory line (default 1).

- axes:

  Integer vector of length 3. Which columns to use as x, y, z
  coordinates. Default is `1:3`.

- main:

  Character. Plot title. If NULL, a default title is generated.

- ...:

  Additional arguments passed to plotly.

## Value

A plotly htmlwidget.

## Details

This visualization is particularly useful for:

- Understanding the temporal structure of the attractor

- Identifying transient behavior or regime transitions

- Visualizing how trajectories move through different regions

## References

Takens F (1981). Detecting strange attractors in turbulence. In: Lecture
Notes in Mathematics, 898. Springer.
[doi:10.1007/BFb0091924](https://doi.org/10.1007/BFb0091924)

Hegger R, Kantz H, Schreiber T (1999). Practical implementation of
nonlinear time series methods: the TISEAN package. Chaos 9:413-435.
[doi:10.1063/1.166424](https://doi.org/10.1063/1.166424)

## See also

[`plot.takens_embedding`](https://robustecologies.github.io/chamaeleon/reference/plot.takens_embedding.md)
for point-based attractor plots;
[`takens_embed`](https://robustecologies.github.io/chamaeleon/reference/takens_embed.md)
for the underlying embedding construction.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(42)
x <- sin(seq(0, 50, by = 0.01)) + 0.5 * sin(2.3 * seq(0, 50, by = 0.01))
embedded <- takens_embed(x, m = 3, tau = 20)
plot_trajectory_3d(embedded)
} # }
```
