# Plot a von Kármán embedded attractor

Visualizes the reconstructed attractor from peak embedding in 2D or 3D.

## Usage

``` r
# S3 method for class 'vonkarman_attractor'
plot(x, type = c("3d", "2d"), ...)
```

## Arguments

- x:

  Object of class "vonkarman_attractor" from
  [`vonkarman_attractor`](https://robustecologies.github.io/tuRbulence/reference/vonkarman_attractor.md).

- type:

  Type of plot: "3d" (default, interactive plotly) or "2d" (ggplot2
  projection).

- ...:

  Additional arguments passed to internal functions.

## Value

Returns the plot object invisibly (plotly for 3D, ggplot2 for 2D).

## Details

The 3D visualization uses plotly for interactive exploration of the
attractor structure. Points are colored by sequence index to show the
flow direction on the attractor.

## See also

[`vonkarman_attractor()`](https://robustecologies.github.io/tuRbulence/reference/vonkarman_attractor.md),
[`print.vonkarman_attractor()`](https://robustecologies.github.io/tuRbulence/reference/print.vonkarman_attractor.md),
[`summary.vonkarman_attractor()`](https://robustecologies.github.io/tuRbulence/reference/summary.vonkarman_attractor.md),
[`vonkarman_sim()`](https://robustecologies.github.io/tuRbulence/reference/vonkarman_sim.md),
[`create_attractor_3d()`](https://robustecologies.github.io/tuRbulence/reference/create_attractor_3d.md)

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- vonkarman_sim(mu = 0.3, n_steps = 100000, seed = 42)
peaks <- vonkarman_peaks(sim)
attr <- vonkarman_attractor(peaks, embed_dim = 3)

# 3D interactive plot
plot(attr)

# 2D projection
plot(attr, type = "2d")
} # }
```
