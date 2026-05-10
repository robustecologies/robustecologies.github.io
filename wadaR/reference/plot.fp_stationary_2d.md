# Plot method for fp_stationary_2d objects

Visualizes the 2D Fokker-Planck stationary density as a filled contour
plot or an interactive 3D surface.

## Usage

``` r
# S3 method for class 'fp_stationary_2d'
plot(x, type = c("density", "surface", "both"), ...)
```

## Arguments

- x:

  An `fp_stationary_2d` object.

- type:

  Character string: `"density"` (default filled contour), `"surface"`
  (3D plotly surface), or `"both"` (side by side, renders the contour
  and prints the surface).

- ...:

  Additional arguments (ignored).

## Value

A ggplot2 or plotly object.
