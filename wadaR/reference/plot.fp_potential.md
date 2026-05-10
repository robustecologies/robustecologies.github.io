# Plot method for fp_potential objects

Visualizes the drift potential landscape with wells and barriers marked.
When the quasi-potential is available, it can be shown alongside or
overlaid for comparison. The `"landscape3d"` type renders the
quasi-potential as an interactive 3D surface using plotly, with the
stationary density mapped to color and markers at potential minima.

## Usage

``` r
# S3 method for class 'fp_potential'
plot(x, type = c("potential", "quasi", "both", "landscape3d"), ...)
```

## Arguments

- x:

  An `fp_potential` object.

- type:

  Character string selecting the plot type: `"potential"` (default),
  `"quasi"`, `"both"`, or `"landscape3d"`.

- ...:

  Additional arguments (ignored).

## Value

A ggplot2 or plotly object.
