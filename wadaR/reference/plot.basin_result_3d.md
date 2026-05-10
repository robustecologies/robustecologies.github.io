# Plot method for 3D basin results

Creates an interactive 3D plotly scatter plot of basin assignments,
colored by attractor membership.

## Usage

``` r
# S3 method for class 'basin_result_3d'
plot(x, slice = NULL, slice_index = NULL, caption = TRUE, ...)
```

## Arguments

- x:

  A basin_result_3d object.

- slice:

  Character. Axis to slice along: "x", "y", or "z". If NULL, plots a 3D
  scatter of all boundary points.

- slice_index:

  Integer. Index of the slice along the chosen axis.

- caption:

  Logical. If TRUE (default), render a one-line caption with the
  function name on the 2D slice plot. The 3D plotly branch ignores this
  argument because plotly widgets carry no caption layer.

- ...:

  Additional arguments (ignored).

## Value

A plotly widget (3D) or ggplot2 object (2D slice).

## References

Daza, A., Wagemakers, A., Georgeot, B., Guery-Odelin, D., & Sanjuan, M.
A. F. (2016). Basin entropy: a new tool to analyze uncertainty in
dynamical systems. *Scientific Reports*, 6, 31416.
[doi:10.1038/srep31416](https://doi.org/10.1038/srep31416)
