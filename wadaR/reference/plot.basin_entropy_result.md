# Plot method for basin_entropy_result

Creates a visualization of the local entropy distribution using ggplot2.

## Usage

``` r
# S3 method for class 'basin_entropy_result'
plot(
  x,
  basins = NULL,
  show_boundary = FALSE,
  plotly = FALSE,
  caption = TRUE,
  ...
)
```

## Arguments

- x:

  A `basin_entropy_result` object.

- basins:

  Optional `wada_basins` object to overlay entropy on basins.

- show_boundary:

  Logical. If TRUE, highlights boundary boxes.

- plotly:

  Logical. If TRUE, returns an interactive plotly plot instead of a
  static ggplot2 plot. Default is FALSE.

- caption:

  Logical. If TRUE (default), render a one-line caption with the
  function name and primary citation. Set to FALSE to suppress the
  caption when the figure is composed with other panels.

- ...:

  Additional arguments (ignored).

## Value

A `ggplot2` object (if `plotly = FALSE`) or a `plotly` object (if
`plotly = TRUE`) that can be further customized.
