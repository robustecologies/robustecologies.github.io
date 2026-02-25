# Plot method for basin_entropy_result

Creates a visualization of the local entropy distribution using ggplot2.

## Usage

``` r
# S3 method for class 'basin_entropy_result'
plot(x, basins = NULL, show_boundary = FALSE, ...)
```

## Arguments

- x:

  A `basin_entropy_result` object.

- basins:

  Optional `wada_basins` object to overlay entropy on basins.

- show_boundary:

  Logical. If TRUE, highlights boundary boxes.

- ...:

  Additional arguments (ignored).

## Value

A ggplot2 object.
