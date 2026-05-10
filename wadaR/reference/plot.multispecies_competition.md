# Plot method for multispecies competition objects

Displays the consumption matrix or resource supply vector as heatmaps,
providing an overview of the ecological system configuration.

## Usage

``` r
# S3 method for class 'multispecies_competition'
plot(x, type = c("matrix", "resources"), caption = TRUE, ...)
```

## Arguments

- x:

  A multispecies_competition object.

- type:

  Character. Type of plot: "matrix" for consumption matrix, "resources"
  for resource supply vector. Default is "matrix".

- caption:

  Logical. If TRUE (default), render a one-line caption with the
  function name and primary citation. Set to FALSE to suppress the
  caption when the figure is composed with other panels.

- ...:

  Additional arguments (ignored).

## Value

A ggplot2 object, invisibly.
