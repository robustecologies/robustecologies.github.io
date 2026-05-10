# Plot method for bifurcation_result

Creates visualizations of bifurcation analysis results showing how basin
properties change with the parameter.

## Usage

``` r
# S3 method for class 'bifurcation_result'
plot(
  x,
  type = c("entropy", "boundary", "wada", "basins", "combined"),
  indices = NULL,
  palette = "viridis",
  caption = TRUE,
  ...
)
```

## Arguments

- x:

  A bifurcation_result object.

- type:

  Character. Type of plot:

  - "entropy": Basin entropy vs parameter (default)

  - "boundary": Boundary fraction vs parameter

  - "wada": Wada statistic vs parameter

  - "basins": Basin plots at selected parameter values

  - "combined": Multi-panel summary plot

- indices:

  Integer vector. For type="basins", which parameter indices to plot.

- palette:

  Character. Color palette name.

- caption:

  Logical. If TRUE (default), render a one-line caption with the
  function name and primary citation. For `type = "basins"` and
  `type = "combined"`, child panels never carry their own caption to
  avoid duplication under composition; only the outer patchwork
  annotation honours `caption`.

- ...:

  Additional arguments passed to plotting functions.

## Value

A ggplot2 object.

## Examples

``` r
if (FALSE) { # \dontrun{
# Plot entropy evolution
plot(bif_result, type = "entropy")

# Plot selected basin snapshots
plot(bif_result, type = "basins", indices = c(1, 10, 20))

# Combined summary plot
plot(bif_result, type = "combined")
} # }
```
