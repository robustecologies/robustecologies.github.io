# Plot method for bifurcation_sweep objects

Plot method for bifurcation_sweep objects

## Usage

``` r
# S3 method for class 'bifurcation_sweep'
plot(x, state = NULL, title = NULL, ylim = NULL, ...)
```

## Arguments

- x:

  A `bifurcation_sweep` object.

- state:

  Character; which state variable to plot on the y-axis. If `NULL`
  (default), plots all states in facets.

- title:

  Optional plot title.

- ylim:

  Optional numeric vector of length 2 restricting the y-axis range.
  Useful for zooming into biologically relevant region.

- ...:

  Additional arguments (ignored).

## Value

A ggplot2 object.
