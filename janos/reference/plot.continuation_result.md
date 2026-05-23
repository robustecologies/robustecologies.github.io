# Plot method for continuation_result objects

Produces visualizations of a bifurcation continuation result. The
default type shows the bifurcation diagram with the continuation
parameter on the x-axis and state variables on the y-axis, colored by
stability. The eigenvalue type shows how the real parts of the
eigenvalues change along the branch.

## Usage

``` r
# S3 method for class 'continuation_result'
plot(x, type = c("bifurcation", "eigenvalue"), state = NULL, ...)
```

## Arguments

- x:

  A `continuation_result` object

- type:

  Type of plot: "bifurcation" (default) or "eigenvalue"

- state:

  Character vector of state variable names to plot. If `NULL`, plots all
  state variables. For the bifurcation diagram, produces one panel per
  state (faceted).

- ...:

  Additional arguments (ignored)

## Value

A ggplot2 object.
