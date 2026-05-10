# Joint density plot

This function plots the joint kernel density from samples of two
marginal posterior distributions.

## Usage

``` r
joint.density.plot(
  x,
  y,
  Title = NULL,
  contour = TRUE,
  color = FALSE,
  Trace = NULL,
  col = NULL
)
```

## Arguments

- x, y:

  vectors of samples from two marginal posterior distributions, such as
  those output by
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
  in components `Posterior1` or `Posterior2`.

- Title:

  character; the title of the joint density plot.

- contour:

  logical; if `TRUE` (default), contour lines are added to the plot.

- color:

  logical; if `TRUE`, color fill is added to the plot. Defaults to
  `FALSE`.

- Trace:

  a length-2 vector specifying beginning and ending iteration/sample to
  trace the exploration of the joint density, or `NULL` (default) to
  disable tracing.

- col:

  Optional character vector of colors. When `NULL` (default), the RElab
  palette is used.

## Value

The ggplot2 plot object is returned invisibly.

## Details

This function produces either a bivariate scatterplot with kernel
density contour lines, or a bivariate plot with kernel
density-influenced colors. The `Trace` argument allows the user to view
the exploration path of the joint density, such as from MCMC chain
output. An efficient algorithm jumps to random points of the joint
density, while an inefficient algorithm explores more slowly. The
initial point of the trace is plotted with a colored dot.

## See also

[`IAT`](https://robustecologies.github.io/lucifer/reference/IAT.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`PosteriorChecks`](https://robustecologies.github.io/lucifer/reference/PosteriorChecks.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(lucifer)
X <- rmvn(1000, runif(2), diag(2))
joint.density.plot(X[,1], X[,2], Title = "Joint Density Plot",
     contour = TRUE, color = FALSE)
joint.density.plot(X[,1], X[,2], Title = "Joint Density Plot",
     contour = TRUE, color = TRUE)
joint.density.plot(X[,1], X[,2], Title = "Joint Trace Plot",
     contour = FALSE, color = TRUE, Trace = c(1, 10))
} # }
```
