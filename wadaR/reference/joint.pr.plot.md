# Joint probability region plot

Given two vectors, this function creates a scatterplot with ellipses
representing probability regions at specified quantiles.

## Usage

``` r
joint.pr.plot(x, y, quantiles = c(0.25, 0.5, 0.75, 0.95), col = NULL)
```

## Arguments

- x:

  a numeric vector, or a two-column matrix (in which case `y` is
  extracted from the second column).

- y:

  a numeric vector.

- quantiles:

  numeric vector of quantiles for which probability region ellipses are
  estimated. Defaults to `c(0.25, 0.50, 0.75, 0.95)`.

- col:

  Optional character vector of colors. When `NULL` (default), the RElab
  palette is used.

## Value

The ggplot2 plot object is returned invisibly.

## Details

A probability region is also commonly called a credible region (see
[`p.interval`](https://robustecologies.github.io/lucifer/reference/p.interval.md)).
Joint probability regions are plotted for two variables, and the
ellipses assume bivariate normality. The center of the ellipse is
plotted by default. This function is often used to plot posterior
samples from
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
and
[`p.interval`](https://robustecologies.github.io/lucifer/reference/p.interval.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(lucifer)
x <- rnorm(100)
y <- rnorm(100)
joint.pr.plot(x, y)
} # }
```
