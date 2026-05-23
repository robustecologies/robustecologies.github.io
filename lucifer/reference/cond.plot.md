# Conditional plots

This function provides several styles of conditional plots.

## Usage

``` r
cond.plot(x, y, z, Style = "smoothscatter", col = NULL)
```

## Arguments

- x:

  This required argument accepts a numeric vector.

- y:

  This argument accepts a numeric vector, and is only used with some
  styles.

- z:

  This required argument accepts a discrete vector.

- Style:

  This argument specifies the style of plot, and accepts "boxplot",
  "densover" (density overlay), "hist", "scatter", or "smoothscatter".

- col:

  Optional character vector of colors. When `NULL` (default), the RElab
  palette is used.

## Value

Conditional plots are returned.

## Details

The `cond.plot` function provides simple conditional plots. All plot
styles are conditional upon `z`. Up to nine conditional plots are
produced in a panel.

Plots include: boxplot (y ~ x \| z), densover (f(x \| z)), hist (x \|
z), scatter (x, y \| z), and smoothscatter (x, y \| z).

## See also

[`joint.density.plot`](https://robustecologies.github.io/lucifer/reference/joint.density.plot.md)
and
[`joint.pr.plot`](https://robustecologies.github.io/lucifer/reference/joint.pr.plot.md).

## Examples

``` r
if (FALSE) { # \dontrun{
library(lucifer)
x <- rnorm(1000)
y <- runif(1000)
z <- rcat(1000, rep(1/4,4))
cond.plot(x, y, z, Style="smoothscatter")
} # }
```
