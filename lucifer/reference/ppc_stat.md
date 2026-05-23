# PPC test statistic

Convenience wrapper for PPC test statistic plots.

## Usage

``` r
ppc_stat(x, stat_fun = mean, ...)
```

## Arguments

- x:

  A PPC object.

- stat_fun:

  A function computing a scalar test statistic (default: `mean`).

- ...:

  Additional arguments passed to `.plot_ppc`.

## Value

Invisibly returns the plot object.

## Details

Implementation of `ppc_stat`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## See also

[`plot.demonoid.ppc`](https://robustecologies.github.io/lucifer/reference/plot.demonoid.ppc.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving ppc_stat
} # }
```
