# PPC fan chart ribbon

Convenience wrapper for PPC ribbon (fan chart) plots.

## Usage

``` r
ppc_ribbon(x, forecast_start = NULL, ...)
```

## Arguments

- x:

  A PPC object.

- forecast_start:

  Optional time index marking forecast start.

- ...:

  Additional arguments passed to `.plot_ppc`.

## Value

Invisibly returns the plot object.

## Details

Implementation of `ppc_ribbon`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## See also

[`plot.demonoid.ppc`](https://robustecologies.github.io/lucifer/reference/plot.demonoid.ppc.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving ppc_ribbon
} # }
```
