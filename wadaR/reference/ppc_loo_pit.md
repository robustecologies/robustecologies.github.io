# PPC LOO probability integral transform

Convenience wrapper for LOO-PIT ECDF plots.

## Usage

``` r
ppc_loo_pit(x, loo, ...)
```

## Arguments

- x:

  A PPC object.

- loo:

  A LOO object from
  [`LOO`](https://robustecologies.github.io/lucifer/reference/LOO.md).

- ...:

  Additional arguments passed to `.plot_ppc`.

## Value

Invisibly returns the plot object.

## Details

Implementation of `ppc_loo_pit`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## See also

[`plot.demonoid.ppc`](https://robustecologies.github.io/lucifer/reference/plot.demonoid.ppc.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving ppc_loo_pit
} # }
```
