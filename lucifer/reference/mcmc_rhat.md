# MCMC R-hat diagnostic plot

Convenience wrapper for `plot(x, type = "rhat", ...)`.

## Usage

``` r
mcmc_rhat(x, ...)
```

## Arguments

- x:

  An object of class `demonoid`.

- ...:

  Additional arguments passed to `plot.demonoid`.

## Value

Invisibly returns the plot object(s).

## Details

Implementation of `mcmc_rhat`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## See also

[`plot.demonoid`](https://robustecologies.github.io/lucifer/reference/plot.demonoid.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving mcmc_rhat
} # }
```
