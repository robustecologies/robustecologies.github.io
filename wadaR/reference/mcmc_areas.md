# MCMC ridgeline density areas

Convenience wrapper for `plot(x, type = "areas", ...)`.

## Usage

``` r
mcmc_areas(x, Parms = NULL, ground_truth = NULL, ...)
```

## Arguments

- x:

  An object of class `demonoid`.

- Parms:

  Parameter selection (see
  [`plot.demonoid`](https://robustecologies.github.io/lucifer/reference/plot.demonoid.md)).

- ground_truth:

  Named numeric vector of true values.

- ...:

  Additional arguments passed to `plot.demonoid`.

## Value

Invisibly returns the plot object(s).

## Details

Implementation of `mcmc_areas`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## See also

[`plot.demonoid`](https://robustecologies.github.io/lucifer/reference/plot.demonoid.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving mcmc_areas
} # }
```
