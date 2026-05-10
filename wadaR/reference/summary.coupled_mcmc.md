# Summary method for coupled MCMC objects

Produces an extended summary including unbiased posterior mean estimates
via the debiasing formula, naive chain-X posterior summaries (discarding
the first `k` iterations as burn-in), and meeting diagnostics.

## Usage

``` r
# S3 method for class 'coupled_mcmc'
summary(object, h = NULL, ...)
```

## Arguments

- object:

  An object of class `coupled_mcmc`.

- h:

  Optional function applied to each parameter vector before computing
  the debiased estimator. Defaults to the identity (posterior means).

- ...:

  Additional arguments (currently ignored).

## Value

Invisibly returns `object`.

## Details

Implementation of `summary.coupled_mcmc`. Refer to the package vignettes
and the cited references for a complete algorithmic and mathematical
description.

## See also

[`coupled_mcmc`](https://robustecologies.github.io/lucifer/reference/coupled_mcmc.md),
[`plot.coupled_mcmc`](https://robustecologies.github.io/lucifer/reference/plot.coupled_mcmc.md),
[`print.coupled_mcmc`](https://robustecologies.github.io/lucifer/reference/print.coupled_mcmc.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving summary.coupled_mcmc
} # }
```
