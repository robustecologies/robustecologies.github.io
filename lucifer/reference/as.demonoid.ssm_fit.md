# Convert ssm_fit to demonoid

Extracts the parameter posterior from an `ssm_fit` object and wraps it
in a `demonoid` structure for compatibility with the lucifer diagnostic
ecosystem (Rhat, ESS, LOO, etc.).

## Usage

``` r
# S3 method for class 'ssm_fit'
as.demonoid(x, ...)
```

## Arguments

- x:

  An object of class `ssm_fit`.

- ...:

  Additional arguments (unused).

## Value

A `demonoid` object.

## Details

Implementation of `as.demonoid.ssm_fit`. Refer to the package vignettes
and the cited references for a complete algorithmic and mathematical
description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving as.demonoid.ssm_fit
} # }
```
