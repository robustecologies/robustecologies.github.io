# LASSO distribution

Density and random generation for the Bayesian LASSO prior distribution.

See Details.

See Details.

## Usage

``` r
dlasso(x, sigma, tau, lambda, a = 1, b = 1, log = FALSE)

rlasso(n, sigma, tau, lambda, a = 1, b = 1)
```

## Arguments

- x:

  location vector of length \\J\\ at which to evaluate density.

- sigma:

  positive-only scalar hyperparameter \\\sigma\\.

- tau:

  positive-only vector of hyperparameters \\\tau\\ of length \\J\\.

- lambda:

  positive-only scalar hyperhyperparameter \\\lambda\\.

- a, b:

  positive-only scalar hyperhyperhyperparameters for gamma-distributed
  \\\lambda\\.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

- n:

  number of observations.

## Value

`dlasso` gives the density or other results.

See Details.

See Details.

## Details

Implementation of `dlasso`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `rlasso`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

## See also

[`dhs`](https://robustecologies.github.io/lucifer/reference/dist.Horseshoe.md)

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dlasso
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rlasso
} # }
```
