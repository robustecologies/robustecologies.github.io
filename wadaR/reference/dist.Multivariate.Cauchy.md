# Multivariate Cauchy distribution

Density and random generation for the multivariate cauchy distribution.

See Details.

See Details.

## Usage

``` r
dmvc(x, mu, S, log = FALSE)

rmvc(n = 1, mu = rep(0, k), S)
```

## Arguments

- x:

  \\N \times K\\ data matrix or vector of length \\K\\.

- mu:

  mean vector or matrix.

- S:

  \\K \times K\\ positive-definite covariance matrix.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

- n:

  number of observations.

## Value

`dmvc` gives the density, and `rmvc` generates random deviates.

See Details.

See Details.

## Details

Implementation of `dmvc`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `rmvc`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dmvc
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rmvc
} # }
```
