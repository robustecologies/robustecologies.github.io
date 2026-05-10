# Multivariate normal distribution

Density and random generation for the multivariate normal distribution.

See Details.

See Details.

## Usage

``` r
dmvn(x, mu, Sigma, log = FALSE)

rmvn(n = 1, mu = rep(0, k), Sigma)
```

## Arguments

- x:

  data matrix or vector.

- mu:

  mean vector or matrix.

- Sigma:

  \\K \times K\\ positive-definite covariance matrix.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

- n:

  number of observations.

## Value

`dmvn` gives the density, and `rmvn` generates random deviates.

See Details.

See Details.

## Details

Implementation of `dmvn`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `rmvn`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dmvn
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rmvn
} # }
```
