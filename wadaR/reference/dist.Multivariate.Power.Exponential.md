# Multivariate power exponential distribution

Density and random generation for the multivariate power exponential
distribution.

See Details.

See Details.

## Usage

``` r
dmvpe(x = c(0, 0), mu = c(0, 0), Sigma = diag(2), kappa = 1, log = FALSE)

rmvpe(n, mu = c(0, 0), Sigma = diag(2), kappa = 1)
```

## Arguments

- x:

  data matrix or vector.

- mu:

  mean vector or matrix.

- Sigma:

  \\K \times K\\ positive-definite covariance matrix.

- kappa:

  positive kurtosis parameter.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

- n:

  number of observations.

## Value

`dmvpe` gives the density, and `rmvpe` generates random deviates.

See Details.

See Details.

## Details

Implementation of `dmvpe`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `rmvpe`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dmvpe
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rmvpe
} # }
```
