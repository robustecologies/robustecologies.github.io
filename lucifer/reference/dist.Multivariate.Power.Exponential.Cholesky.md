# Multivariate power exponential distribution: Cholesky parameterization

Density and random generation for the multivariate power exponential
distribution: cholesky parameterization.

See Details.

See Details.

## Usage

``` r
dmvpec(x = c(0, 0), mu = c(0, 0), U, kappa = 1, log = FALSE)

rmvpec(n, mu = c(0, 0), U, kappa = 1)
```

## Arguments

- x:

  data matrix or vector.

- mu:

  mean vector or matrix.

- U:

  upper-triangular Cholesky factor.

- kappa:

  positive kurtosis parameter.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

- n:

  number of observations.

## Value

`dmvpec` gives the density, and `rmvpec` generates random deviates.

See Details.

See Details.

## Details

Implementation of `dmvpec`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `rmvpec`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dmvpec
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rmvpec
} # }
```
