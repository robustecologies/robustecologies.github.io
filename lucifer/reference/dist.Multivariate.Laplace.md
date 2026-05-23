# Multivariate Laplace distribution

Density and random generation for the multivariate laplace distribution.

See Details.

See Details.

## Usage

``` r
dmvl(x, mu, Sigma, log = FALSE)

rmvl(n, mu, Sigma)
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

`dmvl` gives the density, and `rmvl` generates random deviates.

See Details.

See Details.

## Details

Implementation of `dmvl`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `rmvl`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dmvl
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rmvl
} # }
```
