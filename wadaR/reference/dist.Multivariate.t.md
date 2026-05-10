# Multivariate t distribution

Density and random generation for the multivariate t distribution.

See Details.

See Details.

## Usage

``` r
dmvt(x, mu, S, df = Inf, log = FALSE)

rmvt(n = 1, mu = rep(0, k), S, df = Inf)
```

## Arguments

- x:

  data matrix or vector.

- mu:

  mean vector or matrix.

- S:

  \\K \times K\\ positive-definite covariance matrix.

- df:

  positive degrees of freedom parameter.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

- n:

  number of observations.

## Value

`dmvt` gives the density, and `rmvt` generates random deviates.

See Details.

See Details.

## Details

Implementation of `dmvt`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `rmvt`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dmvt
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rmvt
} # }
```
