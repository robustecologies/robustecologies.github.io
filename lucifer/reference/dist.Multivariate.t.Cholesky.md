# Multivariate t distribution: Cholesky parameterization

Density and random generation for the multivariate t distribution:
cholesky parameterization.

See Details.

See Details.

## Usage

``` r
dmvtc(x, mu, U, df = Inf, log = FALSE)

rmvtc(n = 1, mu = rep(0, k), U, df = Inf)
```

## Arguments

- x:

  data matrix or vector.

- mu:

  mean vector or matrix.

- U:

  upper-triangular Cholesky factor.

- df:

  positive degrees of freedom parameter.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

- n:

  number of observations.

## Value

`dmvtc` gives the density, and `rmvtc` generates random deviates.

See Details.

See Details.

## Details

Implementation of `dmvtc`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `rmvtc`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dmvtc
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rmvtc
} # }
```
