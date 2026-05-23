# Multivariate t distribution: precision-Cholesky parameterization

Density and random generation for the multivariate t distribution:
precision-cholesky parameterization.

See Details.

See Details.

## Usage

``` r
dmvtpc(x, mu, U, nu = Inf, log = FALSE)

rmvtpc(n = 1, mu, U, nu = Inf)
```

## Arguments

- x:

  data matrix or vector.

- mu:

  mean vector or matrix.

- U:

  upper-triangular Cholesky factor of precision matrix.

- nu:

  positive degrees of freedom parameter.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

- n:

  number of observations.

## Value

`dmvtpc` gives the density, and `rmvtpc` generates random deviates.

See Details.

See Details.

## Details

Implementation of `dmvtpc`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `rmvtpc`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dmvtpc
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rmvtpc
} # }
```
