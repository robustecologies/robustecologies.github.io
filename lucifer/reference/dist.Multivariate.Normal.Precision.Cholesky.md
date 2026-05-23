# Multivariate normal distribution: precision-Cholesky parameterization

Density and random generation for the multivariate normal distribution:
precision-cholesky parameterization.

See Details.

See Details.

## Usage

``` r
dmvnpc(x, mu, U, log = FALSE)

rmvnpc(n = 1, mu = rep(0, k), U)
```

## Arguments

- x:

  data matrix or vector.

- mu:

  mean vector or matrix.

- U:

  upper-triangular Cholesky factor of precision matrix.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

- n:

  number of observations.

## Value

`dmvnpc` gives the density, and `rmvnpc` generates random deviates.

See Details.

See Details.

## Details

Implementation of `dmvnpc`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `rmvnpc`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dmvnpc
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rmvnpc
} # }
```
