# Multivariate Cauchy distribution: precision-Cholesky parameterization

Density and random generation for the multivariate cauchy distribution:
precision-cholesky parameterization.

See Details.

See Details.

## Usage

``` r
dmvcpc(x, mu, U, log = FALSE)

rmvcpc(n = 1, mu, U)
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

`dmvcpc` gives the density, and `rmvcpc` generates random deviates.

See Details.

See Details.

## Details

Implementation of `dmvcpc`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `rmvcpc`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dmvcpc
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rmvcpc
} # }
```
