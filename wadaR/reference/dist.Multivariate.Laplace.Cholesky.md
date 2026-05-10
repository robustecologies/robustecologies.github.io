# Multivariate Laplace distribution: Cholesky parameterization

Density and random generation for the multivariate laplace distribution:
cholesky parameterization.

See Details.

See Details.

## Usage

``` r
dmvlc(x, mu, U, log = FALSE)

rmvlc(n, mu, U)
```

## Arguments

- x:

  data matrix or vector.

- mu:

  mean vector or matrix.

- U:

  upper-triangular Cholesky factor.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

- n:

  number of observations.

## Value

`dmvlc` gives the density, and `rmvlc` generates random deviates.

See Details.

See Details.

## Details

Implementation of `dmvlc`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `rmvlc`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dmvlc
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rmvlc
} # }
```
