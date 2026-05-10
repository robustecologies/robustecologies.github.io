# Multivariate normal distribution: Cholesky parameterization

Density and random generation for the multivariate normal distribution:
cholesky parameterization.

See Details.

See Details.

## Usage

``` r
dmvnc(x, mu, U, log = FALSE)

rmvnc(n = 1, mu = rep(0, k), U)
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

`dmvnc` gives the density, and `rmvnc` generates random deviates.

See Details.

See Details.

## Details

Implementation of `dmvnc`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `rmvnc`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dmvnc
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rmvnc
} # }
```
