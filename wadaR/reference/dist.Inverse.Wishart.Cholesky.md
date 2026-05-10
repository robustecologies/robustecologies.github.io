# Inverse Wishart distribution: Cholesky parameterization

Density and random number generation for the inverse Wishart
distribution with the Cholesky parameterization.

See Details.

See Details.

## Usage

``` r
dinvwishartc(U, nu, S, log = FALSE)

rinvwishartc(nu, S)
```

## Arguments

- U:

  upper-triangular \\k \times k\\ Cholesky factor of covariance matrix
  \\\Sigma\\.

- nu:

  scalar degrees of freedom \\\nu\\.

- S:

  symmetric, positive-semidefinite \\k \times k\\ scale matrix.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

## Value

`dinvwishartc` gives the density or other results.

See Details.

See Details.

## Details

Implementation of `dinvwishartc`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `rinvwishartc`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## See also

[`chol`](https://rdrr.io/r/base/chol.html),
[`dinvwishart`](https://robustecologies.github.io/lucifer/reference/dist.Inverse.Wishart.md),
[`dwishartc`](https://robustecologies.github.io/lucifer/reference/dist.Wishart.Cholesky.md)

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dinvwishartc
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rinvwishartc
} # }
```
