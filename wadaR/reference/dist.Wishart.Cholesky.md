# Wishart distribution: Cholesky parameterization

Density and related functions for the wishart distribution: cholesky
parameterization.

See Details.

See Details.

## Usage

``` r
dwishartc(U, nu, S, log = FALSE)

rwishartc(nu, S)
```

## Arguments

- U:

  upper-triangular \\k \times k\\ Cholesky factor of precision matrix.

- nu:

  scalar degrees of freedom \\\nu\\.

- S:

  symmetric, positive-semidefinite \\k \times k\\ scale matrix.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

## Value

`dwishartc` gives the density, and other functions provide related
computations.

See Details.

See Details.

## Details

Implementation of `dwishartc`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `rwishartc`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dwishartc
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rwishartc
} # }
```
