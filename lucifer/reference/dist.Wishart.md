# Wishart distribution

Density and related functions for the wishart distribution.

See Details.

See Details.

## Usage

``` r
dwishart(Omega, nu, S, log = FALSE)

rwishart(nu, S)
```

## Arguments

- Omega:

  symmetric, positive-definite \\k \times k\\ matrix \\\Omega\\.

- nu:

  scalar degrees of freedom \\\nu\\.

- S:

  symmetric, positive-semidefinite \\k \times k\\ scale matrix.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

## Value

`dwishart` gives the density, and other functions provide related
computations.

See Details.

See Details.

## Details

Implementation of `dwishart`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `rwishart`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dwishart
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rwishart
} # }
```
