# Normal-Wishart distribution

Density and related functions for the normal-wishart distribution.

See Details.

See Details.

## Usage

``` r
dnormwishart(mu, mu0, lambda, Omega, S, nu, log = FALSE)

rnormwishart(n = 1, mu0, lambda, S, nu)
```

## Arguments

- mu:

  data or parameters vector of length \\k\\ or matrix with \\k\\
  columns.

- mu0:

  mean vector \\\mu_0\\.

- lambda:

  positive-only scalar.

- Omega:

  \\k \times k\\ precision matrix.

- S:

  symmetric, positive-semidefinite \\k \times k\\ scale matrix.

- nu:

  scalar degrees of freedom \\\nu\\.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

- n:

  number of random draws.

## Value

`dnormwishart` gives the density, and other functions provide related
computations.

See Details.

See Details.

## Details

Implementation of `dnormwishart`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `rnormwishart`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dnormwishart
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rnormwishart
} # }
```
