# Inverse Wishart distribution

Density and random number generation for the inverse Wishart
distribution.

See Details.

See Details.

## Usage

``` r
dinvwishart(Sigma, nu, S, log = FALSE)

rinvwishart(nu, S)
```

## Arguments

- Sigma:

  symmetric, positive-definite \\k \times k\\ matrix \\\Sigma\\.

- nu:

  scalar degrees of freedom \\\nu\\.

- S:

  symmetric, positive-semidefinite \\k \times k\\ scale matrix.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

## Value

`dinvwishart` gives the density or other results.

See Details.

See Details.

## Details

Implementation of `dinvwishart`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `rinvwishart`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## See also

[`dhuangwand`](https://robustecologies.github.io/lucifer/reference/dist.HuangWand.md),
[`dinvwishartc`](https://robustecologies.github.io/lucifer/reference/dist.Inverse.Wishart.Cholesky.md),
[`dwishart`](https://robustecologies.github.io/lucifer/reference/dist.Wishart.md)

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dinvwishart
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rinvwishart
} # }
```
