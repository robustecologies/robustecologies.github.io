# Inverse matrix gamma distribution

Density and random generation for the inverse matrix gamma distribution.

See Details.

See Details.

## Usage

``` r
dinvmatrixgamma(X, alpha, beta, Psi, log = FALSE)

rinvmatrixgamma(alpha, beta, Psi)
```

## Arguments

- X:

  \\k \times k\\ positive-definite covariance matrix.

- alpha:

  scalar shape parameter (degrees of freedom) \\\alpha\\, must be
  greater than 2.

- beta:

  scalar, positive-only scale parameter \\\beta\\.

- Psi:

  \\k \times k\\ positive-definite scale matrix.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

## Value

`dinvmatrixgamma` gives the density or other results.

See Details.

See Details.

## Details

Implementation of `dinvmatrixgamma`. Refer to the package vignettes and
the cited references for a complete algorithmic and mathematical
description.

Implementation of `rinvmatrixgamma`. Refer to the package vignettes and
the cited references for a complete algorithmic and mathematical
description.

## See also

[`dinvgamma`](https://robustecologies.github.io/lucifer/reference/dist.Inverse.Gamma.md),
[`dinvwishart`](https://robustecologies.github.io/lucifer/reference/dist.Inverse.Wishart.md)

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dinvmatrixgamma
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rinvmatrixgamma
} # }
```
