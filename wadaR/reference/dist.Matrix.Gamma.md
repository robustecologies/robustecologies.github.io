# Matrix gamma distribution

Density and random generation for the matrix gamma distribution.

See Details.

See Details.

## Usage

``` r
dmatrixgamma(X, alpha, beta, Sigma, log = FALSE)

rmatrixgamma(alpha, beta, Sigma)
```

## Arguments

- X:

  \\k \times k\\ positive-definite precision matrix.

- alpha:

  scalar shape parameter (degrees of freedom), must be greater than 2.

- beta:

  scalar, positive-only scale parameter.

- Sigma:

  \\k \times k\\ positive-definite scale matrix.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

## Value

`dmatrixgamma` gives the density or other results.

See Details.

See Details.

## Details

Implementation of `dmatrixgamma`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `rmatrixgamma`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## See also

[`dgamma`](https://rdrr.io/r/stats/GammaDist.html),
[`dwishart`](https://robustecologies.github.io/lucifer/reference/dist.Wishart.md)

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dmatrixgamma
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rmatrixgamma
} # }
```
