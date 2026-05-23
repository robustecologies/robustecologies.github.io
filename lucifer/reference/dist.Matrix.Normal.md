# Matrix normal distribution

Density and random number generation for the matrix normal distribution.

See Details.

See Details.

## Usage

``` r
dmatrixnorm(X, M, U, V, log = FALSE)

rmatrixnorm(M, U, V)
```

## Arguments

- X:

  data or parameters matrix with \\n\\ rows and \\k\\ columns.

- M:

  mean matrix with \\n\\ rows and \\k\\ columns.

- U:

  \\n \times n\\ positive-definite scale matrix.

- V:

  \\k \times k\\ positive-definite scale matrix.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

## Value

`dmatrixnorm` gives the density or other results.

See Details.

See Details.

## Details

Implementation of `dmatrixnorm`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `rmatrixnorm`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## See also

[`dinvmatrixgamma`](https://robustecologies.github.io/lucifer/reference/dist.Inverse.Matrix.Gamma.md),
[`dmatrixgamma`](https://robustecologies.github.io/lucifer/reference/dist.Matrix.Gamma.md),
[`dmvn`](https://robustecologies.github.io/lucifer/reference/dist.Multivariate.Normal.md)

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dmatrixnorm
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rmatrixnorm
} # }
```
