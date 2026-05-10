# Matrix-t distribution

Density for the matrix-variate t distribution.

See Details.

## Usage

``` r
dmatrixt(X, M, Omega, Sigma, nu, log = FALSE)
```

## Arguments

- X:

  an \\n \times p\\ data matrix.

- M:

  an \\n \times p\\ location matrix.

- Omega:

  an \\n \times n\\ positive-definite among-row scale matrix.

- Sigma:

  a \\p \times p\\ positive-definite among-column scale matrix.

- nu:

  degrees of freedom, which must be positive.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

## Value

`dmatrixt` gives the (log) density as a scalar.

See Details.

## Details

The matrix-variate t distribution generalizes the multivariate t
distribution to matrix-valued random variables. An \\n \times p\\ random
matrix \\X\\ follows a matrix-t distribution with location \\M\\,
among-row scale \\\Omega\\, among-column scale \\\Sigma\\, and degrees
of freedom \\\nu\\ if its density is \$\$f(X) \propto \|\Omega\|^{-n/2}
\|\Sigma\|^{-p/2} \|I_n + \frac{1}{\nu} \Omega^{-1}(X -
M)\Sigma^{-1}(X - M)^T\| ^{-(\nu + n + p - 1)/2}\$\$ The matrix-t arises
as a scale mixture of matrix normals with an inverse Wishart mixing
distribution and is used for robust covariance estimation.

Implementation of `dmatrixt`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## References

Gupta, A.K. and Nagar, D.K. (1999). *Matrix Variate Distributions*.
Chapman and Hall/CRC. ISBN 978-1-58488-046-2.

## See also

[`dmatrixnorm`](https://robustecologies.github.io/lucifer/reference/dist.Matrix.Normal.md),
[`dwishart`](https://robustecologies.github.io/lucifer/reference/dist.Wishart.md)

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
X <- matrix(rnorm(6), 2, 3)
M <- matrix(0, 2, 3)
Omega <- diag(2)
Sigma <- diag(3)
d <- dmatrixt(X, M, Omega, Sigma, nu = 5)
d_log <- dmatrixt(X, M, Omega, Sigma, nu = 5, log = TRUE)
} # }

if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dmatrixt
} # }
```
