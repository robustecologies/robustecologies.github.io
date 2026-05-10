# Matern correlation and covariance functions

Compute the Matern correlation for a vector of distances, or the full
Matern covariance matrix from a distance matrix.

S3 method: apply `matern()` to objects of class `corr`.

S3 method: apply `matern()` to objects of class `cov`.

## Usage

``` r
matern.corr(h, nu, rho)

matern.cov(dist.mat, sigma2, nu, rho)
```

## Arguments

- h:

  vector of non-negative distances.

- nu:

  smoothness parameter, which must be positive. Common values: \\\nu =
  0.5\\ gives exponential correlation, \\\nu = 1.5\\ gives the Matern
  3/2 model, and \\\nu \to \infty\\ gives squared exponential
  (Gaussian).

- rho:

  range parameter, which must be positive. Roughly, the distance at
  which the correlation drops to about 0.13.

- dist.mat:

  an \\n \times n\\ distance matrix.

- sigma2:

  marginal variance (sill), which must be positive.

## Value

`matern.corr` returns a vector of correlations. `matern.cov` returns an
\\n \times n\\ covariance matrix.

See Details.

See Details.

## Details

The Matern correlation function is \$\$C(h; \nu, \rho) =
\frac{2^{1-\nu}}{\Gamma(\nu)} \left(\frac{\sqrt{2\nu}\\
h}{\rho}\right)^\nu K\_\nu\\\left(\frac{\sqrt{2\nu}\\
h}{\rho}\right)\$\$ where \\K\_\nu\\ is the modified Bessel function of
the second kind. The function `matern.corr` computes this for a vector
of distances; `matern.cov` constructs the full covariance matrix
\\\sigma^2 C(h\_{ij}; \nu, \rho)\\ from a distance matrix, with OpenMP
parallelization for large matrices.

Implementation of `matern.corr`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `matern.cov`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## References

Stein, M.L. (1999). *Interpolation of Spatial Data: Some Theory for
Kriging*. Springer. ISBN 978-0-387-98629-6.

Rasmussen, C.E. and Williams, C.K.I. (2006). *Gaussian Processes for
Machine Learning*. MIT Press. ISBN 978-0-262-18253-9.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
h <- seq(0, 5, length.out = 100)
c1 <- matern.corr(h, nu = 0.5, rho = 1)
c2 <- matern.corr(h, nu = 1.5, rho = 1)
c3 <- matern.corr(h, nu = 2.5, rho = 1)
plot(h, c1, type = "l", ylab = "Correlation")
lines(h, c2, col = 2)
lines(h, c3, col = 4)

# Covariance matrix for 20 spatial locations
locs <- cbind(runif(20), runif(20))
D <- as.matrix(dist(locs))
Sigma <- matern.cov(D, sigma2 = 1, nu = 1.5, rho = 0.3)
} # }

if (FALSE) { # \dontrun{
## see package vignettes for full examples involving matern.corr
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving matern.cov
} # }
```
