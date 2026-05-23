# Asymmetric multivariate Laplace distribution

Density and random generation for the asymmetric multivariate Laplace
distribution with location and skew parameter \\\mu\\ and covariance
\\\Sigma\\.

See Details.

See Details.

## Usage

``` r
daml(x, mu, Sigma, log = FALSE)

raml(n, mu, Sigma)
```

## Arguments

- x:

  \\N \times K\\ matrix of data, or a vector of length \\K\\.

- mu:

  location and skew parameter \\\mu\\, a \\N \times K\\ matrix or vector
  of length \\K\\.

- Sigma:

  \\K \times K\\ positive-definite covariance matrix \\\Sigma\\.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

- n:

  number of observations, which must be a positive integer that has
  length 1.

## Value

`daml` gives the density, and `raml` generates random deviates.

See Details.

See Details.

## Details

The asymmetric multivariate Laplace distribution of Kotz, Kozubowski,
and Podgorski (2003) is a multivariate extension of the univariate
asymmetric Laplace distribution. When \\\mu=0\\, the density is the
symmetric multivariate Laplace of Anderson (1992).

Implementation of `daml`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `raml`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

## References

Kotz, S., Kozubowski, T.J., and Podgorski, K. (2003). "An Asymmetric
Multivariate Laplace Distribution". Working Paper.

## See also

[`dalaplace`](https://robustecologies.github.io/lucifer/reference/dist.Asymmetric.Laplace.md),
[`dmvl`](https://robustecologies.github.io/lucifer/reference/dist.Multivariate.Laplace.md)

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
x <- daml(c(1,2,3), c(0,1,2), diag(3))
X <- raml(1000, c(0,1,2), diag(3))
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving daml
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving raml
} # }
```
