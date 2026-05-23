# Asymmetric log-Laplace distribution

Density, distribution function, quantile function, and random generation
for the univariate asymmetric log-Laplace distribution with location
parameter `location`, scale parameter `scale`, and asymmetry parameter
`kappa`.

See Details.

See Details.

See Details.

See Details.

## Usage

``` r
dallaplace(x, location = 0, scale = 1, kappa = 1, log = FALSE)

pallaplace(q, location = 0, scale = 1, kappa = 1)

qallaplace(p, location = 0, scale = 1, kappa = 1)

rallaplace(n, location = 0, scale = 1, kappa = 1)
```

## Arguments

- x, q:

  vector of quantiles.

- location:

  location parameter \\\mu\\.

- scale:

  scale parameter \\\lambda\\, which must be positive.

- kappa:

  asymmetry or skewness parameter \\\kappa\\, which must be positive.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

- p:

  vector of probabilities.

- n:

  number of observations, which must be a positive integer that has
  length 1.

## Value

`dallaplace` gives the density, `pallaplace` gives the distribution
function, `qallaplace` gives the quantile function, and `rallaplace`
generates random deviates.

See Details.

See Details.

See Details.

See Details.

## Details

The univariate asymmetric log-Laplace distribution is derived from the
Laplace distribution. Multivariate and symmetric versions also exist.

Implementation of `dallaplace`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `pallaplace`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `qallaplace`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `rallaplace`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## References

Kozubowski, T.J. and Podgorski, K. (2003). "Log-Laplace Distributions".
*International Mathematical Journal*, 3, p. 467–495.

## See also

[`dalaplace`](https://robustecologies.github.io/lucifer/reference/dist.Asymmetric.Laplace.md),
[`dlaplace`](https://robustecologies.github.io/lucifer/reference/dist.Laplace.md),
[`dllaplace`](https://robustecologies.github.io/lucifer/reference/dist.Log.Laplace.md)

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
x <- dallaplace(1, 0, 1, 1)
x <- pallaplace(1, 0, 1, 1)
x <- qallaplace(0.5, 0, 1, 1)
x <- rallaplace(100, 0, 1, 1)
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dallaplace
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving pallaplace
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving qallaplace
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rallaplace
} # }
```
