# Asymmetric Laplace distribution: univariate

Density, distribution function, quantile function, and random generation
for the univariate asymmetric Laplace distribution with location
parameter `location`, scale parameter `scale`, and asymmetry parameter
`kappa`.

See Details.

See Details.

See Details.

See Details.

## Usage

``` r
dalaplace(x, location = 0, scale = 1, kappa = 1, log = FALSE)

palaplace(q, location = 0, scale = 1, kappa = 1)

qalaplace(p, location = 0, scale = 1, kappa = 1)

ralaplace(n, location = 0, scale = 1, kappa = 1)
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

`dalaplace` gives the density, `palaplace` gives the distribution
function, `qalaplace` gives the quantile function, and `ralaplace`
generates random deviates.

See Details.

See Details.

See Details.

See Details.

## Details

The asymmetric Laplace of Kotz, Kozubowski, and Podgorski (2001), also
referred to as AL, is an extension of the univariate symmetric Laplace
distribution to allow for skewness. It is parameterized according to
three parameters: location parameter \\\mu\\, scale parameter
\\\lambda\\, and asymmetry or skewness parameter \\\kappa\\. The special
case \\\kappa=1\\ is the symmetric Laplace distribution. Values of
\\\kappa\\ in the intervals \\(0,1)\\ and \\(1,\infty)\\ correspond to
positive (right) and negative (left) skewness, respectively.

Implementation of `dalaplace`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `palaplace`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `qalaplace`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `ralaplace`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## References

Kotz, S., Kozubowski, T.J., and Podgorski, K. (2001). "The Laplace
Distribution and Generalizations: a Revisit with Applications to
Communications, Economics, Engineering, and Finance". Boston:
Birkhauser.

## See also

[`dlaplace`](https://robustecologies.github.io/lucifer/reference/dist.Laplace.md),
[`dallaplace`](https://robustecologies.github.io/lucifer/reference/dist.Asymmetric.Log.Laplace.md)

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
x <- dalaplace(1, 0, 1, 1)
x <- palaplace(1, 0, 1, 1)
x <- qalaplace(0.5, 0, 1, 1)
x <- ralaplace(100, 0, 1, 1)
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dalaplace
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving palaplace
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving qalaplace
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving ralaplace
} # }
```
