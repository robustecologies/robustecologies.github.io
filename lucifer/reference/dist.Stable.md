# Stable distribution

Density and random generation for the stable distribution using Nolan's
0-parameterization.

See Details.

See Details.

## Usage

``` r
dstable(x, alpha = 2, beta = 0, gamma = 1, delta = 0, log = FALSE)

rstable(n, alpha = 2, beta = 0, gamma = 1, delta = 0)
```

## Arguments

- x:

  vector of quantiles.

- alpha:

  vector of stability parameters in \\(0, 2\]\\.

- beta:

  vector of skewness parameters in \\\[-1, 1\]\\.

- gamma:

  vector of scale parameters, which must be positive.

- delta:

  vector of location parameters.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

- n:

  number of observations, which must be a positive integer of length 1.

## Value

`dstable` gives the density and `rstable` generates random deviates.

See Details.

See Details.

## Details

The stable distribution is a four-parameter family that generalizes the
Gaussian distribution. The stability parameter \\\alpha\\ controls tail
heaviness (smaller values produce heavier tails), \\\beta\\ controls
skewness, \\\gamma\\ is the scale, and \\\delta\\ is the location.
Special cases include the Gaussian (\\\alpha = 2\\), Cauchy (\\\alpha =
1, \beta = 0\\), and Levy (\\\alpha = 0.5, \beta = 1\\) distributions.

The density has no closed form for general \\\alpha\\ and is computed
via numerical integration of the characteristic function. Random
generation uses the Chambers-Mallows-Stuck algorithm.

Implementation of `dstable`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `rstable`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## References

Nolan, J.P. (2020). *Stable Distributions: Models for Heavy-Tailed
Data*. Birkhauser. ISBN 978-3-030-52914-0.

Chambers, J.M., Mallows, C.L., and Stuck, B.W. (1976). "A Method for
Simulating Stable Random Variables". *Journal of the American
Statistical Association*, 71(354), p. 340–344.

## See also

[`dnorm`](https://rdrr.io/r/stats/Normal.html),
[`dcauchy`](https://rdrr.io/r/stats/Cauchy.html)

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Cauchy distribution (alpha=1, beta=0)
x <- seq(-10, 10, length.out = 200)
d <- dstable(x, alpha = 1, beta = 0, gamma = 1, delta = 0)
plot(x, d, type = "l")
# Heavy-tailed stable
r <- rstable(1000, alpha = 1.5, beta = 0.5, gamma = 1, delta = 0)
hist(r, breaks = 100, freq = FALSE)
} # }

if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dstable
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rstable
} # }
```
