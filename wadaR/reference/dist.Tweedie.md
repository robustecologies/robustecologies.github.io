# Tweedie distribution

Density and random generation for the Tweedie compound Poisson-gamma
distribution.

See Details.

See Details.

## Usage

``` r
dtweedie(x, mu = 1, phi = 1, power = 1.5, log = FALSE)

rtweedie(n, mu = 1, phi = 1, power = 1.5)
```

## Arguments

- x:

  vector of quantiles (non-negative).

- mu:

  vector of mean parameters, which must be positive.

- phi:

  vector of dispersion parameters, which must be positive.

- power:

  vector of power parameters, which must be in \\(1, 2)\\.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

- n:

  number of observations, which must be a positive integer of length 1.

## Value

`dtweedie` gives the density and `rtweedie` generates random deviates.

See Details.

See Details.

## Details

The Tweedie distribution with power parameter \\p \in (1, 2)\\ is a
compound Poisson-gamma distribution. A Tweedie random variable \\Y\\ can
be represented as \\Y = \sum\_{j=1}^{N} G_j\\ where \\N \sim
\mathrm{Poisson}(\lambda)\\ and the \\G_j\\ are iid gamma random
variables. The distribution has a point mass at zero (when \\N = 0\\)
and a continuous density on \\(0, \infty)\\.

The density for \\y \> 0\\ is computed via the series expansion of Dunn
and Smyth (2005), using log-sum-exp for numerical stability. Random
generation uses the compound Poisson-gamma representation directly.

The Tweedie distribution is central to generalized linear models for
insurance claims data and ecological count data with exact zeros.

Implementation of `dtweedie`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `rtweedie`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## References

Dunn, P.K. and Smyth, G.K. (2005). "Series Evaluation of Tweedie
Exponential Dispersion Model Densities". *Statistics and Computing*,
15(4), p. 267–280.

Jorgensen, B. (1997). *The Theory of Dispersion Models*. Chapman and
Hall. ISBN 978-0-412-99711-4.

## See also

[`dgamma`](https://rdrr.io/r/stats/GammaDist.html),
[`dpois`](https://rdrr.io/r/stats/Poisson.html)

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
x <- seq(0, 10, length.out = 200)
d <- dtweedie(x, mu = 2, phi = 1, power = 1.5)
plot(x, d, type = "l")
r <- rtweedie(1000, mu = 2, phi = 1, power = 1.5)
hist(r, breaks = 50, freq = FALSE)
} # }

if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dtweedie
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rtweedie
} # }
```
