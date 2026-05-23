# Exponentially-modified Gaussian distribution

Density, distribution function, and random generation for the
exponentially-modified Gaussian (ExGaussian) distribution, formed as the
convolution of a normal and an exponential random variable.

See Details.

See Details.

See Details.

## Usage

``` r
dexgaussian(x, mu = 0, sigma = 1, lambda = 1, log = FALSE)

pexgaussian(q, mu = 0, sigma = 1, lambda = 1)

rexgaussian(n, mu = 0, sigma = 1, lambda = 1)
```

## Arguments

- x, q:

  vector of quantiles.

- mu:

  vector of means for the Gaussian component.

- sigma:

  vector of standard deviations for the Gaussian component, which must
  be positive.

- lambda:

  vector of rate parameters for the exponential component, which must be
  positive.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

- n:

  number of observations, which must be a positive integer of length 1.

## Value

`dexgaussian` gives the density, `pexgaussian` gives the distribution
function, and `rexgaussian` generates random deviates.

See Details.

See Details.

See Details.

## Details

The ExGaussian distribution arises as the sum of independent normal and
exponential random variables: if \\X \sim N(\mu, \sigma^2)\\ and \\Y
\sim \mathrm{Exp}(\lambda)\\, then \\Z = X + Y\\ follows an ExGaussian
distribution with density \$\$f(z) = \frac{\lambda}{2}
\exp\left(\frac{\lambda}{2}(2\mu + \lambda\sigma^2 - 2z)\right)
\mathrm{erfc}\left(\frac{\mu + \lambda\sigma^2 - z}{\sigma\sqrt{2}}
\right)\$\$ where \\\mathrm{erfc}\\ is the complementary error function.
The mean is \\\mu + 1/\lambda\\ and variance is \\\sigma^2 +
1/\lambda^2\\. The distribution is widely used in reaction time modeling
and chromatographic peak fitting.

Implementation of `dexgaussian`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `pexgaussian`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `rexgaussian`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## References

Luce, R.D. (1986). *Response Times: Their Role in Inferring Elementary
Mental Organization*. Oxford University Press.

Grushka, E. (1972). "Characterization of Exponentially Modified Gaussian
Peaks in Chromatography". *Analytical Chemistry*, 44(11), p. 1733–1738.

## See also

[`dnorm`](https://rdrr.io/r/stats/Normal.html)

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
x <- seq(-2, 10, length.out = 200)
d <- dexgaussian(x, mu = 0, sigma = 1, lambda = 1)
plot(x, d, type = "l")
r <- rexgaussian(1000, mu = 0, sigma = 1, lambda = 1)
hist(r, breaks = 50, freq = FALSE)
} # }

if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dexgaussian
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving pexgaussian
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rexgaussian
} # }
```
