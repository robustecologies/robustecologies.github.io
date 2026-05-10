# Generalized Poisson distribution

Density function for the univariate, discrete, generalized Poisson
distribution with location parameter \\\lambda\\ and scale parameter
\\\omega\\.

See Details.

## Usage

``` r
dgpois(x, lambda = 0, omega = 0, log = FALSE)
```

## Arguments

- x:

  vector of quantiles.

- lambda:

  parameter \\\lambda\\.

- omega:

  parameter \\\omega\\, which should be in the interval \[0,1) for
  positive counts.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

## Value

`dgpois` gives the density.

See Details.

## Details

The generalized Poisson distribution (Consul, 1989) is also called the
Lagrangian Poisson distribution. The simple Poisson distribution is a
special case when \\\omega = 0\\.

Implementation of `dgpois`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

## References

Consul, P. (1989). "Generalized Poisson Distribution: Properties and
Applications". Marcel Decker: New York, NY.

Ntzoufras, I., Katsis, A., and Karlis, D. (2005). "Bayesian Assessment
of the Distribution of Insurance Claim Counts using Reversible Jump
MCMC". *North American Actuarial Journal*, 9, p. 90–108.

## See also

[`dnbinom`](https://rdrr.io/r/stats/NegBinomial.html),
[`dpois`](https://rdrr.io/r/stats/Poisson.html)

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
y <- rpois(100, 5)
lambda <- rpois(100, 5)
x <- dgpois(y, lambda, 0.5)
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dgpois
} # }
```
