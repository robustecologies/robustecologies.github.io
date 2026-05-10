# Integrated autocorrelation time

The `IAT` function estimates integrated autocorrelation time, which is
the computational inefficiency of a continuous chain or MCMC sampler.
IAT is also called the IACT, ACT, autocorrelation time, autocovariance
time, correlation time, or inefficiency factor. A lower value of `IAT`
is better. `IAT` is an MCMC diagnostic that estimates the number of
iterations, on average, for an independent sample to be drawn, given a
continuous chain or Markov chain. Put another way, `IAT` is the number
of correlated samples with the same variance-reducing power as one
independent sample.

IAT is a univariate function. A multivariate form is not included.

## Usage

``` r
IAT(x)
```

## Arguments

- x:

  a vector of samples from a chain.

## Value

The integrated autocorrelation time of a chain.

## Details

`IAT` is an MCMC diagnostic that is often used to compare continuous
chains of MCMC samplers for computational inefficiency, where the
sampler with the lowest `IAT`s is the most efficient sampler. Otherwise,
chains may be compared within a model, such as with the output of
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
to learn about the inefficiency of the continuous chain. For more
information on comparing MCMC algorithmic inefficiency, see the
[`Arena`](https://robustecologies.github.io/lucifer/reference/Arena.md)
function.

`IAT` is also estimated in the
[`PosteriorChecks`](https://robustecologies.github.io/lucifer/reference/PosteriorChecks.md)
function. `IAT` is usually applied to a stationary, continuous chain
after discarding burn-in iterations (see
[`burnin`](https://robustecologies.github.io/lucifer/reference/burnin.md)
for more information). The `IAT` of a continuous chain correlates with
the variability of the mean of the chain, and relates to Effective
Sample Size
([`ESS`](https://robustecologies.github.io/lucifer/reference/ESS.md))
and Monte Carlo Standard Error
([`MCSE`](https://robustecologies.github.io/lucifer/reference/MCSE.md)).

`IAT` and
[`ESS`](https://robustecologies.github.io/lucifer/reference/ESS.md) are
inversely related, though not perfectly, because each is estimated a
little differently. Given \\N\\ samples and taking autocorrelation into
account,
[`ESS`](https://robustecologies.github.io/lucifer/reference/ESS.md)
estimates a reduced number of \\M\\ samples. Conversely, `IAT` estimates
the number of autocorrelated samples, on average, required to produce
one independently drawn sample.

The `IAT` function is similar to the `IAT` function in the `Rtwalk`
package of Christen and Fox (2010), which is currently unavailable on
CRAN.

## References

Christen, J.A. and Fox, C. (2010). "A General Purpose Sampling Algorithm
for Continuous Distributions (the t-walk)". *Bayesian Analysis*, 5(2),
p. 263–282.

## See also

[`burnin`](https://robustecologies.github.io/lucifer/reference/burnin.md),
[`Compare`](https://rdrr.io/r/methods/S4groupGeneric.html),
[`ESS`](https://robustecologies.github.io/lucifer/reference/ESS.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`MCSE`](https://robustecologies.github.io/lucifer/reference/MCSE.md),
and
[`PosteriorChecks`](https://robustecologies.github.io/lucifer/reference/PosteriorChecks.md).

## Examples

``` r
if (FALSE) { # \dontrun{
theta <- rnorm(100)
IAT(theta)
} # }
```
