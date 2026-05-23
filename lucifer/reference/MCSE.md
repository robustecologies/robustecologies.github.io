# Monte Carlo standard error

Monte Carlo Standard Error (MCSE) is an estimate of the inaccuracy of
Monte Carlo samples, usually regarding the expectation of posterior
samples, \\\mathrm{E}(\theta)\\, from Monte Carlo or Markov chain Monte
Carlo (MCMC) algorithms, such as with the
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
function. MCSE approaches zero as the number of independent posterior
samples approaches infinity. MCSE is essentially a standard deviation
around the posterior mean of the samples, \\\mathrm{E}(\theta)\\, due to
uncertainty associated with using an MCMC algorithm, or Monte Carlo
methods in general.

The acceptable size of the MCSE depends on the acceptable uncertainty
associated around the marginal posterior mean, \\\mathrm{E}(\theta)\\,
and the goal of inference. It has been argued that MCSE is generally
unimportant when the goal of inference is \\\theta\\ rather than
\\\mathrm{E}(\theta)\\ (Gelman et al., 2004, p. 277), and that a
sufficient
[`ESS`](https://robustecologies.github.io/lucifer/reference/ESS.md) is
more important. Others perceive MCSE to be a vital part of reporting any
Bayesian model, and as a stopping rule (Flegal et al., 2008).

In
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
MCSE is part of the posterior summaries because it is easy to estimate,
and Lucifer prefers to continue updating until each MCSE is less than
6.27% of its associated marginal posterior standard deviation (for more
information on this stopping rule, see the
[`Consort`](https://robustecologies.github.io/lucifer/reference/Consort.md)
function), since MCSE has been demonstrated to be an excellent stopping
rule.

Acceptable error may be specified, if known, in the `MCSS` (Monte Carlo
Sample Size) function to estimate the required number of posterior
samples.

`MCSE` is a univariate function that is often applied to each marginal
posterior distribution. A multivariate form is not included. By chance
alone due to multiple independent tests, 5% of the parameters should
indicate unacceptable MCSEs, even when acceptable. Assessing convergence
is difficult.

## Usage

``` r
MCSE(x, method = "IMPS", batch.size = "sqrt", warn = FALSE)
```

## Arguments

- x:

  a vector of posterior samples for which MCSE will be estimated.

- method:

  the method of MCSE estimation. Defaults to Geyer's `"IMPS"` method.
  Optional methods include `"sample.variance"` and `"batch.mean"`. Note
  that `"batch.mean"` is recommended only when the number of posterior
  samples is at least 1,000.

- batch.size:

  either the size of the batches (accepting a numerical argument) or the
  method of creating the size of batches, which is either `"sqrt"` or
  `"cuberoot"`, and refers to the length of `x`. The default argument is
  `"sqrt"`. This argument corresponds only with `method="batch.means"`.

- warn:

  logical. If `warn=TRUE`, then a warning is provided with
  `method="batch.means"` whenever posterior sample size is less than
  1,000, or a warning is produced when more autocovariance is
  recommended with `method="IMPS"`.

## Value

The estimated Monte Carlo standard error.

## Details

The default method for estimating MCSE is Geyer's Initial Monotone
Positive Sequence (IMPS) estimator (Geyer, 1992), which takes the
asymptotic variance into account and is time-series based. This method
goes by other names, such as Initial Positive Sequence (IPS).

The simplest method for estimating MCSE is to modify the formula for
standard error, \\\sigma(\textbf{x}) / \sqrt{N}\\, to account for
non-independence in the sequence \\\textbf{x}\\ of posterior samples.
Non-independence is estimated with the `ESS` function for Effective
Sample Size (see the
[`ESS`](https://robustecologies.github.io/lucifer/reference/ESS.md)
function for more details), where \\M = ESS(\textbf{x})\\, and MCSE is
\\\sigma(\textbf{x}) / \sqrt{M}\\. Although this is the fastest and
easiest method of estimation, it does not incorporate an estimate of the
asymptotic variance of \\\textbf{x}\\.

The batch means method (Jones et al., 2006; Flegal et al., 2008)
separates elements of \\\textbf{x}\\ into batches and estimates MCSE as
a function of multiple batches. This method is excellent, but is not
recommended when the number of posterior samples is less than 1,000.
These journal articles also assert that MCSE is a better stopping rule
than MCMC convergence diagnostics.

The `MCSS` function estimates the required number of posterior samples,
given the user-specified acceptable error, posterior samples `x`, and
the observed variance (rather than asymptotic variance). Due to the
observed variance, this is a rough estimate.

## References

Flegal, J.M., Haran, M., and Jones, G.L. (2008). "Markov chain Monte
Carlo: Can We Trust the Third Significant Figure?". *Statistical
Science*, 23, p. 250–260.

Gelman, A., Carlin, J., Stern, H., and Rubin, D. (2004). "Bayesian Data
Analysis, Texts in Statistical Science, 2nd ed.". Chapman and Hall,
London.

Geyer, C.J. (1992). "Practical Markov Chain Monte Carlo". *Statistical
Science*, 7, 4, p. 473–483.

Jones, G.L., Haran, M., Caffo, B.S., and Neath, R. (2006). "Fixed-Width
Output Analysis for Markov chain Monte Carlo". *Journal of the American
Statistical Association*, 101(1), p. 1537–1547.

## See also

[`Consort`](https://robustecologies.github.io/lucifer/reference/Consort.md),
[`ESS`](https://robustecologies.github.io/lucifer/reference/ESS.md), and
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

## Examples

``` r
if (FALSE) { # \dontrun{
x <- rnorm(1000)
MCSE(x)
MCSE(x, method="batch.means")
MCSE(x, method="sample.variance")
} # }
```
