# Gelman and Rubin's MCMC convergence diagnostic

Gelman and Rubin (1992) proposed a general approach to monitoring
convergence of MCMC output in which \\m \> 1\\ parallel chains are
updated with initial values that are overdispersed relative to each
target distribution, which must be normally distributed. Convergence is
diagnosed when the chains have 'forgotten' their initial values and the
output from all chains is indistinguishable. The `Gelman.Diagnostic`
function makes a comparison of within-chain and between-chain variances,
and is similar to a classical analysis of variance. A large deviation
between these two variances indicates non-convergence.

This diagnostic is popular as a stopping rule, though it requires
parallel chains. The
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
function supports parallel chains directly via the `Chains` argument. As
an alternative, the popular single-chain stopping rule is based on
[`MCSE`](https://robustecologies.github.io/lucifer/reference/MCSE.md).

## Usage

``` r
Gelman.Diagnostic(x, confidence = 0.95, transform = FALSE)
```

## Arguments

- x:

  an object of class `demonoid` with multiple chains (i.e., containing a
  `$Chains` component), or a list of multiple objects of class
  `demonoid`, where the number of components in the list is the number
  of chains.

- confidence:

  the coverage probability of the confidence interval for the potential
  scale reduction factor (PSRF).

- transform:

  logical. If `TRUE`, then marginal posterior distributions in `x` may
  be transformed to improve the normality of the distribution, which is
  assumed. A log-transform is applied to marginal posterior
  distributions in the interval \\(0, \infty\]\\, or a logit-transform
  is applied to marginal posterior distributions in the interval
  \\(0,1)\\.

## Value

A list with the following components:

- PSRF:

  a list containing the point-estimates of the potential scale reduction
  factor (labelled `Point Est.`) and the associated upper confidence
  limits (labelled `Upper C.I.`).

- MPSRF:

  the point-estimate of the multivariate potential scale reduction
  factor.

## Details

To use the `Gelman.Diagnostic` function, the user must first have
multiple MCMC chains for the same model, and three chains is usually
sufficient. The easiest way to obtain multiple chains is with the
`Chains` argument of the
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
function.

Although the
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
function does not simultaneously update multiple MCMC chains, it is easy
enough to obtain multiple chains, and if the computer has multiple
processors (which is common), then multiple chains may be obtained
simultaneously as follows. The model file may be opened in separate,
concurrent R sessions, and it is recommended that a maximum number of
sessions is equal to the number of processors, minus one. Each session
constitutes its own chain, and the code is identical, except the initial
values should be randomized with the
[`GIV`](https://robustecologies.github.io/lucifer/reference/GIV.md)
function so the chains begin in different places. The resulting object
of class `demonoid` for each chain is saved, all objects are read into
one session, put into a list, and passed to the `Gelman.Diagnostic`
function.

Initial values must be overdispersed with respect to each target
distribution, though these distributions are unknown in the beginning.
Since the `Gelman.Diagnostic` function relies heavily on overdispersion
with respect to the target distribution, the user should consider using
MCMC twice, first to estimate the target distributions, and secondly to
overdisperse initial values with respect to them. This may help identify
multimodal target distributions. If multiple modes are found, it remains
possible that more modes exist. When multiple modes are found, and if
chains are combined with the
[`Combine`](https://robustecologies.github.io/lucifer/reference/Combine.md)
function, each mode is probably not represented in a proportion correct
to the distribution.

The 'potential scale reduction factor' (PSRF) is an estimated factor by
which the scale of the current distribution for the target distribution
might be reduced if the simulations were continued for an infinite
number of iterations. Each PSRF declines to 1 as the number of
iterations approaches infinity. PSRF is also often represented as R-hat.
PSRF is calculated for each marginal posterior distribution in `x`,
together with upper and lower confidence limits. Approximate convergence
is diagnosed when the upper limit is close to 1. The recommended
proximity of each PSRF to 1 varies with each problem, but a general goal
is to achieve PSRF \< 1.1. PSRF is an estimate of how much narrower the
posterior might become with an infinite number of iterations. When PSRF
= 1.1, for example, it may be interpreted as a potential reduction of
10% in posterior interval width, given infinite iterations. The
multivariate form bounds above the potential scale reduction factor for
any linear combination of the (possibly transformed) variables.

The confidence limits are based on the assumption that the target
distribution is stationary and normally distributed. The `transform`
argument may be used to improve the normal approximation.

A large PSRF indicates that the between-chain variance is substantially
greater than the within-chain variance, so that longer simulation is
needed. If a PSRF is close to 1, then the associated chains are likely
to have converged to one target distribution. A large PSRF (perhaps
generally when a PSRF \> 1.2) indicates convergence failure, and can
indicate the presence of a multimodal marginal posterior distribution in
which different chains may have converged to different local modes (see
[`is.multimodal`](https://robustecologies.github.io/lucifer/reference/Mode.md)),
or the need to update the associated chains longer, because burn-in (see
[`burnin`](https://robustecologies.github.io/lucifer/reference/burnin.md))
has yet to be completed.

The `Gelman.Diagnostic` is essentially the same as the `gelman.diag`
function in the `coda` package, but here it is programmed to work with
objects of class `demonoid`.

There are two ways to estimate the variance of the stationary
distribution: the mean of the empirical variance within each chain,
\\W\\, and the empirical variance from all chains combined, which can be
expressed as

\$\$\widehat{\sigma}^2 = \frac{(n-1) W}{n} + \frac{B}{n}\$\$

where \\n\\ is the number of iterations and \\B/n\\ is the empirical
between-chain variance.

If the chains have converged, then both estimates are unbiased.
Otherwise the first method will *underestimate* the variance, since the
individual chains have not had time to range all over the stationary
distribution, and the second method will *overestimate* the variance,
since the initial values were chosen to be overdispersed (and this
assumes the target distribution is known, see above).

This convergence diagnostic is based on the assumption that each target
distribution is normal. A Bayesian probability interval (see
[`p.interval`](https://robustecologies.github.io/lucifer/reference/p.interval.md))
can be constructed using a t-distribution with mean

\$\$\widehat{\mu}=\mbox{Sample mean of all chains combined,}\$\$

variance

\$\$\widehat{V} = \widehat{\sigma}^2 + \frac{B}{mn},\$\$

and degrees of freedom estimated by the method of moments

\$\$d = \frac{2\widehat{V}^2}{\mbox{Var}(\widehat{V})}\$\$

Use of the t-distribution accounts for the fact that the mean and
variance of the posterior distribution are estimated. The convergence
diagnostic itself is

\$\$R=\sqrt{\frac{(d+3) \widehat{V}}{(d+1)W}}\$\$

Values substantially above 1 indicate lack of convergence. If the chains
have not converged, then Bayesian probability intervals based on the
t-distribution are too wide, and have the potential to shrink by this
factor if the MCMC run is continued.

The multivariate version of Gelman and Rubin's diagnostic was proposed
by Brooks and Gelman (1998). Unlike the univariate proportional scale
reduction factor, the multivariate version does not include an
adjustment for the estimated number of degrees of freedom.

## References

Brooks, S.P. and Gelman, A. (1998). "General Methods for Monitoring
Convergence of Iterative Simulations". *Journal of Computational and
Graphical Statistics*, 7, p. 434–455.

Gelman, A. and Rubin, D.B. (1992). "Inference from Iterative Simulation
using Multiple Sequences". *Statistical Science*, 7, p. 457–511.

## See also

[`Combine`](https://robustecologies.github.io/lucifer/reference/Combine.md),
[`GIV`](https://robustecologies.github.io/lucifer/reference/GIV.md),
[`is.multimodal`](https://robustecologies.github.io/lucifer/reference/Mode.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`MCSE`](https://robustecologies.github.io/lucifer/reference/MCSE.md),
and
[`p.interval`](https://robustecologies.github.io/lucifer/reference/p.interval.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## After updating with multiple chains via lucifer(..., Chains = 4):
# Gelman.Diagnostic(Fit)
} # }
```
