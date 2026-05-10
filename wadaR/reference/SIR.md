# Sampling importance resampling

The `SIR` function performs Sampling Importance Resampling, also called
Sequential Importance Resampling, and uses a multivariate normal
proposal density.

## Usage

``` r
SIR(Model, Data, mu, Sigma, n = 1000, CPUs = 1, Type = "PSOCK")
```

## Arguments

- Model:

  A model specification function. For more information, see
  [`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md).

- Data:

  A list of data. For more information, see
  [`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md).

- mu:

  A mean vector, \\\mu\\, for a multivariate normal distribution,
  usually the posterior means from an object of class `iterquad` (from
  [`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md))
  or class `vb` (from
  [`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md)),
  or the posterior modes from an object of class `laplace` (from
  [`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md)).

- Sigma:

  A covariance matrix, \\\Sigma\\, for a multivariate normal
  distribution, usually the `Covar` component of an object of class
  `iterquad`, `laplace`, or `vb`.

- n:

  The number of samples to be drawn from the posterior distribution.

- CPUs:

  An integer that specifies the number of central processing units
  (CPUs) of the multicore computer or computer cluster. This argument
  defaults to `CPUs=1`, in which parallel processing does not occur.

- Type:

  The type of parallel processing to perform, accepting either
  `Type="PSOCK"` or `Type="MPI"`.

## Value

A matrix of samples drawn from the posterior distribution.

## Details

Sampling Importance Resampling (SIR) was introduced in Gordon et al.
(1993), and is the original particle filtering algorithm (and this
family of algorithms is also known as Sequential Monte Carlo). A
distribution is approximated with importance weights, which are
approximations to the relative posterior densities of the particles, and
the sum of the weights is one. In this terminology, each sample in the
distribution is a "particle". SIR is a sequential or recursive form of
importance sampling. As in importance sampling, the expectation of a
function can be approximated as a weighted average. The optimal proposal
distribution is the target distribution.

In the `lucifer` package, the main use of the `SIR` function is to
produce posterior samples for iterative quadrature, Laplace
Approximation, or Variational Bayes, and `SIR` is called
behind-the-scenes by the
[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
or
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md)
function.

Iterative quadrature estimates the posterior mean and the associated
covariance matrix. Assuming normality, this output characterizes the
marginal posterior distributions. However, it is often useful to have
posterior samples, in which case the `SIR` function is used to draw
samples. The number of samples, `n`, should increase with the number and
intercorrelations of the parameters. Otherwise, multimodal posterior
distributions may occur.

Laplace Approximation estimates the posterior mode and the associated
covariance matrix. Assuming normality, this output characterizes the
marginal posterior distributions. However, it is often useful to have
posterior samples, in which case the `SIR` function is used to draw
samples. The number of samples, `n`, should increase with the number and
intercorrelations of the parameters. Otherwise, multimodal posterior
distributions may occur.

Variational Bayes estimates both the posterior mean and variance.
Assuming normality, this output characterizes the marginal posterior
distributions. However, it is often useful to have posterior samples, in
which case the `SIR` function is used to draw samples. The number of
samples, `n`, should increase with the number of intercorrelations of
the parameters. Otherwise, multimodal posterior distributions may occur.

SIR is also commonly used when considering a mild change in a prior
distribution. For example, suppose a model was updated in
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
and it had a least-informative prior distribution, but the statistician
would like to estimate the impact of changing to a weakly-informative
prior distribution. The change is made in the model specification
function, and the posterior means and covariance are supplied to the
`SIR` function. The returned samples are estimates of the posterior,
given the different prior distribution. This is akin to sensitivity
analysis (see
[`RobustBayes`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md)
for a comprehensive framework).

In other contexts (for which this function is not designed), SIR is used
with dynamic linear models (DLMs) and state-space models (SSMs) for
state filtering.

Parallel processing may be performed when the user specifies `CPUs` to
be greater than one, implying that the specified number of CPUs exists
and is available. Parallelization may be performed on a multicore
computer or a computer cluster. Either a Simple Network of Workstations
(SNOW) or Message Passing Interface (MPI) is used. With small data sets
and few samples, parallel processing may be slower, due to computer
network communication. With larger data sets and more samples, the user
should experience a faster run-time.

This function was adapted from the `sir` function in the `LearnBayes`
package.

## References

Gordon, N.J., Salmond, D.J., and Smith, A.F.M. (1993). "Novel Approach
to Nonlinear/Non-Gaussian Bayesian State Estimation". *IEEE Proceedings
F on Radar and Signal Processing*, 140(2), p. 107–113.

## See also

[`dmvn`](https://robustecologies.github.io/lucifer/reference/dist.Multivariate.Normal.md),
[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`RobustBayes`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md),
and
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving SIR
} # }
```
