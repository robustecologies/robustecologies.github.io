# Logarithm of the marginal likelihood

This function approximates the logarithm of the marginal likelihood
(LML), where the marginal likelihood (also called the integrated
likelihood or the prior predictive distribution) is \\p(\textbf{y}) =
\int p(\textbf{y} \| \Theta) p(\Theta) d\Theta\\. Several methods of
approximation are available.

## Usage

``` r
LML(
  Model = NULL,
  Data = NULL,
  Modes = NULL,
  theta = NULL,
  LL = NULL,
  Covar = NULL,
  method = "NSIS"
)
```

## Arguments

- Model:

  the model specification function. Used only with the `LME` method.

- Data:

  a list of data passed to the model specification. Used only with the
  `LME` method.

- Modes:

  a vector of posterior modes (or medians for MCMC). Used with the `GD`
  or `LME` methods.

- theta:

  a matrix of posterior samples (parameters only). Used with the `GD`,
  `HME`, or `NSIS` methods.

- LL:

  a vector of log-likelihood samples. Used with the `GD`, `HME`, or
  `NSIS` methods.

- Covar:

  the covariance matrix of the posterior modes. Used with the `GD` or
  `LME` methods.

- method:

  character; one of `"GD"` (Gelfand-Dey), `"HME"` (Harmonic Mean
  Estimator), `"LME"` (Laplace-Metropolis Estimator), or `"NSIS"`
  (nonparametric self-normalized importance sampling, the default).

## Value

A list with two components:

- LML:

  the approximation of the logarithm of the marginal likelihood, or `NA`
  if the method fails.

- VarCov:

  the variance-covariance matrix (negative inverse of the Hessian), if
  estimated; otherwise `NA`.

## Details

Generally a user of
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md), or
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md)
does not need to call `LML` directly, since these functions already
include it. However, `LML` may be called to estimate the marginal
likelihood with a different method or with non-stationary chains.

The `GD` method (Gelfand and Dey, 1994) is a stable modification of the
harmonic mean estimator. The `HME` method (Newton and Raftery, 1994) is
fast but unreliable and should be avoided. The `LME` method uses the
Laplace-Metropolis Estimator and is inappropriate with hierarchical
models. The `NSIS` method uses nonparametric self-normalized importance
sampling and requires at least 301 stationary samples with no more
parameters than half the number of samples.

## References

Gelfand, A.E. and Dey, D.K. (1994). "Bayesian Model Choice: Asymptotics
and Exact Calculations". *Journal of the Royal Statistical Society*,
Series B 56, p. 501–514.

Lewis, S.M. and Raftery, A.E. (1997). "Estimating Bayes Factors via
Posterior Simulation with the Laplace-Metropolis Estimator". *Journal of
the American Statistical Association*, 92, p. 648–655.

Newton, M.A. and Raftery, A.E. (1994). "Approximate Bayesian Inference
by the Weighted Likelihood Bootstrap". *Journal of the Royal Statistical
Society*, Series B 3, p. 3–48.

## See also

[`BayesFactor`](https://robustecologies.github.io/lucifer/reference/BayesFactor.md),
[`is.proper`](https://robustecologies.github.io/lucifer/reference/is.proper.md),
[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Applying HME to an object of class demonoid:
LML(LL = Fit$Deviance * (-1/2), method = "HME")

# Applying LME to an object of class demonoid:
LML(Model, MyData,
    Modes = apply(Fit$Posterior1, 2, median), method = "LME")

# Applying NSIS to an object of class demonoid:
LML(theta = Fit$Posterior1, LL = Fit$Deviance * -(1/2),
    method = "NSIS")
} # }
```
