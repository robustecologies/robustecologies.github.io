# Bayesian quadrature for posterior integration

Performs Bayesian quadrature using a Gaussian process surrogate to model
the integrand and compute a posterior distribution over the integral
value, rather than a mere point estimate. Active learning selects
evaluation points that maximally reduce integral uncertainty.
`BayesianQuadrature` implements four modern Bayesian quadrature methods
that model the integrand with a Gaussian process and produce
uncertainty-calibrated integral estimates.

## Usage

``` r
BayesianQuadrature(
  Model,
  parm,
  Data,
  Covar = NULL,
  Iterations = 100,
  Algorithm = "BQ",
  Specs = NULL,
  Samples = 1000,
  sir = TRUE,
  Stop.Tolerance = 0.001,
  CPUs = 1,
  Type = "PSOCK"
)
```

## Arguments

- Model:

  A function taking parameters `parm` and `Data`, returning a list with
  components `LP` (log-posterior), `Dev` (deviance), `Monitor`, `yhat`,
  and `parm`.

- parm:

  A vector of initial parameter values, ideally near the posterior mode.

- Data:

  A list containing at least `mon.names` and `parm.names`, plus user
  data required by `Model`.

- Covar:

  An optional \\J \times J\\ covariance matrix specifying the Gaussian
  reference measure. Defaults to a scaled identity matrix.

- Iterations:

  Maximum number of active learning iterations (default 100).

- Algorithm:

  Character string specifying the BQ method: `"BQ"`, `"WSABI-L"`,
  `"FWBQ"`, or `"BatchBQ"`.

- Specs:

  An optional list of algorithm-specific parameters:

  BQ

  :   `list(n_initial = 2*d, refit_every = 5)`

  WSABI-L

  :   `list(n_initial = 2*d)`

  FWBQ

  :   `list(n_initial = 1)`

  BatchBQ

  :   `list(n_initial = 2*d, batch_size = CPUs)`

- Samples:

  Number of posterior samples to draw via SIR (default 1000).

- sir:

  Logical; draw posterior samples using sampling importance resampling
  (default `TRUE`).

- Stop.Tolerance:

  Convergence tolerance for the relative integral standard deviation
  (default 1e-3).

- CPUs:

  Number of CPUs for parallel model evaluation (default 1). Used by
  `BatchBQ` and for SIR sampling.

- Type:

  Cluster type for parallel: `"PSOCK"` or `"MPI"`.

## Value

An object of class `bayesquad` containing:

- Algorithm:

  Character string naming the algorithm used.

- Call:

  The matched function call.

- Converged:

  Logical; whether the integral variance met the stopping criterion.

- Covar:

  Final covariance matrix.

- Deviance:

  Vector of deviance values at each evaluation point.

- GP:

  List containing the final GP state: evaluation points `X`, function
  values `y`, kernel matrix `K`, hyperparameters `sigma_f` and
  `lengthscales`.

- History:

  Matrix of parameter values across iterations.

- Initial.Values:

  The initial parameter vector.

- Integral:

  List with `mean` (posterior mean of integral), `variance` (posterior
  variance), `log.mean` (log scale), and `CI` (95% credible interval).

- Iterations:

  Number of iterations completed.

- LML:

  Log marginal likelihood estimate.

- LP.Final:

  Final log-posterior value.

- LP.Initial:

  Initial log-posterior value.

- Minutes:

  Runtime in minutes.

- Monitor:

  Matrix of monitored values from posterior samples.

- N:

  Number of GP evaluation points.

- Posterior:

  Matrix of posterior samples (if `sir = TRUE` and converged).

- Summary1:

  Point-estimate summary (Mean, SD, LB, UB).

- Summary2:

  Posterior sample summary (if available).

- Tolerance.Final:

  Final relative integral standard deviation.

- Tolerance.Stop:

  The stopping tolerance.

## Details

Bayesian quadrature, introduced by O'Hagan (1991), treats numerical
integration as a statistical inference problem. The integrand is modeled
as a realization of a Gaussian process, and the integral becomes a
linear functional of that GP, inheriting a Gaussian posterior
distribution. This contrasts with classical quadrature (as in
[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md))
which produces only a point estimate.

The squared exponential kernel \$\$k(\theta, \theta') = \sigma_f^2
\exp\left(-\frac{1}{2} \sum\_{d=1}^D \frac{(\theta_d -
\theta'\_d)^2}{\ell_d^2}\right)\$\$ is used with automatic relevance
determination (ARD) lengthscales. The kernel mean embedding against the
Gaussian reference measure \\\pi(\theta) = \mathcal{N}(\mu_0,
\Sigma_0)\\ yields a closed-form posterior over the integral value.

Four algorithms are available:

**BQ** (vanilla Bayesian quadrature): the baseline method of Rasmussen
and Ghahramani (2003). Models \\f(\theta) = \exp(\text{LP}(\theta))\\
with a GP and selects nodes by maximizing integral variance reduction.

**WSABI-L** (warped sequential active Bayesian integration, linearized):
Gunter et al. (2014). Applies a square-root warping \\g = \sqrt{2(f -
\alpha)}\\ to enforce non-negativity of the likelihood integrand, then
linearizes the inverse warping for moment propagation. Recommended for
marginal likelihood estimation.

**FWBQ** (Frank-Wolfe Bayesian quadrature): Briol et al. (2015). Uses
Frank-Wolfe optimization for node selection and weight computation. Less
sensitive to GP hyperparameters with provable convergence guarantees
(exponential rate under regularity conditions).

**BatchBQ**: greedy batch selection of \\B\\ nodes for parallel model
evaluation. Each batch is selected to jointly maximize integral variance
reduction, enabling efficient use of the `CPUs` argument.

Bayesian quadrature is practical for problems with \\d \le 20\\
parameters where each model evaluation is moderately expensive. For
higher-dimensional problems, use
[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md)
with algorithm `"CAGH"`, or MCMC via
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

## References

O'Hagan, A. (1991). Bayes-Hermite quadrature. *Journal of the Royal
Statistical Society: Series B*, 53(1), 145-169.

Rasmussen, C.E. and Ghahramani, Z. (2003). Bayesian Monte Carlo.
*Advances in Neural Information Processing Systems*, 15.

Gunter, T., Osborne, M.A., Garnett, R., Hennig, P., and Roberts, S.J.
(2014). Sampling for inference in probabilistic models with fast
Bayesian quadrature. *Advances in Neural Information Processing
Systems*, 27.

Briol, F.-X., Oates, C.J., Girolami, M., and Osborne, M.A. (2015).
Frank-Wolfe Bayesian quadrature: probabilistic integration with
theoretical guarantees. *Advances in Neural Information Processing
Systems*, 28.

Briol, F.-X., Oates, C.J., Girolami, M., Osborne, M.A., and Sejdinovic,
D. (2019). Probabilistic integration: a role in statistical computation?
*Statistical Science*, 34(1), 1-22.
[doi:10.1214/18-STS660](https://doi.org/10.1214/18-STS660)

## See also

[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Normal-normal conjugate model
Model <- function(parm, Data) {
  mu <- parm[1]
  LL <- sum(dnorm(Data$y, mu, 1, log = TRUE))
  LP <- LL + dnorm(mu, 0, 1000, log = TRUE)
  list(LP = LP, Dev = -2 * LL,
       Monitor = LP, yhat = rep(mu, length(Data$y)),
       parm = parm)
}
set.seed(42)
y <- rnorm(50, mean = 3)
Data <- list(y = y, mon.names = "LP", parm.names = "mu")

# Vanilla BQ
fit <- BayesianQuadrature(Model, parm = 0, Data = Data,
  Iterations = 50, Algorithm = "BQ")
print(fit)
summary(fit)
plot(fit, Data = Data)

# WSABI-L for marginal likelihood
fit_wsabi <- BayesianQuadrature(Model, parm = 0, Data = Data,
  Iterations = 50, Algorithm = "WSABI-L")

# FWBQ with provable convergence
fit_fwbq <- BayesianQuadrature(Model, parm = 0, Data = Data,
  Iterations = 50, Algorithm = "FWBQ")
} # }
```
