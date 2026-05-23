# Sequential Monte Carlo sampler

Performs particle-based Bayesian inference by annealing from the prior
to the posterior through a sequence of tempered distributions. At each
stage particles are reweighted, optionally resampled when the effective
sample size drops below a threshold, and rejuvenated via an MCMC kernel.
The marginal likelihood estimate is available as a byproduct.

## Usage

``` r
SMC(
  Model,
  Data,
  Initial.Values = NULL,
  Covar = NULL,
  Iterations = 1000,
  N.particles = 500,
  ESS.threshold = 0.5,
  Rejuvenation = "RWM",
  Rejuvenation.steps = 5,
  Schedule = "adaptive",
  Status = 100
)
```

## Arguments

- Model:

  A function with signature `function(parm, Data)` that returns the
  standard five-component list (LP, Dev, Monitor, yhat, parm) used
  throughout lucifer.

- Data:

  A list containing at least `mon.names` and `parm.names`, as required
  by
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

- Initial.Values:

  A numeric vector of initial parameter values, or `NULL` to draw from a
  simple Gaussian prior.

- Covar:

  A covariance matrix for the rejuvenation kernel proposal. If `NULL`,
  the identity matrix scaled by 2.38^2/d is used.

- Iterations:

  Not used directly; included for interface consistency.

- N.particles:

  Integer number of particles (default 500).

- ESS.threshold:

  Scalar in (0, 1). Resampling is triggered when the normalized ESS
  falls below `ESS.threshold * N.particles`.

- Rejuvenation:

  Character string for the MCMC kernel. Currently only `"RWM"`
  (random-walk Metropolis) is supported.

- Rejuvenation.steps:

  Integer number of MCMC steps per rejuvenation (default 5).

- Schedule:

  Character or numeric. Character `"adaptive"` (ESS-based selection of
  the next temperature, default), character `"ibis"` (sequential data
  incorporation; requires `Data$batches`), or any other character
  (equally spaced fallback). A numeric vector of inverse temperatures in
  (0, 1\] is also accepted as a user-supplied fixed schedule, in which
  case tempering steps through the supplied betas in order.

- Status:

  Integer controlling how often progress is printed.

## Value

An object of class `smc`, a list containing:

- Posterior:

  Matrix of final weighted particle positions (N.particles x LIV).

- Weights:

  Normalized final importance weights.

- Log.Marginal.Likelihood:

  Estimate of log p(y), computed as the sum of log normalizing constants
  across tempering stages.

- Schedule:

  The temperature sequence used.

- ESS.history:

  ESS at each tempering stage.

- Acceptance.rate:

  Rejuvenation acceptance rate per stage.

- Summary:

  Matrix of posterior summaries (mean, SD, quantiles).

- call:

  The matched call.

## Details

Sequential Monte Carlo sampler

The algorithm constructs a sequence of tempered distributions
\\\pi_t(\theta) \propto p(\theta) p(y\|\theta)^{\beta_t}\\ where \\0 =
\beta_0 \< \beta_1 \< \ldots \< \beta_T = 1\\. When
`Schedule = "adaptive"`, each \\\beta\_{t+1}\\ is chosen so that the ESS
of the incremental weights equals `ESS.threshold * N.particles`, using
bisection search. This follows the adaptive tempering strategy of Del
Moral, Doucet, and Jasra (2006).

## References

Del Moral, P., Doucet, A., and Jasra, A. (2006). "Sequential Monte Carlo
samplers." *Journal of the Royal Statistical Society: Series B*, 68(3),
p. 411–436.
[doi:10.1111/j.1467-9868.2006.00553.x](https://doi.org/10.1111/j.1467-9868.2006.00553.x)

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Normal-Normal conjugate model
Model <- function(parm, Data) {
    mu <- parm[1]
    LL <- sum(dnorm(Data$y, mu, 1, log = TRUE))
    LP <- LL + dnorm(mu, 0, 10, log = TRUE)
    list(LP = LP, Dev = -2 * LL, Monitor = LP,
         yhat = rnorm(Data$n, mu, 1), parm = parm)
}
Data <- list(y = rnorm(50, 3, 1), n = 50,
             mon.names = "LP", parm.names = "mu")
fit <- SMC(Model, Data, Initial.Values = 0, N.particles = 200)
print(fit)
summary(fit)
plot(fit)
} # }
```
