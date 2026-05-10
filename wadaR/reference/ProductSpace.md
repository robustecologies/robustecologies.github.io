# Product space method for Bayesian model selection

Computes Bayes factors and posterior model probabilities by constructing
a product-space supermodel that encompasses \\K\\ candidate models
simultaneously. A discrete model indicator \\M \sim
\mathrm{Categorical}(p_1, \ldots, p_K)\\ selects the active model at
each MCMC iteration. Parameters of inactive models are drawn from
pseudoprior distributions that approximate their posteriors. Posterior
model probabilities are estimated by the proportion of MCMC iterations
spent in each model, and pairwise Bayes factors follow from the
posterior-to-prior odds ratio.

## Usage

``` r
ProductSpace(
  Models,
  Data,
  Initial.Values = NULL,
  Prior.Model.Probabilities = NULL,
  Pseudopriors = NULL,
  Pilot.Iterations = 5000,
  Pilot.Algorithm = "CHARM",
  Iterations = 50000,
  Thinning = 10,
  Chains = 3L,
  CPUs = 3L,
  Bisection = TRUE,
  Bisection.Tolerance = 0.1,
  Bisection.Max.Rounds = 10,
  Bisection.Iterations = 5000,
  Production.Algorithm = "RAM",
  Slice.Width = 1,
  Status = 1000,
  Verbose = TRUE
)
```

## Arguments

- Models:

  A named list of Model functions, each following the standard lucifer
  contract: `function(parm, Data)` returning
  `list(LP, Dev, Monitor, yhat, parm)`. Models may have different
  numbers of parameters, different likelihoods, and different priors;
  they need not be nested.

- Data:

  Either a single Data list shared across all models, or a named list of
  Data lists (one per model). Each Data list must contain at least
  `parm.names` and `mon.names`.

- Initial.Values:

  A named list of numeric vectors (one per model), or `NULL` to use
  zeros as starting values. Each vector must match the length of
  `parm.names` in the corresponding Data list.

- Prior.Model.Probabilities:

  Numeric vector of length \\K\\ summing to 1, specifying prior
  probabilities for each model. If `NULL` (default), a uniform prior
  \\1/K\\ is used.

- Pseudopriors:

  A named list of lists, one per model. Each inner list contains `mean`
  and `sd` numeric vectors matching the model's parameter count. If
  `NULL` (default), pseudopriors are estimated automatically from pilot
  MCMC runs.

- Pilot.Iterations:

  Integer. Number of MCMC iterations for each pilot fit used to estimate
  pseudopriors. Default 5000.

- Pilot.Algorithm:

  Character. MCMC algorithm for pilot fits. Default `"CHARM"`.

- Iterations:

  Integer. Number of product-space MCMC iterations for the production
  run. Default 50000.

- Thinning:

  Integer. Thinning interval for the production MCMC. Default 10.

- Chains:

  Integer. Number of independent product-space chains. Default 3.

- CPUs:

  Integer. Number of CPUs for parallel execution. Default 3.

- Bisection:

  Logical. Whether to run the bisection algorithm (Lodewyckx et
  al., 2011) to calibrate prior model probabilities for balanced model
  activation. Default `TRUE`.

- Bisection.Tolerance:

  Numeric. Convergence tolerance for the bisection algorithm. The
  algorithm stops when the maximum absolute difference between observed
  and target posterior activation probabilities is below this value.
  Default 0.10.

- Bisection.Max.Rounds:

  Integer. Maximum number of bisection iterations. Default 10.

- Bisection.Iterations:

  Integer. Number of MCMC iterations per bisection round. Default 5000.

- Production.Algorithm:

  Character. MCMC algorithm for the production run. `"RAM"` (default)
  uses a hybrid sampler with Robust Adaptive Metropolis (Vihola, 2012)
  for active model parameters, exact enumeration for the discrete model
  indicator, and direct sampling from pseudopriors for inactive
  parameters. This reduces the cost per iteration from \\\sim 10 \sum
  d_k + K\\ model evaluations (slice sampling) to \\K + 1\\ (enumeration
  plus one block proposal). `"Gibbs"` uses the original auto-FC Gibbs
  sampler with univariate slice sampling.

- Slice.Width:

  Numeric scalar or vector. Slice width(s) for continuous parameters in
  the product-space Gibbs sampler (used only when
  `Production.Algorithm = "Gibbs"`). If scalar, applied to all
  parameters. Default 1.0.

- Status:

  Integer. Progress reporting interval. Default 1000.

- Verbose:

  Logical. Whether to print progress messages. Default `TRUE`.

## Value

An object of class `product_space`, a list containing:

- Posterior.Model.Probabilities:

  Named numeric vector of corrected posterior model probabilities.

- Bayes.Factors:

  Matrix of pairwise Bayes factors (row vs. column).

- Log.Bayes.Factors:

  Matrix of log Bayes factors.

- Per.Model.Posteriors:

  Named list of matrices, one per model, containing posterior samples
  conditional on that model being active.

- Per.Model.Summaries:

  Named list of summary matrices (mean, SD, 2.5%, 50%, 97.5%) for each
  model's parameters.

- Model.Indicator:

  Integer vector of model indicator values across all retained
  iterations.

- Transition.Matrix:

  K x K matrix of model-switch transition probabilities.

- Switching.Rate:

  Proportion of consecutive iterations where the model indicator
  changed.

- M.ESS:

  Effective sample size of the model indicator trace.

- Complexity:

  A list of posterior model space complexity metrics: `Shannon.Entropy`
  (raw, in nats), `Shannon.Normalized` (scaled to \\\[0,1\]\\),
  `Effective.K` (exp of entropy, Hill number of order 1),
  `Simpson.Diversity` (\\1 - \sum p_k^2\\), and `Gini.Effective.K` (Hill
  number of order 2).

- Prior.Model.Probabilities:

  The (possibly calibrated) prior model probabilities used in the
  production run.

- Pseudopriors:

  The pseudoprior specifications used.

- Pilot.Fits:

  List of pilot demonoid objects (one per model).

- Bisection.History:

  Data frame of bisection iterations, or NULL.

- MCMC.Fit:

  The raw demonoid object from the product-space MCMC.

- Production.Algorithm:

  Character; the algorithm used (`"RAM"` or `"Gibbs"`).

- Acceptance.Rates:

  (RAM only) Named numeric vector of per-model RAM acceptance rates.
  When multi-chain, averaged across chains.

- Adaptation.States:

  (RAM only) List of final per-model RAM adaptation states (proposal
  covariance, iteration counts).

- Chains.n:

  (Multi-chain RAM only) Integer number of chains, or `NULL` for
  single-chain runs.

- Per.Chain.Trans:

  (Multi-chain RAM only) List of per-chain transition matrices.

- Per.Chain.Switch:

  (Multi-chain RAM only) Numeric vector of per-chain switching rates.

- Per.Chain.M.ESS:

  (Multi-chain RAM only) Numeric vector of per-chain model indicator ESS
  values.

- Chain.Model.Probs:

  (Multi-chain RAM only) Matrix of per-chain raw model proportions (rows
  = chains, columns = models).

- Chain.Model.Probs.SD:

  (Multi-chain RAM only) Named numeric vector of cross-chain standard
  deviations of model proportions.

- M.Rhat:

  (Multi-chain RAM only) Split R-hat on the model indicator across
  chains.

- Model.Names:

  Character vector of model names.

- Param.Map:

  List mapping model names to super-parameter indices.

- K:

  Number of models.

- Minutes:

  Computation time in minutes.

- call:

  The matched call.

## Details

The product space method (Carlin and Chib, 1995) constructs a single
augmented parameter space \\\boldsymbol{\theta} = (\theta_1, \ldots,
\theta_K, M)\\ by concatenating all model-specific parameter vectors
with a discrete model indicator. When model \\M = k\\ is active, the
joint density is

\$\$p(\theta, M{=}k \mid y) \propto p(y \mid \theta_k, M{=}k)\\
p(\theta_k \mid M{=}k)\\ \prod\_{j \neq k} q(\theta_j \mid M{=}k)\\
p(M{=}k)\$\$

where \\q(\theta_j \mid M{=}k)\\ is the pseudoprior for model \\j\\'s
parameters when model \\k\\ is active. The pseudopriors integrate out to
1 and do not affect the marginal likelihood, but their choice critically
affects MCMC mixing: pseudopriors that approximate the true
model-specific posteriors produce frequent model switches and reliable
Bayes factor estimates.

The posterior model probability is estimated as \\\hat{p}(M{=}k \mid y)
= n_k / N\\ where \\n_k\\ is the number of iterations with \\M = k\\.
When non-uniform prior model probabilities are used (e.g., from
bisection calibration), the corrected posterior probabilities are

\$\$\tilde{p}(M{=}k \mid y) = \frac{\hat{p}(M{=}k \mid y) /
p(M{=}k)}{\sum_j \hat{p}(M{=}j \mid y) / p(M{=}j)}\$\$

and the Bayes factor between models \\a\\ and \\b\\ is \\B\_{ab} =
\tilde{p}(M{=}a \mid y) / \tilde{p}(M{=}b \mid y)\\.

The function proceeds in stages: (1) validate inputs; (2) run pilot MCMC
fits for each model in parallel to estimate pseudopriors; (3) construct
the product-space supermodel; (4) optionally run bisection calibration
to equalize posterior model activation; (5) run production MCMC; (6)
post-process to extract model probabilities, Bayes factors, per-model
posteriors, and diagnostics.

The default production algorithm (`"RAM"`) exploits the conditional
independence structure of the product space. When model \\M = k\\ is
active, inactive parameters \\\theta_j\\ (\\j \neq k\\) are
conditionally independent of the data and their full conditional is the
pseudoprior, so they are sampled directly with zero model evaluations.
The active model parameters are updated with a block Robust Adaptive
Metropolis (RAM) proposal (Vihola, 2012), and the discrete model
indicator is sampled via exact enumeration. This reduces the cost per
iteration from approximately \\10 \sum d_k + K\\ model evaluations
(univariate slice sampling) to \\K + 1\\, and eliminates the
correlation-blind univariate bottleneck for active model parameters.

When `Chains > 1` and `Production.Algorithm = "RAM"`, the hybrid sampler
runs multiple independent chains with jittered initial values and
independent adaptation states. Chains are dispatched in parallel via
`callr` (if available) or run sequentially. The combined posterior is
formed by row-binding all chain posteriors, and multi-chain convergence
diagnostics (split R-hat on M, cross-chain model probability SD,
per-chain transition matrices and ESS) are computed to assess
reliability. Bisection calibration also uses the RAM sampler when
`Production.Algorithm = "RAM"`, reducing calibration time
proportionally.

## References

Carlin, B.P. and Chib, S. (1995). "Bayesian model choice via Markov
chain Monte Carlo methods." *Journal of the Royal Statistical Society:
Series B*, 57(3), p. 473–484.

Lodewyckx, T., Kim, W., Lee, M.D., Tuerlinckx, F., Kuppens, P., and
Wagenmakers, E.-J. (2011). "A tutorial on Bayes factor estimation with
the product space method." *Journal of Mathematical Psychology*, 55(5),
p. 331–347.
[doi:10.1016/j.jmp.2011.06.001](https://doi.org/10.1016/j.jmp.2011.06.001)

Tenan, S., O'Hara, R.B., Hendriks, I., and Tavecchia, G. (2014).
"Bayesian model selection: the steepest mountain to climb." *Ecological
Modelling*, 283, p. 62–69.
[doi:10.1016/j.ecolmodel.2014.03.017](https://doi.org/10.1016/j.ecolmodel.2014.03.017)

Vihola, M. (2012). "Robust adaptive Metropolis algorithm with coerced
acceptance rate." *Statistics and Computing*, 22(5), p. 997–1008.
[doi:10.1007/s11222-011-9269-5](https://doi.org/10.1007/s11222-011-9269-5)

## See also

[`BayesFactor`](https://robustecologies.github.io/lucifer/reference/BayesFactor.md)
for Bayes factor computation from fitted objects,
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
for single-model MCMC sampling,
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md)
and
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md)
for deterministic approximations,
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md) and
[`SMC`](https://robustecologies.github.io/lucifer/reference/SMC.md) for
particle-based inference,
[`WAIC`](https://robustecologies.github.io/lucifer/reference/WAIC.md)
and [`LOO`](https://robustecologies.github.io/lucifer/reference/LOO.md)
for information-criteria model comparison,
[`print.product_space`](https://robustecologies.github.io/lucifer/reference/print.product_space.md),
[`summary.product_space`](https://robustecologies.github.io/lucifer/reference/summary.product_space.md),
[`plot.product_space`](https://robustecologies.github.io/lucifer/reference/plot.product_space.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(lucifer)

# Compare two models for normal data:
# M1: known mean (mu = 0), estimate sigma
# M2: estimate both mu and sigma
# Data generated under M1 (true mean = 0), so M1 should be preferred
# by the parsimony principle (Occam's razor via Bayes factors).

set.seed(42)
y <- rnorm(30, mean = 0, sd = 2)

Model1 <- function(parm, Data) {
    sigma <- interval(parm[1], 1e-100, Inf)
    parm[1] <- sigma
    LL <- sum(dnorm(Data$y, 0, sigma, log = TRUE))
    LP <- LL + dhalfcauchy(sigma, 5, log = TRUE)
    list(LP = LP, Dev = -2 * LL, Monitor = LP,
         yhat = rnorm(Data$n, 0, sigma), parm = parm)
}

Model2 <- function(parm, Data) {
    mu <- parm[1]
    sigma <- interval(parm[2], 1e-100, Inf)
    parm[2] <- sigma
    LL <- sum(dnorm(Data$y, mu, sigma, log = TRUE))
    LP <- LL + dnorm(mu, 0, 10, log = TRUE) +
          dhalfcauchy(sigma, 5, log = TRUE)
    list(LP = LP, Dev = -2 * LL, Monitor = LP,
         yhat = rnorm(Data$n, mu, sigma), parm = parm)
}

Data1 <- list(y = y, n = length(y),
              parm.names = "sigma", mon.names = "LP")
Data2 <- list(y = y, n = length(y),
              parm.names = c("mu", "sigma"), mon.names = "LP")

ps <- ProductSpace(
    Models = list(M1 = Model1, M2 = Model2),
    Data = list(M1 = Data1, M2 = Data2),
    Iterations = 20000, Pilot.Iterations = 3000,
    Chains = 2, CPUs = 2
)

print(ps)
summary(ps)
plot(ps)
plot(ps, type = "trace")
plot(ps, type = "bayes_factors")
plot(ps, type = "complexity")
} # }
```
