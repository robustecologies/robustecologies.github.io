# Comprehensive posterior sensitivity analysis

Performs five complementary sensitivity analyses on Bayesian posterior
distributions: power-scaling sensitivity, cross-model divergence,
prior-data conflict diagnostics, observation influence analysis, and
robust posterior aggregation. Accepts any lucifer fit object or a named
list of fit objects and returns a unified `robust_bayes` S3 object with
print, summary, and plot methods.

## Usage

``` r
RobustBayes(
  x,
  Model = NULL,
  Data = NULL,
  modules = "all",
  alpha = seq(0.5, 1.5, by = 0.05),
  sensitivity.threshold = 0.05,
  divergence = c("js", "wasserstein", "overlap"),
  log_lik = NULL,
  aggregate.method = c("stacking", "bma", "pseudobma"),
  verbose = TRUE
)
```

## Arguments

- x:

  A fit object of class `demonoid`, `laplace`, `iterquad`, `pmc`, `vb`,
  `smc`, `abc`, `sbi`, or `bayesquad`, or a named list of such objects.
  A single fit enables power-scaling, conflict, and influence analyses.
  Two or more fits additionally enable cross-model comparison and
  aggregation.

- Model:

  The Model function used for fitting, following the lucifer interface
  (returns a list with `LP`, `Dev`, `Monitor`, `yhat`, `parm`). Required
  for the `"power"`, `"conflict"`, and `"influence"` modules.

- Data:

  The Data list passed to `Model`. Required when `Model` is supplied.

- modules:

  Character vector specifying which analyses to run. Options are
  `"power"`, `"compare"`, `"conflict"`, `"influence"`, `"aggregate"`,
  and `"all"` (default). Modules whose prerequisites are not met are
  skipped with a message.

- alpha:

  Numeric vector of perturbation exponents for power-scaling
  sensitivity. Values \< 1 weaken the component, values \> 1 strengthen
  it. Default `seq(0.5, 1.5, by = 0.05)`.

- sensitivity.threshold:

  Numeric threshold for flagging parameters as sensitive in the
  power-scaling module. Default 0.05.

- divergence:

  Character string specifying the divergence measure for cross-model
  comparison. One of `"js"` (Jensen-Shannon, default), `"wasserstein"`
  (Wasserstein-2), or `"overlap"` (overlap coefficient).

- log_lik:

  Optional \\N \times S\\ matrix of pointwise log-likelihoods for the
  influence module. If `NULL` and `Model`/`Data` are supplied,
  recomputed from the Model.

- aggregate.method:

  Character string specifying the primary aggregation method. One of
  `"stacking"` (default), `"bma"`, or `"pseudobma"`.

- verbose:

  Logical; if `TRUE` (default), prints progress messages.

## Value

An object of class `robust_bayes`, which is a list with components:

- call:

  The matched function call.

- n_models:

  Number of models analyzed.

- n_params:

  Number of parameters (from the first model).

- n_obs:

  Number of observations, if determinable.

- param_names:

  Character vector of parameter names.

- model_names:

  Character vector of model names.

- fit_classes:

  Character vector of fit object classes.

- power:

  Power-scaling results (including stored posterior samples and
  log-prior/log-likelihood decomposition for density plotting via
  `plot(x, type = "power_density")`), or `NULL`.

- compare:

  Cross-model comparison results including marginal divergence array,
  joint energy distance matrix, stochastic dominance probabilities, and
  per-parameter sensitivity ranking, or `NULL`.

- conflict:

  Prior-data conflict results, or `NULL`.

- influence:

  Observation influence results, or `NULL`.

- aggregate:

  Robust aggregation results, or `NULL`.

- minutes:

  Elapsed computation time in minutes.

## Details

`RobustBayes()` exploits a key property of the lucifer Model interface:
every Model evaluation returns both `LP` (log-posterior) and `Dev` (\\-2
\times\\ log-likelihood), enabling automatic decomposition into
log-prior and log-likelihood without requiring the user to separate
these components: \$\$\log\pi(\theta) = \mathrm{LP}(\theta) +
\mathrm{Dev}(\theta)/2\$\$ \$\$\ell(\theta) =
-\mathrm{Dev}(\theta)/2\$\$

The five modules address distinct sources of inferential fragility:

**Power-scaling sensitivity** (Kallioinen et al. 2024) evaluates
sensitivity to the prior and likelihood by computing the power-scaled
posterior \\\tilde{p}(\theta\|y,\alpha) \propto \pi(\theta)^\alpha
L(y\|\theta)\\ via PSIS importance reweighting, without refitting the
model. For each perturbation exponent \\\alpha\\ in the grid specified
by the `alpha` argument, the function computes importance weights
\\w_s(\alpha) = \pi(\theta_s)^{\alpha - 1}\\, smooths them via
Pareto-smoothed importance sampling, and estimates the Jensen-Shannon
divergence between the original and perturbed posterior marginals. The
cumulative sensitivity score for each parameter is obtained by
trapezoidal integration of the JS divergence curve over \\\alpha\\.
Parameters exceeding the `sensitivity.threshold` are classified as
prior-sensitive, likelihood-sensitive, or conflicted depending on which
perturbation direction produces larger divergence.

**Cross-model comparison** quantifies posterior divergence across a
family of \\K \geq 2\\ models using three complementary metrics on the
marginal posteriors: Jensen-Shannon divergence (a symmetric, bounded
measure of distributional discrepancy), Wasserstein-2 distance (the
optimal transport cost between distributions, sensitive to location
shifts), and the overlap coefficient (the area of intersection between
two density estimates, ranging from 0 for disjoint distributions to 1
for identical ones). Additionally, stochastic dominance probabilities
\\P(\theta_i \> \theta_j)\\ are computed for every model pair and every
shared parameter, giving the directional probability that one model's
posterior exceeds the other's. These complement the symmetric divergence
measures by indicating the sign of the shift. Joint multivariate
dissimilarity is assessed via the energy distance, a kernel-based
two-sample statistic that captures differences in both location and
shape without assuming a parametric form.

**Prior-data conflict** (Box 1980; Evans and Moshonov 2006) assesses
whether the observed data are surprising under the prior predictive
distribution. A prior generative function (`Data$PGF`) is required. The
module draws \\B\\ parameter vectors from the prior, evaluates the model
at each draw, and computes a Box p-value as the proportion of
prior-predictive deviances exceeding the observed deviance.
Per-parameter Evans-Moshonov calibration compares each posterior
marginal against its prior predictive counterpart. Small p-values signal
genuine conflict between the prior and the data, suggesting the prior is
too concentrated in a region that the data disfavor.

**Observation influence** uses PSIS-LOO importance reweighting to
quantify each observation's impact on the posterior without refitting.
For each observation \\i\\, leave-one-out importance weights are
computed from the pointwise log-likelihoods, and the shift in each
posterior marginal moment is estimated. Observations whose Pareto-k
diagnostic exceeds 0.7 are flagged as high-leverage points that may
disproportionately drive the posterior; k \> 1.0 indicates the
importance sampling approximation is unreliable and refitting is
advised.

**Robust aggregation** (Yao et al. 2018) combines multiple models into a
single predictive distribution using model weights. Three methods are
available: Bayesian model averaging (BMA), which weights by marginal
likelihoods; stacking (Wolpert 1992; Yao et al. 2018), which optimizes
pointwise predictive performance via leave-one-out cross-validation; and
pseudo-BMA+ (Yao et al. 2018), which applies a Bayesian bootstrap to
stabilize LOO-based BMA weights. Stacking is the default because it
targets predictive accuracy rather than posterior probability of model
correctness, and is therefore more robust to M-open settings where no
candidate model is the true data-generating process.

## References

Kallioinen, N., Paananen, T., Buerkner, P.-C., and Vehtari, A. (2024).
"Detecting and diagnosing prior and likelihood sensitivity with
power-scaling." *Statistics and Computing*, 34, 57.
[doi:10.1007/s11222-023-10366-5](https://doi.org/10.1007/s11222-023-10366-5)

Yao, Y., Vehtari, A., Simpson, D., and Gelman, A. (2018). "Using
stacking to average Bayesian predictive distributions." *Bayesian
Analysis*, 13(3), 917–1007.
[doi:10.1214/17-BA1091](https://doi.org/10.1214/17-BA1091)

Box, G.E.P. (1980). "Sampling and Bayes' inference in scientific
modelling and robustness." *Journal of the Royal Statistical Society,
Series A*, 143(4), 383–430.
[doi:10.2307/2982063](https://doi.org/10.2307/2982063)

Evans, M. and Moshonov, H. (2006). "Checking for prior-data conflict."
*Bayesian Analysis*, 1(4), 893–914.
[doi:10.1214/06-BA129](https://doi.org/10.1214/06-BA129)

Vehtari, A., Simpson, D., Gelman, A., Yao, Y., and Gabry, J. (2024).
"Pareto smoothed importance sampling." *Journal of Machine Learning
Research*, 25(72), 1–58.

Berger, J.O. (1994). "An overview of robust Bayesian analysis." *Test*,
3, 5–124. [doi:10.1007/BF02562676](https://doi.org/10.1007/BF02562676)

## See also

[`LOO`](https://robustecologies.github.io/lucifer/reference/LOO.md),
[`PSIS`](https://robustecologies.github.io/lucifer/reference/PSIS.md),
[`BayesFactor`](https://robustecologies.github.io/lucifer/reference/BayesFactor.md),
[`Arena`](https://robustecologies.github.io/lucifer/reference/Arena.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate data for a simple linear regression
set.seed(42)
N <- 30
x <- rnorm(N)
y <- 0.5 + 1.2 * x + rnorm(N, 0, 0.8)

# Data list with prior generative function for conflict module
Data <- list(
  mon.names = "LP",
  parm.names = c("beta0", "beta1", "log.sigma"),
  y = y, x = x, N = N,
  PGF = function(Data) {
    c(rnorm(1, 0, 10), rnorm(1, 0, 10), log(abs(rcauchy(1, 0, 5))))
  }
)

# Model with wide N(0,10) priors
Model <- function(parm, Data) {
  beta0 <- parm[1]; beta1 <- parm[2]; sigma <- exp(parm[3])
  mu <- beta0 + beta1 * Data$x
  LL <- sum(dnorm(Data$y, mu, sigma, log = TRUE))
  LP <- LL + sum(dnorm(parm[1:2], 0, 10, log = TRUE)) +
             dhalfcauchy(sigma, 5, log = TRUE)
  yhat <- rnorm(Data$N, mu, sigma)
  list(LP = LP, Dev = -2 * LL, Monitor = LP, yhat = yhat, parm = parm)
}

# Fit with MCMC
fit <- lucifer(Model, Data, Initial.Values = c(0, 0, 0),
               Algorithm = "NUTS", Iterations = 5000, Status = 1000)

# Single-model sensitivity (power-scaling, conflict, influence)
rb <- RobustBayes(fit, Model = Model, Data = Data)
print(rb)
summary(rb)
plot(rb)
plot(rb, type = "power")
plot(rb, type = "influence")

# Multi-model: compare wide vs tight priors
Model_tight <- function(parm, Data) {
  beta0 <- parm[1]; beta1 <- parm[2]; sigma <- exp(parm[3])
  mu <- beta0 + beta1 * Data$x
  LL <- sum(dnorm(Data$y, mu, sigma, log = TRUE))
  LP <- LL + sum(dnorm(parm[1:2], 0, 1, log = TRUE)) +
             dhalfcauchy(sigma, 5, log = TRUE)
  yhat <- rnorm(Data$N, mu, sigma)
  list(LP = LP, Dev = -2 * LL, Monitor = LP, yhat = yhat, parm = parm)
}
fit2 <- lucifer(Model_tight, Data, Initial.Values = c(0, 0, 0),
                Algorithm = "NUTS", Iterations = 5000, Status = 1000)

rb2 <- RobustBayes(list(wide = fit, tight = fit2),
                    Model = Model, Data = Data)
plot(rb2, type = "divergence")
plot(rb2, type = "weights")
} # }
```
