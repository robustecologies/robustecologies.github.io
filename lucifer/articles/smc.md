# Sequential Monte Carlo samplers in lucifer

## Introduction

Markov chain Monte Carlo methods produce correlated samples from the
posterior by constructing an ergodic chain whose stationary distribution
is \\p(\theta \mid y)\\. Sequential Monte Carlo (SMC) samplers offer a
fundamentally different strategy: they propagate a population of
weighted particles through a sequence of intermediate distributions that
bridge the prior to the posterior, combining importance sampling with
resampling and MCMC rejuvenation to maintain particle diversity
[\[1\]](#ref1). This approach has several practical advantages over
single-chain MCMC. SMC provides an unbiased estimate of the marginal
likelihood \\p(y)\\ as a natural byproduct, which is essential for
Bayesian model comparison. It is inherently parallel, since each
particle can be evaluated independently. And it avoids the burn-in and
convergence diagnostics that complicate MCMC, because the particle
population is initialized from the prior and gradually guided toward the
posterior without requiring stationarity.

The
[`SMC()`](https://robustecologies.github.io/lucifer/reference/SMC.md)
function in **lucifer** implements the SMC sampler framework of Del
Moral, Doucet, and Jasra [\[1\]](#ref1) with adaptive tempering,
systematic resampling, and random-walk Metropolis rejuvenation. The
computationally intensive components, including incremental weight
evaluation, systematic resampling, and the rejuvenation kernel, are
implemented in C++ with Rcpp for efficiency.

## The tempering approach

### Tempered distributions

The core idea is to construct a sequence of tempered distributions that
interpolate between the prior and the posterior:

\\ \pi_t(\theta) \propto p(\theta) \\ p(y \mid \theta)^{\beta_t}, \qquad
0 = \beta_0 \< \beta_1 \< \cdots \< \beta_T = 1. \\

At \\\beta_0 = 0\\, the tempered distribution equals the prior
\\p(\theta)\\; at \\\beta_T = 1\\, it equals the full posterior
\\p(\theta \mid y)\\. Intermediate values of \\\beta\\ define
distributions that are “softened” versions of the posterior, where the
likelihood is raised to a fractional power. This tempering scheme
ensures that adjacent distributions \\\pi_t\\ and \\\pi\_{t+1}\\ are
sufficiently close that importance weights remain well-behaved, avoiding
the catastrophic weight degeneracy that would occur if one attempted to
importance-sample directly from the prior to the posterior.

The temperature schedule \\\\\beta_t\\\\ can be fixed in advance (e.g.,
equally spaced on \\\[0, 1\]\\) or chosen adaptively. The adaptive
strategy, following [\[1\]](#ref1), selects each \\\beta\_{t+1}\\ so
that the effective sample size (ESS) of the incremental weights does not
drop below a user-specified threshold. This is implemented via bisection
search: the algorithm finds the largest \\\beta\_{t+1}\\ such that the
ESS of the weights \\w_i \propto p(y \mid \theta_i)^{\beta\_{t+1} -
\beta_t}\\ remains above `ESS.threshold * N.particles`. Adaptive
tempering concentrates stages where the likelihood surface is most
informative, typically inserting many small \\\beta\\ increments early
in the sequence (where the likelihood first begins to dominate the
prior) and larger steps near the end.

### Marginal likelihood estimation

A key advantage of SMC over MCMC is that it produces an unbiased
estimate of the marginal likelihood. At each tempering stage \\t\\, the
incremental weights \\w_i^{(t)} = p(y \mid \theta_i)^{\beta_t -
\beta\_{t-1}}\\ have a normalizing constant that estimates the ratio
\\Z_t / Z\_{t-1}\\ of successive normalizing constants. The log-marginal
likelihood is then

\\ \log \hat{p}(y) = \sum\_{t=1}^{T} \log \left(\frac{1}{N}
\sum\_{i=1}^{N} w_i^{(t)}\right) \\

where the sum telescopes across all tempering stages. This estimate is
unbiased on the natural (non-log) scale [\[2\]](#ref2), which means it
can be used directly for Bayes factor computation between competing
models.

## The SMC algorithm

### Initialization

The algorithm begins by initializing \\N\\ particles from the prior. In
**lucifer**, particles are created by jittering the user-supplied
initial values with Gaussian noise: \\\theta_i = \theta_0 +
\varepsilon_i\\ where \\\varepsilon_i \sim \mathcal{N}(0,
\text{diag}((0.1 \cdot \|\theta_0\| + 0.01)^2))\\. Each particle is
validated through the model function to ensure it produces a finite
log-posterior; invalid particles are reset to the initial values. The
initial weights are uniform: \\w_i = 1/N\\.

### Reweighting

At each tempering stage, the temperature increases from \\\beta_t\\ to
\\\beta\_{t+1}\\ and the incremental log-weights are computed as

\\ \Delta \log w_i = (\beta\_{t+1} - \beta_t) \cdot \ell(\theta_i) \\

where \\\ell(\theta_i) = \log p(y \mid \theta_i)\\ is the
log-likelihood, extracted from the model function as \\-\text{Dev}/2\\.
These incremental weights are added to the running log-weights, which
are then normalized using the log-sum-exp trick to avoid numerical
overflow. The ESS is computed as \\\text{ESS} = 1 / \sum_i
\tilde{w}\_i^2\\ where \\\tilde{w}\_i\\ are the normalized weights.

### Resampling

When the ESS falls below the threshold `ESS.threshold * N.particles`,
the particle population is resampled to eliminate particles with
negligible weight and duplicate particles with large weight. **lucifer**
uses systematic resampling [\[3\]](#ref3), which partitions the unit
interval into \\N\\ equal segments and uses a single uniform random
number to select ancestors:

\\ u_i = \frac{U + i - 1}{N}, \quad U \sim \text{Uniform}(0, 1/N), \quad
i = 1, \ldots, N. \\

The ancestor of particle \\i\\ is \\a_i = \min\\j : \sum\_{k=1}^{j}
\tilde{w}\_k \geq u_i\\\\. Systematic resampling has lower variance than
multinomial resampling and is \\O(N)\\ in computational cost. After
resampling, the weights are reset to \\1/N\\.

### Rejuvenation

Resampling introduces duplicate particles, which reduces the diversity
of the population. The rejuvenation step applies an MCMC kernel that
targets the current tempered distribution \\\pi_t(\theta)\\ to each
particle, diversifying the population while preserving the target
distribution. In **lucifer**, the rejuvenation kernel is a random-walk
Metropolis (RWM) algorithm that proposes \\\theta^\* = \theta_i +
\varepsilon\\ where \\\varepsilon \sim \mathcal{N}(0, \sigma^2 \Sigma)\\
with \\\sigma^2 = 2.38^2 / d\\ (the optimal scaling for Gaussian targets
[\[4\]](#ref4)) and \\\Sigma\\ is the weighted covariance of the current
particle population. The acceptance probability for the tempered target
is

\\ \alpha(\theta_i, \theta^\*) = \min\left(1, \frac{p(\theta^\*) \\ p(y
\mid \theta^\*)^{\beta_t}}{p(\theta_i) \\ p(y \mid
\theta_i)^{\beta_t}}\right). \\

The number of rejuvenation steps per particle is controlled by
`Rejuvenation.steps` (default 5). Using multiple steps per stage
improves mixing but increases computational cost linearly; in practice,
3 to 10 steps provide a good tradeoff.

### Covariance adaptation

After resampling (and before rejuvenation), the weighted covariance of
the particle population is computed and used as the proposal covariance
for the rejuvenation kernel. This adaptation ensures that the proposal
distribution tracks the evolving shape of the tempered posterior as
\\\beta\\ increases. The covariance is only updated when the number of
unique particles exceeds \\d + 1\\, where \\d\\ is the parameter
dimension, to avoid degenerate estimates.

## Practical examples

### Setting up a simple model

The
[`SMC()`](https://robustecologies.github.io/lucifer/reference/SMC.md)
function shares the same model specification interface as
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
and
[`VariationalBayes()`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).
The model function must return a list with `LP` (log-posterior), `Dev`
(deviance), `Monitor`, `yhat`, and `parm`. The following example uses a
normal-normal conjugate model where the posterior is analytically
available for validation.

``` r

library(lucifer)

Model <- function(parm, Data) {
  mu <- parm[1]
  sigma <- interval(parm[2], 1e-100, Inf)
  parm[2] <- sigma
  LL <- sum(dnorm(Data$y, mu, sigma, log = TRUE))
  LP <- LL + dnorm(mu, 0, 10, log = TRUE) +
    dhalfcauchy(sigma, 25, log = TRUE)
  yhat <- rnorm(Data$n, mu, sigma)
  list(LP = LP, Dev = -2 * LL, Monitor = mu,
       yhat = yhat, parm = parm)
}

set.seed(42)
y <- rnorm(100, mean = 3, sd = 1.5)
MyData <- list(
  y = y, n = length(y),
  mon.names = "mu",
  parm.names = c("mu", "sigma")
)
```

### Running SMC with adaptive tempering

The simplest call uses the default adaptive schedule with 500 particles:

``` r

set.seed(123)
fit <- SMC(Model, MyData, Initial.Values = c(0, 1),
           N.particles = 500, ESS.threshold = 0.5,
           Rejuvenation.steps = 5, Status = 100)
```

The `print` method provides a concise overview of the fit:

``` r

print(fit)
```

The `summary` method provides extended diagnostics including the
tempering schedule, ESS trajectory, and rejuvenation acceptance rates:

``` r

summary(fit)
```

### Visualizing the results

The `plot` method supports three views. The default `type = "posterior"`
shows weighted marginal histograms of the particle positions:

``` r

plot(fit, type = "posterior")
```

The `type = "ess"` option shows how the effective sample size evolved
across tempering stages, with the resampling threshold indicated by a
dashed line:

``` r

plot(fit, type = "ess")
```

The `type = "schedule"` option displays the adaptive temperature
schedule, revealing where the algorithm concentrated its tempering
effort:

``` r

plot(fit, type = "schedule")
```

### Comparing adaptive vs. fixed schedules

The adaptive schedule concentrates tempering stages where the likelihood
is most informative, while a fixed schedule distributes them uniformly.
The following comparison illustrates the difference:

``` r

set.seed(123)
fit_adaptive <- SMC(Model, MyData, Initial.Values = c(0, 1),
                    N.particles = 500, Schedule = "adaptive",
                    ESS.threshold = 0.5, Status = 1000)

set.seed(123)
fit_fixed <- SMC(Model, MyData, Initial.Values = c(0, 1),
                 N.particles = 500, Schedule = "fixed",
                 ESS.threshold = 0.5, Status = 1000)

knitr::kable(data.frame(
  Schedule = c("Adaptive", "Fixed"),
  Stages = c(length(fit_adaptive$Schedule) - 1,
             length(fit_fixed$Schedule) - 1),
  `Log ML` = c(round(fit_adaptive$Log.Marginal.Likelihood, 2),
               round(fit_fixed$Log.Marginal.Likelihood, 2)),
  check.names = FALSE
), caption = "Adaptive vs. fixed tempering schedules.")
```

In most cases, the adaptive schedule achieves a comparable or better
marginal likelihood estimate with fewer tempering stages, because it
avoids wasting stages in regions where the tempered distribution changes
slowly.

### Using SMC for model comparison

Since
[`SMC()`](https://robustecologies.github.io/lucifer/reference/SMC.md)
produces an estimate of the log-marginal likelihood, it can be used
directly for Bayesian model comparison via Bayes factors. Suppose we
want to compare a model with unknown variance (as above) against a model
that fixes \\\sigma = 1\\:

``` r
Model_fixed_sigma <- function(parm, Data) {
  mu <- parm[1]
  LL <- sum(dnorm(Data$y, mu, 1, log = TRUE))
  LP <- LL + dnorm(mu, 0, 10, log = TRUE)
  yhat <- rnorm(Data$n, mu, 1)
  list(LP = LP, Dev = -2 * LL, Monitor = mu,
       yhat = yhat, parm = parm)
}

MyData2 <- list(
  y = y, n = length(y),
  mon.names = "mu",
  parm.names = "mu"
)

set.seed(123)
fit_fixed_sigma <- SMC(Model_fixed_sigma, MyData2,
                       Initial.Values = 0, N.particles = 500,
                       Status = 1000)

log_BF <- fit$Log.Marginal.Likelihood -
           fit_fixed_sigma$Log.Marginal.Likelihood

interpretation <- if (log_BF > 3) "Strong evidence for free sigma"
    else if (log_BF > 1) "Moderate evidence for free sigma"
    else if (log_BF > 0) "Weak evidence for free sigma"
    else "Evidence favors fixed sigma"

knitr::kable(data.frame(
  Metric = c("Log Bayes factor (free vs fixed sigma)", "Interpretation"),
  Value = c(round(log_BF, 2), interpretation)
), caption = "Bayes factor model comparison via SMC marginal likelihoods.")
```

## Tuning and diagnostics

### Number of particles

The number of particles \\N\\ controls the resolution of the posterior
approximation and the variance of the marginal likelihood estimate.
Increasing \\N\\ reduces the Monte Carlo error in the weight
computations and produces smoother posterior histograms, but increases
computational cost linearly. As a rule of thumb, \\N = 500\\ to \\N =
2000\\ is adequate for low-dimensional problems (up to roughly 10
parameters); higher-dimensional problems may require more particles or a
more efficient rejuvenation kernel.

### ESS threshold

The `ESS.threshold` parameter (default 0.5) controls the tradeoff
between the number of tempering stages and particle diversity. A higher
threshold (closer to 1) triggers resampling more frequently and produces
smaller temperature increments, which reduces weight degeneracy but
increases the total number of stages. A lower threshold (closer to 0)
allows more weight imbalance before resampling, producing fewer stages
but potentially losing effective particles. The value 0.5 is a widely
used default that balances these concerns [\[1\]](#ref1).

### Rejuvenation steps

More rejuvenation steps per stage improve the diversity of particles
after resampling, at the cost of additional model evaluations. With too
few steps, resampled duplicates remain close together and the particle
population degenerates into clusters. With too many steps, the
computation time per stage becomes the bottleneck. Monitoring the
rejection acceptance rate (reported in
[`summary()`](https://rdrr.io/r/base/summary.html)) provides guidance:
if the mean acceptance rate is below 0.1, the proposal covariance may be
too wide or the temperature increment too large; if it is above 0.5, the
rejuvenation steps could potentially be reduced.

### Diagnostics

The ESS trajectory (accessible via `fit$ESS.history` or
`plot(fit, type = "ess")`) is the primary diagnostic for SMC. A healthy
run shows ESS values that fluctuate near the threshold rather than
collapsing to near-zero. The tempering schedule (via `fit$Schedule` or
`plot(fit, type = "schedule")`) reveals whether the algorithm adapted
the temperature increments sensibly: a concave schedule (many small
steps early, fewer large steps later) is typical for models where the
likelihood is highly informative relative to the prior.

The acceptance rate history (`fit$Acceptance.rate`) tracks the
rejuvenation kernel’s mixing across stages. Declining acceptance rates
as \\\beta\\ approaches 1 are normal, since the tempered posterior
becomes more concentrated; consistently low rates (below 0.05) suggest
that the proposal covariance is poorly adapted or the problem
dimensionality is too high for random-walk proposals.

## IBIS: online sequential learning

The tempering approach described above treats the full dataset as a
fixed object and constructs a path from prior to posterior by raising
the likelihood to a fractional power \\\beta \in \[0, 1\]\\. Iterated
Batch Importance Sampling (IBIS), introduced by Chopin [\[2\]](#ref2),
takes a fundamentally different perspective: instead of tempering the
likelihood exponent, it incorporates observations sequentially in
batches, so that the sequence of intermediate distributions corresponds
to posteriors conditioned on progressively larger subsets of the data.
At stage \\b\\, the target is \\p(\theta \mid y\_{1:b}) \propto
p(\theta) \prod\_{k=1}^{b} p(y_k \mid \theta)\\, where \\y_k\\ denotes
the \\k\\-th batch of observations. The transition from stage \\b-1\\ to
stage \\b\\ is driven by a single new batch rather than by a change in
tempering exponent, which makes IBIS a natural framework for online and
streaming Bayesian inference where data arrives incrementally over time.

This sequential data-incorporation perspective has practical
consequences that distinguish IBIS from standard SMC tempering. Because
each batch shifts the posterior by an amount that depends on the
informativeness of the new observations, IBIS automatically adapts to
heterogeneous data: uninformative batches produce near-uniform weights
and require no resampling, while highly informative batches trigger
resampling and rejuvenation exactly when the posterior changes
substantially. The method also enables sequential forecasting, since
after processing batch \\b\\ the weighted particle population represents
a valid approximation to \\p(\theta \mid y\_{1:b})\\ from which
predictive distributions for future observations can be computed without
refitting the model from scratch.

### The IBIS algorithm

At each stage \\b = 1, \ldots, B\\, the particle population
\\\\\theta_i, w_i\\\_{i=1}^N\\ is updated to reflect the new batch
\\y_b\\. Each particle receives an incremental weight proportional to
the likelihood of the new batch under its current parameter value:

\\ \tilde{w}\_i^{(b)} = w_i^{(b-1)} \cdot \frac{p(\theta_i \mid
y\_{1:b})}{p(\theta_i \mid y\_{1:b-1})} = w_i^{(b-1)} \cdot p(y_b \mid
\theta_i, y\_{1:b-1}). \\

In practice, for models where observations are conditionally independent
given parameters, this simplifies to \\\tilde{w}\_i^{(b)} = w_i^{(b-1)}
\cdot \exp(\text{LP}\_{\text{new}} - \text{LP}\_{\text{old}})\\, where
\\\text{LP}\_{\text{new}}\\ and \\\text{LP}\_{\text{old}}\\ are the
log-posteriors evaluated with the cumulative data \\y\_{1:b}\\ and
\\y\_{1:b-1}\\ respectively. After reweighting, the ESS is computed from
the normalized weights; if it falls below the threshold, the particles
are resampled and rejuvenated using MCMC moves that target the current
cumulative posterior \\p(\theta \mid y\_{1:b})\\. The marginal
likelihood estimate accumulates across batches as \\\log \hat{p}(y) =
\sum\_{b=1}^{B} \log \left( \frac{1}{N} \sum\_{i=1}^{N}
\tilde{w}\_i^{(b)} / w_i^{(b-1)} \right)\\, which telescopes to give an
unbiased estimate of \\p(y\_{1:B})\\ on the natural scale, just as in
the tempering case.

### Usage in lucifer

IBIS mode is activated by setting `Schedule = "ibis"` in the call to
[`SMC()`](https://robustecologies.github.io/lucifer/reference/SMC.md).
The user must supply a list of data batches in `Data$batches`, where
each element contains the observations for that batch along with any
auxiliary information the model function requires (sample size, monitor
names, parameter names). The model function receives the cumulative data
at each stage, which
[`SMC()`](https://robustecologies.github.io/lucifer/reference/SMC.md)
constructs internally by concatenating batches. The per-batch weight
computation is handled by the C++ function
`smc_ibis_batch_weights_cpp()`, which evaluates the model at all
particle positions and computes the log-weight increments in a single
vectorized loop, avoiding R’s per-particle overhead. The following
example fits a normal model to five sequential batches of observations:

``` r

set.seed(42)
y_full <- rnorm(100, mean = 3, sd = 1.5)

Model <- function(parm, Data) {
  mu <- parm[1]
  sigma <- interval(parm[2], 1e-100, Inf)
  parm[2] <- sigma
  LL <- sum(dnorm(Data$y, mu, sigma, log = TRUE))
  LP <- LL + dnorm(mu, 0, 10, log = TRUE) +
    dhalfcauchy(sigma, 25, log = TRUE)
  yhat <- rnorm(Data$n, mu, sigma)
  list(LP = LP, Dev = -2 * LL, Monitor = mu,
       yhat = yhat, parm = parm)
}

batches <- lapply(split(y_full, cut(seq_along(y_full), 5)), function(yb) {
  list(y = yb, n = length(yb))
})

MyData <- list(
  y = y_full, n = length(y_full),
  mon.names = "mu",
  parm.names = c("mu", "sigma"),
  batches = batches
)

set.seed(123)
fit_ibis <- SMC(Model, MyData, Initial.Values = c(0, 1),
                N.particles = 1000, Schedule = "ibis",
                ESS.threshold = 0.5, Rejuvenation.steps = 10,
                Status = 100)

print(fit_ibis)
summary(fit_ibis)
plot(fit_ibis, type = "posterior")
```

### Practical considerations

IBIS typically requires more rejuvenation steps per stage than
tempering-based SMC, because a single batch of real observations can
shift the posterior more abruptly than a small increment in the
tempering exponent \\\beta\\. Setting `Rejuvenation.steps` to 10 or
higher is a reasonable starting point, and the acceptance rate history
from [`summary()`](https://rdrr.io/r/base/summary.html) should be
monitored to ensure adequate mixing after each batch incorporation.
Using more particles (1000 or above) also helps, since the weight
variance across batches can be substantial when individual batches are
highly informative. For sequential forecasting applications, the
weighted particle population after processing batch \\b\\ can be used
directly to compute predictive expectations: \\\hat{E}\[f(\theta) \mid
y\_{1:b}\] = \sum_i w_i^{(b)} f(\theta_i)\\, enabling rolling
predictions without the need to refit the model as new data arrives.

## Method selection guide

| Criterion | SMC | MCMC |
|:---|:---|:---|
| Marginal likelihood needed | Yes (unbiased) | Not directly |
| Embarrassingly parallel | Yes (per particle) | With multiple chains |
| No burn-in required | Yes | No (burn-in needed) |
| Multimodal posterior | Partial (if particles span modes) | Difficult (chain trapping) |
| High dimension (d \> 50) | Challenging (RWM kernel) | Yes (with HMC/NUTS) |
| Computational budget | Moderate to high | Low to moderate |

When to prefer SMC over MCMC in lucifer. {.table}

The primary use cases for SMC in **lucifer** are model comparison (where
the marginal likelihood estimate is essential), problems where MCMC
chains struggle with multimodality or poor mixing, and settings where
parallel evaluation of particles is available. For routine posterior
inference on unimodal, moderate-dimensional problems,
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
with an appropriate MCMC algorithm (NUTS, HMCDA, or adaptive Metropolis)
will typically be more efficient per model evaluation.

SMC output can also serve as initialization for MCMC: the weighted
particle mean and covariance from `fit$Summary` and
`cov.wt(fit$Posterior, fit$Weights)$cov` provide informed starting
values and a proposal covariance matrix for
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

## References

**\[1\]** Del Moral, P., Doucet, A., and Jasra, A. (2006). *Sequential
Monte Carlo samplers*. Journal of the Royal Statistical Society: Series
B, 68(3), 411-436. [DOI:
10.1111/j.1467-9868.2006.00553.x](https://doi.org/10.1111/j.1467-9868.2006.00553.x)

**\[2\]** Chopin, N. (2002). *A sequential particle filter method for
static models*. Biometrika, 89(3), 539-552. [DOI:
10.1093/biomet/89.3.539](https://doi.org/10.1093/biomet/89.3.539)

**\[3\]** Carpenter, J., Clifford, P., and Fearnhead, P. (1999).
*Improved particle filter for nonlinear problems*. IEE Proceedings:
Radar, Sonar and Navigation, 146(1), 2-7. [DOI:
10.1049/ip-rsn:19990255](https://doi.org/10.1049/ip-rsn:19990255)

**\[4\]** Roberts, G.O., Gelman, A., and Gilks, W.R. (1997). *Weak
convergence and optimal scaling of random walk Metropolis algorithms*.
Annals of Applied Probability, 7(1), 110-120. [DOI:
10.1214/aoap/1034625254](https://doi.org/10.1214/aoap/1034625254)

**\[5\]** Doucet, A., de Freitas, N., and Gordon, N. (2001). *Sequential
Monte Carlo methods in practice*. Springer. [DOI:
10.1007/978-1-4757-3437-9](https://doi.org/10.1007/978-1-4757-3437-9)

**\[6\]** Chopin, N., Papaspiliopoulos, O., and others (2020). *An
introduction to sequential Monte Carlo*. Springer. [DOI:
10.1007/978-3-030-47845-2](https://doi.org/10.1007/978-3-030-47845-2)
