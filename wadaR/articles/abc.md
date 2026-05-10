# Approximate Bayesian Computation in lucifer

## Introduction

Standard Bayesian inference requires evaluating the likelihood function
\\p(y \mid \theta)\\ at every proposed parameter value. For many complex
models, particularly in population genetics, ecology, epidemiology, and
agent-based modeling, the likelihood is either analytically intractable
or computationally prohibitive to evaluate, even though simulating data
from the model is straightforward. Approximate Bayesian Computation
(ABC) circumvents this barrier by replacing likelihood evaluation with
simulation: parameters are proposed, synthetic data are generated from
the model, and the proposal is accepted if the simulated data are “close
enough” to the observed data according to some distance metric on
summary statistics [\[1\]](#ref1).

The trade-off is explicit: ABC produces samples from an approximate
posterior \\p\_\varepsilon(\theta \mid y) \propto p(\theta) \\
\Pr(\rho(S(y\_\text{sim}), S(y\_\text{obs})) \leq \varepsilon)\\ rather
than the exact posterior, where \\\rho\\ is a distance function, \\S\\
computes summary statistics, and \\\varepsilon\\ is the tolerance. As
\\\varepsilon \to 0\\ and the summary statistics become sufficient, this
approximate posterior converges to the exact posterior; in practice, the
quality of the approximation depends critically on the choice of summary
statistics and the tolerance threshold.

The
[`ABC()`](https://robustecologies.github.io/lucifer/reference/ABC.md)
function in **lucifer** implements four ABC algorithms: rejection
sampling, MCMC-ABC [\[2\]](#ref2), SMC-ABC [\[3\]](#ref3),
[\[4\]](#ref4), and Simulated Annealing ABC (SA-ABC) [\[9\]](#ref9).
Unlike the other inference functions in the package,
[`ABC()`](https://robustecologies.github.io/lucifer/reference/ABC.md)
uses a different model specification where the model function returns
simulated data rather than a log-posterior.

  

## The ABC framework

### The rejection principle

The simplest ABC algorithm proceeds by rejection sampling. For each
iteration: draw \\\theta^\* \sim p(\theta)\\ from the prior, simulate a
dataset \\y^\* \sim p(y \mid \theta^\*)\\ from the model, compute
summary statistics \\S(y^\*)\\, and accept \\\theta^\*\\ if
\\\rho(S(y^\*), S(y\_\text{obs})) \leq \varepsilon\\. The accepted
samples are draws from the approximate posterior \\p\_\varepsilon(\theta
\mid S(y\_\text{obs}))\\.

The acceptance probability under rejection ABC equals
\\\Pr(\rho(S(y\_\text{sim}), S(y\_\text{obs})) \leq \varepsilon \mid
\theta)\\ integrated against the prior. When \\\varepsilon\\ is small,
the acceptance rate becomes extremely low, particularly in moderate to
high dimensions where the prior places negligible mass near the
posterior. This motivates the MCMC and SMC extensions.

### Summary statistics

The choice of summary statistics \\S(y)\\ is arguably the most
consequential design decision in ABC. Sufficient statistics, when they
exist, yield an approximate posterior that converges to the exact
posterior as \\\varepsilon \to 0\\. In practice, sufficient statistics
are rarely available for the complex models that motivate ABC in the
first place. The typical strategy is to use a vector of scientifically
motivated quantities: sample moments, quantiles, autocorrelations,
spectral features, or other domain-specific summaries.

Using too few statistics loses information and widens the approximate
posterior; using too many increases the effective dimensionality of the
comparison and makes the distance criterion harder to satisfy, driving
down acceptance rates. The `Summary.Stats` argument in
[`ABC()`](https://robustecologies.github.io/lucifer/reference/ABC.md)
accepts any function that maps a dataset to a numeric vector, giving the
user full control over this choice.

### Distance and tolerance

The distance function \\\rho\\ in **lucifer** is the Euclidean norm of
the difference between simulated and observed summary statistics:
\\\rho(S(y^\*), S(y\_\text{obs})) = \\S(y^\*) - S(y\_\text{obs})\\\_2\\.
When summary statistics have different scales, users should standardize
them within the `Summary.Stats` function (e.g., dividing each statistic
by its expected standard deviation under the prior predictive) to
prevent high-variance statistics from dominating the distance.

The tolerance \\\varepsilon\\ controls the accuracy-efficiency tradeoff.
If `epsilon = NULL`, **lucifer** performs a pilot run (1000
prior-predictive simulations for rejection ABC, 200 for MCMC-ABC) and
sets \\\varepsilon\\ to the 1st percentile of the pilot distance
distribution. This adaptive calibration provides a reasonable starting
point, but users should verify that the resulting acceptance rate is not
too low (below 0.001) or too high (above 0.1), adjusting \\\varepsilon\\
accordingly.

  

## Rejection ABC

### Algorithm

Rejection ABC is the simplest and most parallelizable ABC method. The
implementation in **lucifer** proceeds as follows:

1.  If `epsilon` is `NULL`, run a pilot of \\\min(1000, N)\\
    prior-predictive simulations and set \\\varepsilon\\ to the 1st
    percentile of the distance distribution.
2.  Repeat until \\N\\ samples are accepted or the total number of
    proposals exceeds \\1000 \times N\\ (a safety limit):
    - Draw \\\theta^\* \sim p(\theta)\\ via the user-supplied `Prior()`
      function.
    - Simulate \\y^\* = \text{Model}(\theta^\*, \text{Data})\\.
    - Compute \\S(y^\*)\\ via `Summary.Stats()` and the distance \\d =
      \\S(y^\*) - S(y\_\text{obs})\\\_2\\.
    - Accept \\\theta^\*\\ if \\d \leq \varepsilon\\.
3.  Return the matrix of accepted parameters, their distances, and the
    acceptance rate.

The computational cost is \\N / r\\ model simulations where \\r\\ is the
acceptance rate. For low-dimensional problems with informative priors,
rejection ABC can be efficient; for high-dimensional or diffuse-prior
problems, the acceptance rate drops precipitously and more sophisticated
methods are needed.

  

## MCMC-ABC

### The Marjoram algorithm

MCMC-ABC, introduced by Marjoram et al. [\[2\]](#ref2), replaces the
independent prior sampling of rejection ABC with a Markov chain that
explores the parameter space more efficiently. At each iteration, a
proposal \\\theta^\*\\ is drawn from a random-walk kernel centered at
the current state \\\theta\\, data are simulated from the model, and the
proposal is accepted if the simulated summary statistics fall within
\\\varepsilon\\ of the observed ones. The key insight is that the
standard Metropolis-Hastings acceptance criterion simplifies
dramatically when using a uniform kernel (indicator function) for the
ABC likelihood: the acceptance probability reduces to

\\ \alpha(\theta, \theta^\*) = \begin{cases} 1 & \text{if }
\rho(S(y^\*), S(y\_\text{obs})) \leq \varepsilon \\ 0 & \text{otherwise}
\end{cases} \\

combined with the prior ratio (which equals 1 for a symmetric proposal
and uniform prior support). The chain therefore moves to \\\theta^\*\\
whenever the simulated data are sufficiently close to the observations,
and stays at \\\theta\\ otherwise.

### Implementation details

The **lucifer** implementation uses a random-walk Metropolis proposal
with covariance \\\sigma^2 \Sigma\\, where \\\sigma^2 = 2.38^2 / d\\ is
the optimal scaling factor [\[5\]](#ref5) and \\\Sigma\\ is either
user-supplied via `Covar` or defaults to the identity matrix. The
Cholesky factor of \\\sigma^2 \Sigma\\ is precomputed for efficient
sampling.

MCMC-ABC shares a fundamental limitation with all MCMC methods: it
produces correlated samples and requires care with burn-in. The first
portion of the chain should be discarded, and thinning may be necessary
if the acceptance rate is very low (since the chain then spends long
stretches at the same parameter value). The advantage over rejection ABC
is that proposals are made near the current accepted state rather than
from the diffuse prior, which dramatically improves the acceptance rate
for problems where the posterior is concentrated.

  

## SMC-ABC

### Population-based refinement

SMC-ABC [\[3\]](#ref3), [\[4\]](#ref4) combines the population-based
strategy of SMC with the likelihood-free approach of ABC. Instead of
targeting a single tolerance, it constructs a sequence of decreasing
tolerances \\\varepsilon_1 \> \varepsilon_2 \> \cdots \> \varepsilon_T\\
and evolves a particle population through these intermediate ABC
posteriors, analogous to the tempering approach used in SMC samplers
(see
[`vignette("smc")`](https://robustecologies.github.io/lucifer/articles/smc.md)).

### Algorithm

The **lucifer** implementation proceeds through the following stages:

1.  **Initialization.** Draw \\N\\ particles \\\\\theta_i\\\\ from the
    prior and compute distances \\d_i = \rho(S(\text{Model}(\theta_i,
    \text{Data})), S(y\_\text{obs}))\\ for each.

2.  **Tolerance schedule.** If `epsilon` is `NULL`, construct a
    decreasing schedule by taking quantiles of the initial distance
    distribution: \\\varepsilon_t \in \\\text{quantile}(d, 0.5), \ldots,
    \text{quantile}(d, 0.01)\\\\ over 10 stages. If `epsilon` is
    provided, the schedule interpolates linearly from the median initial
    distance to the target.

3.  **For each stage** \\t = 1, \ldots, T\\:

    - **Reweight.** Assign weight 1 to particles with \\d_i \leq
      \varepsilon_t\\ and weight 0 to the rest; normalize.
    - **Resample.** Draw \\N\\ particles with replacement according to
      the weights (multinomial resampling).
    - **Rejuvenate.** For each particle, propose a move using a
      random-walk kernel with covariance estimated from the current
      particle population (twice the empirical covariance), simulate
      data from the proposed parameters, and accept if the distance is
      below \\\varepsilon_t\\. This is identical to a single step of
      MCMC-ABC at tolerance \\\varepsilon_t\\.
    - **Check degeneracy.** If fewer than 2 particles satisfy the
      tolerance, terminate early to avoid degenerate populations.

4.  **Return** the final particle positions, distances, acceptance rates
    across stages, and the posterior summary.

### Comparison with MCMC-ABC

SMC-ABC has several advantages over single-chain MCMC-ABC. It avoids the
burn-in problem, since the population is initialized from the prior and
progressively refined. It naturally adapts the tolerance schedule to the
problem. And it provides a population of samples rather than a single
correlated chain, which gives a more reliable picture of the posterior.
The cost is higher per stage (since all \\N\\ particles must be
simulated), but the total number of model simulations is often
comparable to MCMC-ABC because the gradually decreasing tolerance
prevents the low acceptance rates that plague MCMC-ABC with a tight
fixed tolerance.

  

## SA-ABC (Simulated Annealing)

### Background and motivation

Simulated annealing (SA) was introduced by Kirkpatrick, Gelatt, and
Vecchi [\[10\]](#ref10) as a stochastic optimization method inspired by
the physical process of cooling a molten material into a crystalline
state. The algorithm explores a landscape by accepting both downhill
moves (which improve the objective) and uphill moves (which worsen it),
with the probability of accepting uphill moves controlled by a
temperature parameter \\T\\ that decreases over time. As \\T \to 0\\,
the algorithm converges to a local (ideally global) minimum. In its
classical form, SA is an optimizer, not a sampler: it produces a point
estimate, not a distribution.

Albert, Kuensch, and Scheidegger [\[9\]](#ref9) adapted the SA idea to
ABC by replacing the hard indicator acceptance of MCMC-ABC with a soft,
temperature-modulated Boltzmann kernel \\\exp(-d(\theta) /
\varepsilon)\\. Their SABC algorithm, implemented in the **EasyABC**
package [\[14\]](#ref14) (co-authored by Albert himself), operates on a
population of particles and uses an adaptive cooling schedule derived
from the principle of minimal entropy production: the temperature
decreases as a function of the ensemble’s mean energy, governed by an
implicit equation \\\varepsilon^2 + v\varepsilon^{3/2} = \bar{\rho}^2\\
where \\v\\ controls the annealing speed and \\\bar{\rho}\\ is the mean
transformed distance. The algorithm includes periodic importance
resampling of the ensemble (analogous to SMC-ABC), adaptive proposal
covariance estimation from the full particle population, and a CDF
reparametrization that maps raw distances through the empirical
distribution function of prior-predictive simulations to standardize
them to \\\[0, 1\]\\. An “informative” variant extends the framework to
a two-temperature system \\(T_1, T_2)\\ that separately controls the
data-fit tolerance and a prior-density correction, tracked via tensor
algebra. The full Albert et al. algorithm requires an evaluable prior
density \\\pi(\theta)\\ for the Metropolis ratio, not merely a prior
sampler.

A related but distinct approach is simulated tempering [\[11\]](#ref11),
[\[12\]](#ref12), where the temperature itself becomes a random variable
sampled jointly with the parameters, enabling transitions between hot
(exploratory) and cold (concentrated) regimes within a single chain.

### Relationship to the Albert et al. algorithm

The SA-ABC implementation in **lucifer** departs from the full Albert et
al. algorithm in several deliberate ways. It uses a single chain rather
than a particle ensemble, a fixed cooling schedule rather than adaptive
entropy-minimizing temperature control, and raw Euclidean distances
rather than CDF-transformed energies. Most importantly, it does not
require an evaluable prior density function: the
[`ABC()`](https://robustecologies.github.io/lucifer/reference/ABC.md)
interface accepts only a `Prior()` sampler, consistent with the other
three ABC methods in the package. These simplifications are not
incidental; they reflect a different design goal.

The full SABC algorithm maintains a soft Boltzmann kernel
\\\exp(-d/\varepsilon)\\ throughout the entire run. This kernel defines
a proper target distribution \\\pi\_\varepsilon(\theta) \propto
p(\theta) \exp(-d(\theta)/\varepsilon)\\ at every temperature level, and
the population structure implicitly estimates the normalizing constant
ratios that a single-chain approach cannot access. The resulting
posterior, however, is a soft-kernel ABC posterior that differs from the
hard-threshold posterior \\p\_\varepsilon(\theta \mid S(y\_\text{obs}))
\propto p(\theta) \\ I(d \leq \varepsilon)\\ targeted by rejection ABC,
MCMC-ABC, and SMC-ABC. In practice, the exponential kernel assigns
non-negligible weight to proposals well beyond the tolerance: at \\d =
2\varepsilon\\, the weight is still \\\exp(-2) \approx 0.14\\, compared
to zero under the hard indicator. This makes posterior summaries from a
pure soft-kernel SABC not directly comparable to those from
hard-threshold methods.

The **lucifer** implementation resolves this with a two-phase hybrid
design that borrows the SA acceptance formula from Albert et al. for the
exploratory phase, then switches to the same hard indicator kernel used
by the other ABC methods for the sampling phase. This follows the
“anneal then equilibrate” pattern established in computational physics
and molecular dynamics, where a system is first annealed to locate the
energy basin and then simulated at constant conditions to sample the
equilibrium distribution. The tradeoff is explicit: **lucifer**
sacrifices the theoretical coherence of a single soft-kernel target (and
the adaptive cooling, ensemble structure, and CDF reparametrization that
support it) in exchange for posterior comparability across all four ABC
methods and a simpler interface that requires no prior density
evaluator.

### Acceptance criterion

During the annealing phase, the acceptance probability for a proposed
move from \\\theta\\ to \\\theta^\*\\ is

\\ \alpha(\theta, \theta^\*) = \min\\\left(1, \\
\exp\\\left(-\frac{d(\theta^\*) - d(\theta)}{T_t}\right)\right) \\

where \\d(\cdot) = \\S(y^\*) - S(y\_\text{obs})\\\_2\\ is the Euclidean
distance between simulated and observed summary statistics and \\T_t\\
is the temperature at iteration \\t\\. Moves that reduce the distance
(downhill) are always accepted; moves that increase it (uphill) are
accepted with a probability that decreases as temperature drops. When
\\T_t\\ is large, the chain moves freely through the distance landscape;
as \\T_t\\ approaches \\\varepsilon\\, the criterion becomes
increasingly selective.

### Two-phase design

The algorithm operates in two phases:

- **Phase 1** (annealing, first 30% of iterations): the temperature
  cools from \\T_0\\ to \\T\_\text{final} = \varepsilon\\, driving the
  chain toward the low-distance region. The soft acceptance criterion
  allows the chain to escape local optima via uphill moves at high
  temperature, traversing barriers in the distance landscape that would
  permanently trap a hard-threshold MCMC-ABC chain. The proposal kernel
  is adapted continuously during annealing (every 200 iterations) using
  the empirical covariance of recent samples.

- **Phase 2** (sampling, remaining 70%): the algorithm switches to the
  standard hard ABC threshold \\d(\theta^\*) \leq \varepsilon\\,
  identical to MCMC-ABC. The chain now targets the standard ABC
  posterior \\p\_\varepsilon(\theta \mid S(y\_\text{obs})) \propto
  p(\theta) \\ I(\rho(S(y^\*), S(y\_\text{obs})) \leq \varepsilon)\\ and
  produces samples that are directly comparable to those from the other
  three methods.

At the transition between phases, the chain is reinitialized at a point
satisfying \\d \leq \varepsilon\\. If the annealing phase produced any
samples within the tolerance, the best one is used; otherwise, a short
burst of rejection sampling from the prior locates a valid starting
point. The proposal kernel is adapted using the empirical covariance of
annealing samples that fell within \\\varepsilon\\, calibrating the
proposal scale for the target region rather than the broad landscape
explored during annealing.

### Cooling schedules

**lucifer** provides two schedules for Phase 1:

- **Geometric cooling:** \\T_t = T_0 \cdot r^t\\ where \\r \in (0,1)\\
  is the cooling rate. This produces exponential decay with most of the
  cooling concentrated in early iterations.
- **Linear cooling:** \\T_t\\ interpolates linearly from \\T_0\\ to
  \\\varepsilon\\ over the annealing phase.

By default, both the initial temperature \\T_0\\ and the cooling rate
\\r\\ are auto-calibrated from a pilot run of 1000 prior-predictive
simulations. \\T_0\\ is set to the median pilot distance, which ensures
that at the start approximately 37% of uphill moves with distance equal
to the median are accepted (since \\\exp(-1) \approx 0.37\\). For
geometric cooling, the rate \\r\\ is computed so the temperature reaches
\\\varepsilon\\ at the end of Phase 1: \\r = (\varepsilon /
T_0)^{1/n\_\text{anneal}}\\, clamped to \\\[0.9, 0.9999\]\\ for
numerical stability.

### Posterior extraction and comparability

The `Posterior` field returned by SA-ABC contains only Phase 2 samples
(unlike MCMC-ABC, where the user must discard burn-in manually). Because
Phase 2 uses the same hard indicator kernel as MCMC-ABC, the SA-ABC
posterior targets the same approximate distribution
\\p\_\varepsilon(\theta \mid S(y\_\text{obs}))\\ as rejection ABC,
MCMC-ABC, and SMC-ABC. Phase 1 serves as a structured,
temperature-guided burn-in that replaces the random initialization of
MCMC-ABC with a directed search through the distance landscape, while
the reinitialization at the phase boundary guarantees that Phase 2
begins from a valid state.

This design trades the theoretical elegance of a single coherent target
distribution (as in the full Albert et al. population-based SABC) for
practical comparability across all four ABC methods in the package. A
user can meaningfully compare posterior summaries from
`Method = "rejection"`, `"MCMC"`, `"SMC"`, and `"SA"` knowing that all
four target the same ABC posterior, differing only in how they explore
the parameter space to reach it.

### Alternative approaches not implemented

Two temperature-based strategies were considered and deliberately
excluded from **lucifer**.

**Simulated tempering** [\[11\]](#ref11), [\[12\]](#ref12) treats the
temperature as an auxiliary random variable sampled jointly with the
parameters, so the chain can revisit high temperatures at any point
during sampling rather than only during a fixed burn-in phase. This
would allow the chain to escape local optima discovered after the
annealing phase has ended, a genuine limitation of the two-phase SA-ABC
design. The obstacle is the normalizing constant problem: simulated
tempering requires the ratios \\Z(T_i) / Z(T_j)\\ between partition
functions at each temperature level, which in ABC involve a
simulation-based kernel with no closed-form normalizer. Adaptive
estimation strategies such as Wang-Landau [\[13\]](#ref13) or the
Geyer-Thompson exchange algorithm introduce their own convergence
complications. Albert et al.’s population-based SABC [\[9\]](#ref9)
sidesteps this problem because the ensemble’s importance weights
implicitly estimate the normalizing constant ratios; a single-chain
simulated tempering implementation does not have this structural
advantage.

**Parallel tempering** (replica exchange) would run \\K\\ chains at
different temperatures and swap adjacent pairs periodically, requiring
no normalizing constant estimation. However, it demands \\K\\
simultaneous model evaluations per iteration, which is expensive for ABC
where each evaluation requires a full forward simulation, and its
benefits overlap substantially with those of SMC-ABC.

**The full Albert et al. SABC** (as implemented in **EasyABC**
[\[14\]](#ref14)) would provide a theoretically complete
population-based approach with adaptive cooling, CDF reparametrization,
and importance resampling. Adopting it would, however, essentially
duplicate the SMC-ABC functionality already in **lucifer** (population +
resampling + adaptive tolerance) with a different kernel, require an
evaluable prior density (breaking the ABC interface contract that only
demands a sampler), and produce posterior samples from a soft kernel not
directly comparable to the other three methods.

The four methods currently in **lucifer** cover the ABC design space
without redundancy: rejection sampling (simple, embarrassingly
parallel), MCMC-ABC (efficient for concentrated posteriors), SMC-ABC
(robust, adaptive, population-based), and SA-ABC (directed
initialization for landscapes with local traps).

  

## Practical examples

### Inferring the mean and variance of a normal distribution

This example demonstrates all three ABC methods on a simple problem
where the posterior is analytically available. The model simulates
normal data, and the summary statistics are the sample mean and standard
deviation.

``` r

library(lucifer)

# Simulator: returns a vector of simulated data
# sigma <= 0 produces extreme data, so proposals with invalid sigma
# are automatically rejected by the distance criterion
Simulator <- function(parm, Data) {
  mu <- parm[1]
  sigma <- parm[2]
  if (sigma <= 0) return(rep(1e10, Data$n))
  rnorm(Data$n, mu, sigma)
}

# Prior: mu ~ N(0, 10), sigma ~ Half-Cauchy(0, 5)
Prior <- function() {
  c(rnorm(1, 0, 10), abs(rcauchy(1, 0, 5)))
}

# Summary statistics: mean and SD
Summary.Stats <- function(data) {
  c(mean(data), sd(data))
}

# Observed data
set.seed(42)
obs_data <- rnorm(100, mean = 3, sd = 1.5)
Observed.Stats <- Summary.Stats(obs_data)

Data <- list(n = 100, parm.names = c("mu", "sigma"))
```

### Rejection ABC

``` r

set.seed(123)
fit_rej <- ABC(Simulator, Data, Prior, Summary.Stats, Observed.Stats,
               Method = "rejection", N = 2000, Status = 10000)
print(fit_rej)
true_vals <- c(mu = 3, sigma = 1.5)
plot(fit_rej, type = "posterior", true_values = true_vals)
plot(fit_rej, type = "distances")
```

The posterior summary shows the estimated mean and credible intervals.
With 100 observations and the sample mean and SD as summary statistics
(which are jointly sufficient for this model), the ABC posterior should
closely match the exact posterior.

### MCMC-ABC

``` r

set.seed(123)
fit_mcmc <- ABC(Simulator, Data, Prior, Summary.Stats, Observed.Stats,
                Method = "MCMC", MCMC.iterations = 10000,
                Covar = diag(c(0.3, 0.2)), Status = 5000)
print(fit_mcmc)
plot(fit_mcmc, type = "posterior", true_values = true_vals)
plot(fit_mcmc, type = "distances")
```

MCMC-ABC typically achieves a higher effective sample size per model
simulation than rejection ABC, since proposals are local rather than
from the prior. However, the chain produces correlated samples and
requires burn-in removal. The first portion of the chain should be
discarded before computing posterior summaries; for the comparison below
we discard the first 20% of iterations.

### SMC-ABC

``` r

set.seed(123)
fit_smc <- ABC(Simulator, Data, Prior, Summary.Stats, Observed.Stats,
               Method = "SMC", N = 1000, Status = 100)
print(fit_smc)
plot(fit_smc, type = "posterior", true_values = true_vals)
plot(fit_smc, type = "distances")
```

### SA-ABC

``` r

set.seed(123)
fit_sa <- ABC(Simulator, Data, Prior, Summary.Stats, Observed.Stats,
              Method = "SA", SA.iterations = 5000, Status = 1000)
print(fit_sa)
```

SA-ABC avoids the hard threshold problem of MCMC-ABC by allowing uphill
moves early in the chain when temperature is high, then progressively
concentrating on low-distance regions. The `Posterior` field contains
only Phase 2 (sampling phase) samples, so no manual burn-in removal is
needed. The `type = "temperature"` plot shows the cooling schedule
alongside the distance trace:

``` r

plot(fit_sa, type = "temperature")
plot(fit_sa, type = "posterior", true_values = true_vals)
plot(fit_sa, type = "distances")
```

### Comparing methods

``` r

# Apply burn-in for MCMC (discard first 20%)
burnin <- floor(0.2 * nrow(fit_mcmc$Posterior))
mcmc_post <- fit_mcmc$Posterior[(burnin + 1):nrow(fit_mcmc$Posterior), ,
                                drop = FALSE]

# SA: Posterior already contains Phase 2 (sampling phase) only
sa_post <- fit_sa$Posterior

comparison <- data.frame(
  Method = c("True", "Rejection", "MCMC", "SMC", "SA"),
  mu.mean = c(3,
              mean(fit_rej$Posterior[, "mu"]),
              mean(mcmc_post[, "mu"]),
              mean(fit_smc$Posterior[, "mu"]),
              mean(sa_post[, "mu"])),
  mu.sd = c(NA,
            sd(fit_rej$Posterior[, "mu"]),
            sd(mcmc_post[, "mu"]),
            sd(fit_smc$Posterior[, "mu"]),
            sd(sa_post[, "mu"])),
  sigma.mean = c(1.5,
                 mean(fit_rej$Posterior[, "sigma"]),
                 mean(mcmc_post[, "sigma"]),
                 mean(fit_smc$Posterior[, "sigma"]),
                 mean(sa_post[, "sigma"])),
  sigma.sd = c(NA,
               sd(fit_rej$Posterior[, "sigma"]),
               sd(mcmc_post[, "sigma"]),
               sd(fit_smc$Posterior[, "sigma"]),
               sd(sa_post[, "sigma"])),
  Accept.rate = c(NA,
                  fit_rej$Acceptance.rate,
                  fit_mcmc$Acceptance.rate,
                  fit_smc$Acceptance.rate,
                  fit_sa$Acceptance.rate),
  row.names = NULL
)
knitr::kable(comparison, digits = 3,
  caption = "Posterior estimates across ABC methods.")
```

### Visualizing posteriors

The `plot` method displays marginal posterior histograms by default:

``` r

plot(fit_rej, type = "posterior", true_values = true_vals)
```

The `type = "distances"` option shows the distribution of distances for
accepted samples, with the tolerance threshold marked:

``` r

plot(fit_rej, type = "distances")
```

### A non-trivial example: MA(2) process

ABC is most useful when the likelihood is intractable. A classic example
from the ABC literature is inference for a moving average process of
order 2, where the likelihood involves a sum over all possible latent
states. The summary statistics are the sample autocovariances at lags 1
and 2, which are informative but not sufficient.

``` r

# MA(2) simulator
MA2_sim <- function(parm, Data) {
  theta1 <- parm[1]
  theta2 <- parm[2]
  n <- Data$n
  e <- rnorm(n + 2)
  y <- numeric(n)
  for (t in 1:n) {
    y[t] <- e[t + 2] + theta1 * e[t + 1] + theta2 * e[t]
  }
  y
}

# Prior: uniform on the invertibility triangle
MA2_prior <- function() {
  repeat {
    theta <- runif(2, -2, 2)
    # Invertibility constraints: theta2 + theta1 > -1,
    # theta2 - theta1 > -1, |theta2| < 1
    if (theta[2] + theta[1] > -1 &&
        theta[2] - theta[1] > -1 &&
        abs(theta[2]) < 1) return(theta)
  }
}

# Summary: autocovariances at lags 1, 2
MA2_stats <- function(y) {
  n <- length(y)
  ybar <- mean(y)
  g1 <- sum((y[-1] - ybar) * (y[-n] - ybar)) / n
  g2 <- sum((y[-(1:2)] - ybar) * (y[-((n-1):n)] - ybar)) / n
  c(g1, g2)
}

# True parameters: theta1 = 0.6, theta2 = 0.2
set.seed(42)
obs_y <- MA2_sim(c(0.6, 0.2),
                 list(n = 500))
obs_stats <- MA2_stats(obs_y)
MA2_Data <- list(n = 500,
                 parm.names = c("theta1", "theta2"))
ma2_true <- c(theta1 = 0.6, theta2 = 0.2)
```

We now fit all four ABC methods to this problem and compare their
posterior estimates against the true parameter values.

``` r

set.seed(123)
fit_ma2_rej <- ABC(MA2_sim, MA2_Data, MA2_prior, MA2_stats, obs_stats,
                   Method = "rejection", N = 2000, Status = 10000)
```

``` r

set.seed(123)
fit_ma2_mcmc <- ABC(MA2_sim, MA2_Data, MA2_prior, MA2_stats, obs_stats,
                    Method = "MCMC", MCMC.iterations = 10000,
                    Covar = diag(c(0.1, 0.1)), Status = 5000)
```

``` r

set.seed(123)
fit_ma2_smc <- ABC(MA2_sim, MA2_Data, MA2_prior, MA2_stats, obs_stats,
                   Method = "SMC", N = 2000, Status = 100)
```

``` r

set.seed(123)
fit_ma2_sa <- ABC(MA2_sim, MA2_Data, MA2_prior, MA2_stats, obs_stats,
                  Method = "SA", SA.iterations = 5000, Status = 1000)
```

``` r

plot(fit_ma2_rej, type = "posterior", true_values = ma2_true)
plot(fit_ma2_mcmc, type = "posterior", true_values = ma2_true)
plot(fit_ma2_smc, type = "posterior", true_values = ma2_true)
plot(fit_ma2_sa, type = "posterior", true_values = ma2_true)
```

``` r

# Apply burn-in for MCMC (discard first 20%)
ma2_mcmc_burn <- floor(0.2 * nrow(fit_ma2_mcmc$Posterior))
ma2_mcmc_post <- fit_ma2_mcmc$Posterior[(ma2_mcmc_burn + 1):nrow(fit_ma2_mcmc$Posterior), ,
                                         drop = FALSE]

# SA: Posterior already contains Phase 2 (sampling phase) only
ma2_sa_post <- fit_ma2_sa$Posterior

ma2_comparison <- data.frame(
  Method = c("True", "Rejection", "MCMC", "SMC", "SA"),
  theta1.mean = c(0.6,
                  mean(fit_ma2_rej$Posterior[, "theta1"]),
                  mean(ma2_mcmc_post[, "theta1"]),
                  mean(fit_ma2_smc$Posterior[, "theta1"]),
                  mean(ma2_sa_post[, "theta1"])),
  theta1.sd = c(NA,
                sd(fit_ma2_rej$Posterior[, "theta1"]),
                sd(ma2_mcmc_post[, "theta1"]),
                sd(fit_ma2_smc$Posterior[, "theta1"]),
                sd(ma2_sa_post[, "theta1"])),
  theta2.mean = c(0.2,
                  mean(fit_ma2_rej$Posterior[, "theta2"]),
                  mean(ma2_mcmc_post[, "theta2"]),
                  mean(fit_ma2_smc$Posterior[, "theta2"]),
                  mean(ma2_sa_post[, "theta2"])),
  theta2.sd = c(NA,
                sd(fit_ma2_rej$Posterior[, "theta2"]),
                sd(ma2_mcmc_post[, "theta2"]),
                sd(fit_ma2_smc$Posterior[, "theta2"]),
                sd(ma2_sa_post[, "theta2"])),
  Accept.rate = c(NA,
                  fit_ma2_rej$Acceptance.rate,
                  fit_ma2_mcmc$Acceptance.rate,
                  fit_ma2_smc$Acceptance.rate,
                  fit_ma2_sa$Acceptance.rate),
  row.names = NULL
)
knitr::kable(ma2_comparison, digits = 3,
  caption = "MA(2) posterior estimates across ABC methods.")
```

This example illustrates a case where MCMC with a standard likelihood is
impractical, but ABC provides a principled and computationally feasible
route to posterior inference. The autocovariance summary statistics are
not sufficient for the MA(2) model, so the ABC posteriors are
necessarily wider than the exact posterior; nevertheless, all four
methods recover the true parameter values within their credible
intervals.

  

## Tuning guidelines

### Choosing summary statistics

The selection of summary statistics determines the quality of the ABC
approximation. Insufficient statistics lead to a posterior that is wider
than the exact posterior, a phenomenon known as the “loss of
information” or “summary statistic curse.” In the normal-mean example
above, using only the sample mean (ignoring the sample SD) would produce
correct inference for \\\mu\\ but would provide no information about
\\\sigma\\.

When natural sufficient statistics are unavailable, consider using a
moderate number of scientifically motivated summaries (5 to 15) rather
than either too few or too many. Semi-automatic methods such as the
regression adjustment of Beaumont et al. [\[4\]](#ref4) or neural
network-based summary selection [\[6\]](#ref6) can help, but these are
beyond the scope of the current implementation.

### Setting the tolerance

The tolerance controls accuracy at the expense of acceptance rate. Lower
\\\varepsilon\\ produces a better approximation to the exact posterior
but requires more simulations to achieve the same number of accepted
samples. The adaptive pilot-run calibration in **lucifer** (1st
percentile of the prior-predictive distance distribution) provides a
reasonable default, but users should inspect the resulting acceptance
rate and posterior width.

For rejection ABC, acceptance rates between 0.001 and 0.01 are typical;
rates below 0.0001 suggest that \\\varepsilon\\ is too tight or the
summary statistics are too high-dimensional. For MCMC-ABC, acceptance
rates between 0.05 and 0.3 are healthy. For SMC-ABC, the per-stage
acceptance rate during rejuvenation should be above 0.1.

### When to use each method

| Criterion | Rejection | MCMC | SMC | SA |
|:---|:---|:---|:---|:---|
| Dimensionality | Low (d \<= 3) | Moderate (d \<= 10) | Moderate (d \<= 10) | Moderate (d \<= 10) |
| Acceptance efficiency | Low (prior sampling) | Moderate (local proposals) | High (population-based) | High (soft acceptance) |
| Correlated samples | No (independent) | Yes | Low (resampled population) | Yes |
| Burn-in required | No | Yes | No | No (Phase 2 only) |
| Parallelizable | Yes (trivially) | No (single chain) | Yes (per particle) | No (single chain) |
| Adaptive tolerance | No | No | Yes (decreasing schedule) | Yes (cooling schedule) |

Comparison of ABC methods in lucifer. {.table style="width:100%;"}

Rejection ABC is the simplest option and should be tried first for
low-dimensional problems with an informative prior. When the prior is
diffuse or the parameter space has more than a few dimensions, MCMC-ABC
provides a more efficient exploration. SMC-ABC is the most robust choice
for moderate-dimensional problems and has the advantage of adapting the
tolerance automatically, but it requires the most model simulations per
run. SA-ABC is particularly useful when the distance landscape has
multiple local optima or when MCMC-ABC gets trapped: the annealing phase
uses temperature-modulated acceptance to traverse barriers in the
distance landscape, then the sampling phase switches to standard
hard-threshold MCMC-ABC for proper posterior inference. Because all four
methods target the same ABC posterior \\p\_\varepsilon(\theta \mid
S(y\_\text{obs}))\\, their results are directly comparable.

Regardless of method, ABC output should be interpreted with caution. The
approximate posterior is only as good as the summary statistics and the
tolerance allow. For final inference in critical applications, users
should verify that the results are robust to different choices of
summary statistics and tolerance values. When feasible, comparing ABC
results against a likelihood-based method (e.g.,
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
with MCMC or
[`SMC()`](https://robustecologies.github.io/lucifer/reference/SMC.md))
on a simplified version of the model provides a valuable sanity check.

  

## References

**\[1\]** Marin, J.-M., Pudlo, P., Robert, C.P., and Ryder, R.J. (2012).
*Approximate Bayesian computational methods*. Statistics and Computing,
22(6), 1167-1180. [DOI:
10.1007/s11222-011-9288-2](https://doi.org/10.1007/s11222-011-9288-2)

**\[2\]** Marjoram, P., Molitor, J., Plagnol, V., and Tavare, S. (2003).
*Markov chain Monte Carlo without likelihoods*. Proceedings of the
National Academy of Sciences, 100(26), 15324-15328. [DOI:
10.1073/pnas.0306899100](https://doi.org/10.1073/pnas.0306899100)

**\[3\]** Sisson, S.A., Fan, Y., and Tanaka, M.M. (2007). *Sequential
Monte Carlo without likelihoods*. Proceedings of the National Academy of
Sciences, 104(6), 1760-1765. [DOI:
10.1073/pnas.0607208104](https://doi.org/10.1073/pnas.0607208104)

**\[4\]** Beaumont, M.A., Cornuet, J.-M., Marin, J.-M., and Robert, C.P.
(2009). *Adaptive approximate Bayesian computation*. Biometrika, 96(4),
983-990. [DOI:
10.1093/biomet/asp052](https://doi.org/10.1093/biomet/asp052)

**\[5\]** Roberts, G.O., Gelman, A., and Gilks, W.R. (1997). *Weak
convergence and optimal scaling of random walk Metropolis algorithms*.
Annals of Applied Probability, 7(1), 110-120. [DOI:
10.1214/aoap/1034625254](https://doi.org/10.1214/aoap/1034625254)

**\[6\]** Jiang, B., Wu, T.-Y., Zheng, C., and Wong, W.H. (2017).
*Learning summary statistic for approximate Bayesian computation via
deep neural network*. Statistica Sinica, 27(4), 1595-1618. [DOI:
10.5705/ss.202015.0340](https://doi.org/10.5705/ss.202015.0340)

**\[7\]** Csillery, K., Blum, M.G.B., Gaggiotti, O.E., and Francois, O.
(2010). *Approximate Bayesian Computation (ABC) in practice*. Trends in
Ecology and Evolution, 25(7), 410-418. [DOI:
10.1016/j.tree.2010.04.001](https://doi.org/10.1016/j.tree.2010.04.001)

**\[8\]** Sisson, S.A., Fan, Y., and Beaumont, M.A. (2018). *Handbook of
Approximate Bayesian Computation*. Chapman and Hall/CRC. [DOI:
10.1201/9781315117195](https://doi.org/10.1201/9781315117195)

**\[9\]** Albert, C., Kuensch, H.R., and Scheidegger, A. (2015). *A
simulated annealing approach to approximate Bayes computations*.
Statistics and Computing, 25(6), 1217-1232. [DOI:
10.1007/s11222-014-9507-8](https://doi.org/10.1007/s11222-014-9507-8)

**\[10\]** Kirkpatrick, S., Gelatt, C.D., and Vecchi, M.P. (1983).
*Optimization by simulated annealing*. Science, 220(4598), 671-680.
[DOI:
10.1126/science.220.4598.671](https://doi.org/10.1126/science.220.4598.671)

**\[11\]** Marinari, E. and Parisi, G. (1992). *Simulated tempering: a
new Monte Carlo scheme*. Europhysics Letters, 19(6), 451-458. [DOI:
10.1209/0295-5075/19/6/002](https://doi.org/10.1209/0295-5075/19/6/002)

**\[12\]** Geyer, C.J. and Thompson, E.A. (1995). *Annealing Markov
chain Monte Carlo with applications to ancestral inference*. Journal of
the American Statistical Association, 90(431), 909-920. [DOI:
10.1080/01621459.1995.10476590](https://doi.org/10.1080/01621459.1995.10476590)

**\[13\]** Wang, F. and Landau, D.P. (2001). *Efficient, multiple-range
random walk algorithm to calculate the density of states*. Physical
Review Letters, 86(10), 2050-2053. [DOI:
10.1103/PhysRevLett.86.2050](https://doi.org/10.1103/PhysRevLett.86.2050)

**\[14\]** Jabot, F., Faure, T., Dumoulin, N., and Albert, C. (2015).
*EasyABC: Efficient approximate Bayesian computation sampling schemes*.
R package version 1.5. [CRAN:
EasyABC](https://cran.r-project.org/package=EasyABC)
