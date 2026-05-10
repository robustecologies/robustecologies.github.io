# Approximate Bayesian Computation

Performs likelihood-free Bayesian inference via Approximate Bayesian
Computation (ABC). Unlike the other estimation functions in lucifer, the
model function here acts as a forward simulator that generates synthetic
data rather than evaluating a log-posterior. Acceptance is determined by
whether the Euclidean distance between summary statistics of simulated
and observed data falls below a tolerance \\\varepsilon\\, yielding
samples from the approximate posterior \\p\_\varepsilon(\theta \mid
S(y\_\text{obs})) \propto p(\theta)\\\Pr(\rho(S(y^\*), S(y\_\text{obs}))
\le \varepsilon)\\. Four algorithms are available: rejection sampling,
MCMC-ABC, SMC-ABC, and Simulated Annealing ABC (SA-ABC).

## Usage

``` r
ABC(
  Model,
  Data,
  Prior,
  Summary.Stats,
  Observed.Stats,
  Method = "rejection",
  N = 10000,
  epsilon = NULL,
  Covar = NULL,
  MCMC.iterations = 1000,
  Status = 100,
  SA.iterations = 5000,
  T0 = NULL,
  cooling = "geometric",
  cooling.rate = NULL,
  CPUs = 1L
)
```

## Arguments

- Model:

  A function with signature `function(parm, Data)` that returns a
  simulated dataset as a numeric vector or matrix. This is a forward
  simulator, not a likelihood evaluator: given parameter values `parm`
  and auxiliary information in `Data`, it generates synthetic
  observations from the model. For proposals with invalid parameter
  values (e.g., negative variance), the function should return extreme
  or degenerate data rather than throwing an error, so that the distance
  criterion naturally rejects them.

- Data:

  A list containing at least `parm.names`, a character vector of
  parameter names whose length defines the dimensionality of the
  inference problem. Additional elements (sample size, design matrices,
  hyperparameters, etc.) are passed through to `Model` and may be
  accessed within it.

- Prior:

  A function with no arguments that returns a single draw from the joint
  prior distribution as a numeric vector of length equal to
  `length(Data$parm.names)`. The prior must have support that covers the
  region of non-negligible posterior mass; for constrained parameters
  (e.g., variances), enforce positivity within this function.

- Summary.Stats:

  A function with signature `function(data)` that computes summary
  statistics from a dataset, returning a numeric vector. The choice of
  summary statistics is the most consequential design decision in ABC:
  sufficient statistics (when they exist) yield exact inference as
  \\\varepsilon \to 0\\; too few statistics lose information and widen
  the posterior; too many increase the effective comparison dimension
  and drive down acceptance rates. When statistics have different
  scales, standardize them within this function to prevent any single
  statistic from dominating the Euclidean distance.

- Observed.Stats:

  A numeric vector of pre-computed summary statistics for the observed
  data, obtained by applying `Summary.Stats()` to the observed dataset.
  Pre-computing avoids redundant evaluation at every iteration.

- Method:

  Character string selecting the ABC algorithm. One of `"rejection"`
  (default), `"MCMC"`, `"SMC"`, or `"SA"`. Partial matching is
  supported. See **Details** for algorithm descriptions and guidance on
  when to use each method.

- N:

  Integer controlling the sampling effort. For rejection ABC, this is
  the target number of accepted samples; the algorithm draws up to
  \\1000 \times N\\ proposals to reach this target. For MCMC-ABC and
  SMC-ABC, `N` is passed through to control the total proposal budget
  and particle population size respectively. Default 10000.

- epsilon:

  Numeric tolerance threshold. Proposals are accepted when \\\\S(y^\*) -
  S(y\_\text{obs})\\\_2 \le \varepsilon\\. If `NULL` (default), an
  adaptive tolerance is computed from a pilot run: \\\min(1000, N)\\
  prior-predictive simulations are drawn and \\\varepsilon\\ is set to
  the 1st percentile of the resulting distance distribution (the 5th
  percentile if the 1st is zero). Lower values improve posterior
  accuracy at the cost of acceptance rate.

- Covar:

  Proposal covariance matrix for MCMC-ABC, a \\J \times J\\
  positive-definite matrix where \\J\\ is the number of parameters.
  Scaled internally by \\2.38^2 / J\\ following Roberts et al. (1997).
  Defaults to `NULL`, which uses the identity matrix. For SMC-ABC, the
  rejuvenation kernel uses twice the empirical covariance of the current
  particle population, ignoring this argument.

- MCMC.iterations:

  Integer number of MCMC iterations for `Method = "MCMC"`. The chain
  includes a burn-in phase (typically the first 20% of iterations should
  be discarded before computing posterior summaries). Default 1000.

- Status:

  Integer controlling the frequency of progress messages. A status line
  is printed every `Status` proposals (rejection ABC), iterations
  (MCMC-ABC and SA-ABC), or at scaled intervals (SMC-ABC). Default 100.

- SA.iterations:

  Integer number of SA-ABC iterations for `Method = "SA"`. Default 5000.
  The first 30% of iterations are used for annealing (Phase 1) and the
  remaining 70% for hard-threshold sampling (Phase 2). Only Phase 2
  samples are returned in the `Posterior` field.

- T0:

  Numeric initial temperature for SA-ABC. If `NULL` (default),
  auto-calibrated from the median pilot distance so that approximately
  37% of uphill moves are accepted at the start.

- cooling:

  Character string specifying the cooling schedule for SA-ABC:
  `"geometric"` (default) or `"linear"`. Geometric cooling uses \\T_t =
  T_0 \cdot r^t\\ and linear uses \\T_t = T_0 (1 - t / N)\\.

- cooling.rate:

  Numeric cooling rate for geometric SA-ABC cooling, in \\(0, 1)\\. If
  `NULL` (default), auto-calibrated so that the final temperature equals
  \\\varepsilon\\: \\r = (\varepsilon / T_0)^{1/N}\\, clamped to
  \\\[0.9, 0.9999\]\\.

- CPUs:

  Integer number of parallel workers for rejection ABC. When `CPUs > 1`
  and `Method = "rejection"`, the target sample count is split across
  workers running independent rejection loops. Uses PSOCK clusters
  (works on all platforms). Default 1 (sequential).

## Value

An object of class `abc` with S3 methods
[`print.abc`](https://robustecologies.github.io/lucifer/reference/print.abc.md),
[`summary.abc`](https://robustecologies.github.io/lucifer/reference/summary.abc.md),
and
[`plot.abc`](https://robustecologies.github.io/lucifer/reference/plot.abc.md).
The object is a list containing:

- Posterior:

  Numeric matrix of accepted parameter samples. Column names correspond
  to `Data$parm.names`. For rejection and SMC-ABC, all accepted/final
  particles. For MCMC-ABC, the full chain (burn-in should be discarded
  manually). For SA-ABC, only the Phase 2 (sampling phase) samples.

- Distances:

  Numeric vector of Euclidean distances \\\\S(y^\*) -
  S(y\_\text{obs})\\\_2\\ for each accepted sample (rejection ABC) or
  each iteration (MCMC/SMC-ABC).

- Acceptance.rate:

  Scalar acceptance rate: number of accepted proposals divided by total
  proposals for rejection ABC, accepted moves divided by chain length
  for MCMC-ABC, or the average per-stage rejuvenation acceptance rate
  for SMC-ABC.

- Epsilon:

  The final tolerance \\\varepsilon\\ used, either user-specified or
  adaptively determined from the pilot run.

- Method:

  Character string identifying the algorithm used.

- Summary:

  Numeric matrix with rows for each parameter and columns `Mean`, `SD`,
  `2.5%`, `50%`, `97.5%`, computed from the accepted posterior samples.

- N.total:

  Total number of proposals or iterations executed.

- Temperatures:

  (SA-ABC only) Numeric vector of length `SA.iterations` giving the
  temperature at each iteration. `NULL` for other methods.

- Interrupted:

  Logical flag indicating whether the computation was interrupted by the
  user before completing. If `TRUE`, the `Posterior`, `Distances`, and
  other fields contain partial results up to the point of interruption.

- Minutes:

  Wall-clock time in minutes.

- call:

  The matched call.

## Details

Approximate Bayesian Computation

### Rejection ABC

The simplest and most parallelizable algorithm. For each proposal: draw
\\\theta^\* \sim p(\theta)\\ from the prior, simulate \\y^\* =
\text{Model}(\theta^\*, \text{Data})\\, compute summary statistics
\\S(y^\*)\\, and accept \\\theta^\*\\ if \\\\S(y^\*) -
S(y\_\text{obs})\\\_2 \le \varepsilon\\. The accepted samples are
independent draws from the approximate posterior \\p\_\varepsilon(\theta
\mid S(y\_\text{obs}))\\.

The computational cost is \\N / r\\ model simulations where \\r\\ is the
acceptance rate. For low-dimensional problems (\\J \le 3\\) with
informative priors, rejection ABC is often adequate. Acceptance rates
between 0.001 and 0.01 are typical; rates below 0.0001 suggest that
\\\varepsilon\\ is too tight or the summary statistics are too
high-dimensional.

### MCMC-ABC

The algorithm of Marjoram et al. (2003) replaces independent prior
sampling with a Markov chain that concentrates proposals near the
current accepted state. At each iteration, a candidate \\\theta^\*\\ is
drawn from a random-walk Gaussian kernel centered at the current
\\\theta\\ with covariance \\(2.38^2 / J) \\ \Sigma\\, data are
simulated, and the proposal is accepted if the distance criterion is
satisfied:

\$\$\alpha(\theta, \theta^\*) = \begin{cases} 1 & \text{if }
\rho(S(y^\*), S(y\_\text{obs})) \le \varepsilon \\ 0 & \text{otherwise.}
\end{cases}\$\$

The chain must be initialized at a point satisfying the distance
criterion. The implementation searches pilot draws first, then draws up
to 50000 additional prior samples; if no valid starting point is found,
it uses the closest draw and widens \\\varepsilon\\.

MCMC-ABC produces correlated samples and requires burn-in removal.
Acceptance rates between 0.05 and 0.3 are healthy. The advantage over
rejection ABC is that proposals cluster near the posterior, dramatically
improving the acceptance rate for concentrated posteriors.

### SMC-ABC

Sequential Monte Carlo ABC (Sisson et al., 2007; Beaumont et al., 2009)
evolves a particle population through a sequence of decreasing
tolerances \\\varepsilon_1 \> \cdots \> \varepsilon_T\\:

1.  **Initialization.** Draw \\N\\ particles from the prior and compute
    their distances.

2.  **Tolerance schedule.** If `epsilon` is `NULL`, construct a 10-stage
    schedule from the 50th to the 1st percentile of the initial distance
    distribution; otherwise interpolate linearly from the median to the
    target.

3.  **For each stage:** reweight particles (weight 1 if within
    tolerance, 0 otherwise), resample multinomially, and rejuvenate via
    one step of MCMC-ABC at the current tolerance using twice the
    empirical covariance as the proposal kernel. Terminate early if
    fewer than 2 particles satisfy the tolerance.

SMC-ABC avoids the burn-in problem of MCMC-ABC, adapts the tolerance
schedule automatically, and provides a diverse particle population.
Per-stage acceptance rates during rejuvenation should be above 0.1.

### SA-ABC (Simulated Annealing)

SA-ABC replaces the hard distance threshold of MCMC-ABC with a soft,
temperature-modulated Metropolis criterion inspired by simulated
annealing (Albert et al., 2015). At each iteration, a candidate
\\\theta^\*\\ is proposed from a random-walk Gaussian kernel as in
MCMC-ABC, data are simulated, and the acceptance probability is:

\$\$\alpha(\theta, \theta^\*) = \min\\\Big(1,\\
\exp\\\Big(-\frac{d(\theta^\*) - d(\theta)}{T_t}\Big)\Big)\$\$

where \\d(\cdot)\\ is the Euclidean distance between simulated and
observed summary statistics and \\T_t\\ is the temperature at iteration
\\t\\. When the proposed distance is smaller (downhill), the move is
always accepted; when it is larger (uphill), acceptance probability
decreases as temperature drops.

The algorithm operates in two phases. Phase 1 (annealing, first 30% of
iterations) cools the temperature from \\T_0\\ to \\T\_\text{final} =
\varepsilon\\, driving the chain toward the low-distance region via the
soft acceptance criterion. The proposal kernel is adapted continuously
during annealing. Phase 2 (sampling, remaining 70%) switches to the hard
ABC threshold \\d(\theta^\*) \le \varepsilon\\, identical to MCMC-ABC,
so the chain targets the standard ABC posterior \\p\_\varepsilon(\theta
\mid S(y\_\text{obs}))\\. At the phase boundary, the chain is
reinitialized at a point satisfying \\d \le \varepsilon\\ (from
annealing samples within the tolerance, or via rejection from the prior
if none exist), and the proposal kernel is adapted from those annealing
samples that fell within the tolerance region.

Two cooling schedules are available for Phase 1. Geometric cooling sets
\\T_t = T_0 \cdot r^t\\ where \\r \in (0, 1)\\ is the cooling rate.
Linear cooling interpolates linearly from \\T_0\\ to \\\varepsilon\\. By
default, \\T_0\\ is set to the median pilot distance (ensuring roughly
37% uphill acceptance at the start) and the cooling rate is calibrated
so the temperature reaches \\\varepsilon\\ at the end of Phase 1.

The posterior summary is computed from Phase 2 samples only (the second
half of the chain). SA-ABC is particularly useful when the distance
landscape has multiple local optima, since the annealing phase allows
the chain to escape local traps before the hard- threshold sampling
phase begins.

This implementation differs from the full SABC algorithm of Albert et
al. (2015), which uses a particle ensemble with adaptive cooling based
on minimal entropy production, CDF-reparametrized distances, and a soft
Boltzmann kernel throughout (as implemented in EasyABC). The lucifer
variant uses a single chain with a fixed cooling schedule and switches
to the hard indicator kernel in Phase 2, so that all four ABC methods
target the same posterior \\p\_\varepsilon(\theta \mid
S(y\_\text{obs}))\\. This design requires only a prior sampler (not an
evaluable prior density) and produces directly comparable output across
methods. See
[`vignette("abc")`](https://robustecologies.github.io/lucifer/articles/abc.md)
for a detailed discussion of the design rationale and comparison with
alternative approaches.

### Choosing summary statistics

The quality of the ABC approximation depends critically on the summary
statistics. Sufficient statistics yield exact inference as \\\varepsilon
\to 0\\. When sufficient statistics are unavailable, a moderate number
(5 to 15) of scientifically motivated summaries (moments, quantiles,
autocorrelations, spectral features) is recommended. See Csillery et al.
(2010) for a practical review and the
[`vignette("abc")`](https://robustecologies.github.io/lucifer/articles/abc.md)
for worked examples.

## References

Marin, J.-M., Pudlo, P., Robert, C.P., and Ryder, R.J. (2012).
"Approximate Bayesian computational methods." *Statistics and
Computing*, 22(6), 1167–1180.
[doi:10.1007/s11222-011-9288-2](https://doi.org/10.1007/s11222-011-9288-2)

Marjoram, P., Molitor, J., Plagnol, V., and Tavare, S. (2003). "Markov
chain Monte Carlo without likelihoods." *Proceedings of the National
Academy of Sciences*, 100(26), 15324–15328.
[doi:10.1073/pnas.0306899100](https://doi.org/10.1073/pnas.0306899100)

Sisson, S.A., Fan, Y., and Tanaka, M.M. (2007). "Sequential Monte Carlo
without likelihoods." *Proceedings of the National Academy of Sciences*,
104(6), 1760–1765.
[doi:10.1073/pnas.0607208104](https://doi.org/10.1073/pnas.0607208104)

Beaumont, M.A., Cornuet, J.-M., Marin, J.-M., and Robert, C.P. (2009).
"Adaptive approximate Bayesian computation." *Biometrika*, 96(4),
983–990.
[doi:10.1093/biomet/asp052](https://doi.org/10.1093/biomet/asp052)

Roberts, G.O., Gelman, A., and Gilks, W.R. (1997). "Weak convergence and
optimal scaling of random walk Metropolis algorithms." *Annals of
Applied Probability*, 7(1), 110–120.
[doi:10.1214/aoap/1034625254](https://doi.org/10.1214/aoap/1034625254)

Albert, C., Kuensch, H.R., and Scheidegger, A. (2015). "A simulated
annealing approach to approximate Bayes computations." *Statistics and
Computing*, 25(6), 1217–1232.
[doi:10.1007/s11222-014-9507-8](https://doi.org/10.1007/s11222-014-9507-8)

Csillery, K., Blum, M.G.B., Gaggiotti, O.E., and Francois, O. (2010).
"Approximate Bayesian Computation (ABC) in practice." *Trends in Ecology
and Evolution*, 25(7), 410–418.
[doi:10.1016/j.tree.2010.04.001](https://doi.org/10.1016/j.tree.2010.04.001)

Sisson, S.A., Fan, Y., and Beaumont, M.A. (2018). *Handbook of
Approximate Bayesian Computation*. Chapman and Hall/CRC.
[doi:10.1201/9781315117195](https://doi.org/10.1201/9781315117195)

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
for likelihood-based MCMC inference;
[`SMC`](https://robustecologies.github.io/lucifer/reference/SMC.md) for
likelihood-based Sequential Monte Carlo;
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md)
for deterministic variational approximation;
[`vignette("abc")`](https://robustecologies.github.io/lucifer/articles/abc.md)
for worked examples including the MA(2) intractable-likelihood problem.

## Examples

``` r
if (FALSE) { # \dontrun{
# --- Normal mean-variance inference via ABC ---

# Forward simulator (returns synthetic data, not a log-posterior)
Simulator <- function(parm, Data) {
  mu <- parm[1]; sigma <- parm[2]
  if (sigma <= 0) return(rep(1e10, Data$n))
  rnorm(Data$n, mu, sigma)
}

# Joint prior: mu ~ N(0, 10), sigma ~ Half-Cauchy(0, 5)
Prior <- function() c(rnorm(1, 0, 10), abs(rcauchy(1, 0, 5)))

# Sufficient statistics for the normal model
Summary.Stats <- function(data) c(mean(data), sd(data))

# Observed data
set.seed(42)
obs_data <- rnorm(100, mean = 3, sd = 1.5)
Observed.Stats <- Summary.Stats(obs_data)
Data <- list(n = 100, parm.names = c("mu", "sigma"))

# Rejection ABC
fit_rej <- ABC(Simulator, Data, Prior, Summary.Stats, Observed.Stats,
               Method = "rejection", N = 2000, Status = 10000)
print(fit_rej)
summary(fit_rej)
plot(fit_rej)
plot(fit_rej, type = "distances")

# MCMC-ABC with user-specified proposal covariance
fit_mcmc <- ABC(Simulator, Data, Prior, Summary.Stats, Observed.Stats,
                Method = "MCMC", MCMC.iterations = 10000,
                Covar = diag(c(0.3, 0.2)), Status = 5000)
print(fit_mcmc)
summary(fit_mcmc)
plot(fit_mcmc)

# SMC-ABC with adaptive tolerance schedule
fit_smc <- ABC(Simulator, Data, Prior, Summary.Stats, Observed.Stats,
               Method = "SMC", N = 1000, Status = 100)
print(fit_smc)
summary(fit_smc)
plot(fit_smc)

# SA-ABC with geometric cooling (auto-calibrated)
fit_sa <- ABC(Simulator, Data, Prior, Summary.Stats, Observed.Stats,
              Method = "SA", SA.iterations = 5000, Status = 2500)
print(fit_sa)
summary(fit_sa)
plot(fit_sa)
plot(fit_sa, type = "temperature")
} # }
```
