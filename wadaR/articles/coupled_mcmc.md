# Coupled MCMC for unbiased estimation and convergence diagnostics

## Introduction

Standard Markov chain Monte Carlo produces estimators that are
asymptotically unbiased but carry a transient bias introduced by the
arbitrary choice of initial state. The chain must “forget” its starting
point before its samples can be treated as draws from the target
distribution, and this forgetting period, commonly called burn-in, is
selected by ad hoc inspection of trace plots or summary diagnostics. The
resulting estimates are therefore biased for any finite number of
iterations, and the magnitude of that bias is difficult to quantify.
This problem is not merely academic: it undermines the theoretical
guarantees that make MCMC attractive, and it complicates the aggregation
of results from parallel chains that may have been initialized
differently or run for different lengths.

Jacob, O’Leary, and Atchade [\[1\]](#ref1) showed that this bias can be
eliminated entirely by running two copies of the same MCMC kernel with a
maximal coupling. The key insight is that when two chains are coupled so
that they meet after a random time \\\tau\\ and remain identical
thereafter, a telescoping sum formula produces estimators whose
expectation is exactly the target posterior mean, regardless of how many
burn-in iterations are discarded. The meeting time \\\tau\\ is a random
variable whose distribution encodes the mixing properties of the
underlying kernel; once \\\tau\\ is finite with probability one, the
debiasing estimator is unbiased for any choice of the tuning parameter
\\k\\. This transforms the burn-in problem from an unverifiable modeling
assumption into a tunable bias-variance tradeoff with rigorous
theoretical backing.

The practical implications are substantial. Because each coupled pair of
chains produces an independent, unbiased estimate, the problem of
posterior inference becomes embarrassingly parallel: one can launch
thousands of short coupled runs on separate cores and average the
results, obtaining valid confidence intervals without ever diagnosing
convergence of a single long chain. The distribution of meeting times
provides a rigorous, tuning-free convergence diagnostic that
characterizes how quickly the kernel mixes, complementing or replacing
the heuristic diagnostics (trace plots, \\\hat{R}\\, ESS) that are
standard practice. The
[`coupled_mcmc()`](https://robustecologies.github.io/lucifer/reference/coupled_mcmc.md)
function in **lucifer** implements this framework with
reflection-maximal coupling of Gaussian random-walk Metropolis
proposals.

## The coupling framework

### Coupled Metropolis-Hastings chains

The algorithm runs two Metropolis-Hastings chains \\X_t\\ and \\Y_t\\
that share the same proposal kernel but are coupled through their
proposal mechanism. Chain \\Y\\ is initialized one iteration behind
chain \\X\\, so that at iteration \\t\\ the algorithm has \\X_t\\ and
\\Y\_{t-1}\\. The coupling is constructed so that the proposals for
\\X\\ and \\Y\\ agree with the highest possible probability at each
step, while each marginal chain remains a valid Metropolis-Hastings
chain targeting the posterior \\\pi(\theta) = p(\theta \mid y)\\.

Given current states \\x\\ (for chain \\X\\) and \\y\\ (for chain \\Y\\)
and a shared proposal covariance \\\Sigma\\ with Cholesky factor \\L\\
(so \\\Sigma = LL^\top\\), the coupled proposal proceeds in two stages.
First, draw a standard normal innovation \\z \sim \mathcal{N}(0, I_d)\\
and form the proposal for chain \\X\\ as \\z_X = x + Lz\\. Then evaluate
the coupling probability

\\ \log p = \log \mathcal{N}(z_X; y, \Sigma) - \log \mathcal{N}(z_X; x,
\Sigma) \\

which measures how likely the proposal \\z_X\\ would have been if it had
been generated from the kernel centered at \\y\\ instead of \\x\\. With
probability \\\min(1, \exp(\log p))\\, the maximal coupling succeeds and
the same proposal \\z_X\\ is used for both chains. If the coupling
fails, a reflection of the standardized innovation across the hyperplane
separating \\x\\ and \\y\\ generates the residual proposal for chain
\\Y\\: the unit direction \\e = (x - y) / \\x - y\\\\ defines this
hyperplane, and the reflected innovation is \\z' = z - 2(e^\top z) e\\,
producing the proposal \\z_Y = y + Lz'\\. This reflection-maximal
coupling, introduced in the context of MCMC by [\[1\]](#ref1), maximizes
the probability that both proposals agree while ensuring that the
residual (when they disagree) is a valid draw from the correct
conditional distribution.

### Meeting and coalescence

Each chain independently accepts or rejects its proposal via the
standard Metropolis-Hastings acceptance ratio, so the coupling of
proposals does not guarantee that both chains accept or reject
simultaneously. Nevertheless, when the chains are close in parameter
space, the coupling probability \\p\\ approaches 1 and both proposals
are likely to be identical; if both chains then accept the same
proposal, they meet. Formally, the meeting time is \\\tau = \min\\t \geq
1 : X_t = Y\_{t-1}\\\\, detected in practice when \\\\X_t -
Y\_{t-1}\\\_\infty \< \epsilon\\ for a small tolerance \\\epsilon\\
(default \\10^{-8}\\). After meeting, the chains are set equal and
evolve identically using a single Markov kernel, so they remain coupled
forever. The existence of a finite meeting time with probability one is
guaranteed under mild conditions on the kernel, specifically that the
kernel is geometrically ergodic, which holds for random-walk Metropolis
on compact support or with sufficiently light-tailed targets
[\[1\]](#ref1).

## The debiasing estimator

The central result of [\[1\]](#ref1) is the construction of unbiased
estimators from coupled chains. For any square-integrable function \\h\\
of the parameters, define the single-lag estimator

\\ H_k = h(X_k) + \sum\_{t=k+1}^{\tau - 1} \left\[ h(X_t) - h(Y\_{t-1})
\right\] \\

where the sum is empty (and thus zero) when \\\tau \leq k + 1\\. By
Theorem 1 of [\[1\]](#ref1), \\E\[H_k\] = E\_\pi\[h(\theta)\]\\ exactly,
not asymptotically, provided the meeting time has finite expectation.
The estimator works by a telescoping argument: after the chains meet at
time \\\tau\\, the correction terms \\h(X_t) - h(Y\_{t-1})\\ vanish
because \\X_t = Y\_{t-1}\\ for all \\t \geq \tau\\, so the sum truncates
at \\\tau - 1\\. The terms before meeting form a telescoping series that
cancels the initialization bias in \\h(X_k)\\.

For variance reduction, the time-averaged version

\\ H\_{k:m} = \frac{1}{m - k + 1} \sum\_{l=k}^{m} H_l \\

averages the estimator over a window of iterations from \\k\\ to \\m\\,
where \\m \geq \max(k, \tau - 1)\\. This is analogous to using more
post-burn-in samples in standard MCMC but retains the unbiasedness
property. The parameter \\k\\ controls a bias-variance tradeoff: larger
\\k\\ discards more of the initial transient, reducing the magnitude of
the correction terms (and hence the variance of \\H_k\\), but it also
increases the minimum computational cost per coupled pair because the
chains must run for at least \\k\\ iterations. In practice, \\k\\ is
chosen as the median or upper quantile of the meeting-time distribution
from pilot runs.

## Usage in lucifer

### A simple normal model

The
[`coupled_mcmc()`](https://robustecologies.github.io/lucifer/reference/coupled_mcmc.md)
function shares the same model specification interface as
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
and the other inference engines in the package. The model function must
return a list with `LP` (log-posterior), `Dev` (deviance), `Monitor`,
`yhat`, and `parm`. The main coupling loop runs in a C++ backend
(`coupled_mcmc_cpp` in `src/sampler_other.cpp`), which handles the
reflection-maximal coupling, MH accept/reject for both chains, meeting
detection, and post-meeting extension natively; the R wrapper handles
input validation, Cholesky decomposition, and S3 class construction. The
following example uses a normal likelihood with conjugate-like priors on
the mean and log-transformed standard deviation.

``` r

library(lucifer)

# Simulate data
set.seed(42)
N <- 50
y <- rnorm(N, mean = 3, sd = 1.5)

Data <- list(
  N = N, y = y,
  mon.names = "LP",
  parm.names = c("mu", "log.sigma")
)

Model <- function(parm, Data) {
  mu    <- parm[1]
  sigma <- exp(parm[2])
  LL <- sum(dnorm(Data$y, mu, sigma, log = TRUE))
  LP <- LL + dnorm(mu, 0, 100, log = TRUE) +
    dcauchy(sigma, 0, 25, log = TRUE)
  yhat    <- rep(mu, Data$N)
  Monitor <- LP
  list(LP = LP, Dev = -2 * LL, Monitor = Monitor,
       yhat = yhat, parm = parm)
}
```

### Running a single coupled pair

The simplest call runs a single coupled pair with a specified burn-in
parameter `k`. The algorithm draws proposals from a Gaussian random walk
with covariance \\(2.38)^2 \Sigma / d\\, following the optimal scaling
result of Roberts, Gelman, and Gilks [\[3\]](#ref3), and uses
reflection-maximal coupling to maximize the probability that both chains
receive the same proposal.

``` r

set.seed(123)
fit <- coupled_mcmc(Model, Data, Initial.Values = c(0, 0),
                    max_iterations = 5000, k = 100, m = 500,
                    verbose = TRUE)
```

The `print` method provides a concise one-screen overview of the run,
including whether the chains met, the meeting time, chain dimensions,
and wall-clock time.

``` r

print(fit)
```

The `summary` method extends this with the naive posterior summary from
chain \\X\\ (discarding the first \\k\\ iterations as burn-in) and the
debiased posterior mean estimates computed via the \\H\_{k:m}\\ formula.

``` r

summary(fit)
```

### Trace plots and log-posterior diagnostics

The default `plot` method produces trace plots of both chains with the
meeting time marked by a dashed vertical line. Before meeting, the two
chains explore the parameter space independently; after meeting, they
collapse onto a single trajectory, confirming that the coupling
succeeded.

``` r

plot(fit, type = "trace")
```

The `type = "log_posterior"` option shows the log-posterior traces of
both chains, which is useful for diagnosing whether the chains reached
comparable posterior density levels before meeting.

``` r

plot(fit, type = "log_posterior")
```

## Meeting time diagnostics

### The meeting-time distribution

The distribution of meeting times \\\tau\\ across independent replicate
runs provides a powerful convergence diagnostic that requires no
subjective tuning. The tail behavior of this distribution directly
reflects the mixing properties of the underlying MCMC kernel: a kernel
that mixes rapidly will produce meeting times concentrated at small
values with thin tails, while a slowly mixing kernel will produce
heavy-tailed meeting times. Unlike \\\hat{R}\\ or effective sample size,
which measure properties of a single chain’s output, the meeting-time
distribution characterizes the kernel itself and is invariant to the
length of the run.

To construct this diagnostic, one runs many independent coupled pairs
with the same model and data, collects the meeting times, and examines
their empirical distribution. The
[`plot_meeting_times()`](https://robustecologies.github.io/lucifer/reference/plot_meeting_times.md)
function produces a histogram annotated with the number of runs that
met, the median meeting time, and the maximum meeting time. Large right
tails or a substantial fraction of non-meeting runs signal that the
kernel needs better tuning, perhaps a different proposal covariance, a
reparameterized model, or a more efficient MCMC algorithm.

``` r

# Run 50 independent coupled pairs
set.seed(2025)
fits <- replicate(50, coupled_mcmc(
  Model, Data, Initial.Values = c(0, 0),
  max_iterations = 5000, k = 100
), simplify = FALSE)

# Plot the meeting-time distribution
plot_meeting_times(fits)
```

### Comparing kernels via meeting times

One of the most valuable applications of meeting-time diagnostics is
comparing the mixing efficiency of different MCMC kernels or tuning
configurations on the same problem. Because the meeting time \\\tau\\ is
a summary of kernel quality that does not depend on arbitrary choices
like burn-in or thinning, comparing the meeting-time distributions of
two kernels gives a rigorous answer to the question of which kernel
mixes faster. For instance, one could compare coupled runs using
different proposal covariance matrices, or compare the meeting times of
a random-walk Metropolis kernel against a gradient-based kernel applied
to the same target. The kernel with shorter meeting times and thinner
tails is mixing more efficiently, and the quantitative difference in
median meeting times provides a concrete metric for the improvement.

``` r

# Compare default vs. diagonal covariance
fits_default <- replicate(30, coupled_mcmc(
  Model, Data, Initial.Values = c(0, 0),
  max_iterations = 5000
), simplify = FALSE)

# Use a tighter proposal
Sigma_tight <- diag(c(0.1, 0.1))
fits_tight <- replicate(30, coupled_mcmc(
  Model, Data, Initial.Values = c(0, 0),
  Covar = Sigma_tight, max_iterations = 5000
), simplify = FALSE)

# Extract meeting times
tau_default <- sapply(fits_default, function(f) f$meeting_time)
tau_tight   <- sapply(fits_tight, function(f) f$meeting_time)

knitr::kable(data.frame(
  Proposal = c("Default (2.38^2/d * I)", "Tight (diag(0.1, 0.1))"),
  `Median tau` = c(median(tau_default, na.rm = TRUE),
                   median(tau_tight, na.rm = TRUE)),
  check.names = FALSE
), caption = "Meeting times under different proposal covariances.")
```

## Practical recommendations

The choice of \\k\\ controls the bias-variance tradeoff in the debiasing
estimator. When \\k\\ is small, the correction terms \\h(X_t) -
h(Y\_{t-1})\\ may be large because the chains have not yet approached
stationarity, inflating the variance of \\H_k\\. When \\k\\ is large,
these correction terms are typically small (because the chains are
near-stationary and close to meeting), but the computational cost per
pair increases because each run must last at least \\k\\ iterations. A
practical strategy is to run a pilot batch of 50 to 100 coupled pairs
with \\k = 0\\, examine the meeting-time distribution, and set \\k\\
equal to the median or the 75th percentile of the observed meeting times
for the production runs. The parameter \\m\\ should be set at least as
large as \\k\\; setting \\m \> k\\ provides additional post-burn-in
samples that reduce the variance of \\H\_{k:m}\\ at a modest additional
cost per pair.

For production-quality unbiased estimation, the workflow is to run a
large number \\R\\ of independent coupled pairs (each with its own
random seed), compute \\H\_{k:m}^{(r)}\\ for each pair \\r = 1, \ldots,
R\\, and report the average \\\bar{H} = R^{-1} \sum\_{r=1}^R
H\_{k:m}^{(r)}\\ as the point estimate. Standard errors and confidence
intervals follow from the central limit theorem applied to the i.i.d.
estimates \\H\_{k:m}^{(r)}\\. This procedure is embarrassingly parallel:
each coupled pair runs independently, with no communication between
pairs, making it straightforward to distribute across cores or machines.
The per-pair cost is roughly twice that of a single MCMC run of length
\\\max(m, \tau)\\, but the parallelism and the elimination of burn-in
bias often make this competitive with long single-chain runs,
particularly in settings where convergence of a single chain is
difficult to assess.

``` r

# Production run: 200 coupled pairs, k = median pilot tau
R <- 200
k_prod <- median(sapply(fits, function(f) f$meeting_time), na.rm = TRUE)

set.seed(999)
production <- replicate(R, coupled_mcmc(
  Model, Data, Initial.Values = c(0, 0),
  max_iterations = 10000,
  k = k_prod, m = 2 * k_prod
), simplify = FALSE)

# Unbiased posterior mean estimates with standard errors
estimates <- t(sapply(production, function(f) {
  if (f$met) summary(f)  # triggers .debiased_estimate internally
  # For manual extraction:
  lucifer:::.debiased_estimate(f)
}))

knitr::kable(data.frame(
  Parameter = c("mu", "log_sigma"),
  `Unbiased estimate` = c(mean(estimates[, 1], na.rm = TRUE),
                          mean(estimates[, 2], na.rm = TRUE)),
  `Standard error` = c(sd(estimates[, 1], na.rm = TRUE) / sqrt(R),
                        sd(estimates[, 2], na.rm = TRUE) / sqrt(R)),
  check.names = FALSE
), caption = "Debiased posterior mean estimates from 200 coupled pairs.",
  digits = 4)
```

## References

**\[1\]** Jacob, P. E., O’Leary, J., and Atchade, Y. F. (2020).
*Unbiased Markov chain Monte Carlo methods with couplings*. Journal of
the Royal Statistical Society: Series B, 82(3), 543-600. [DOI:
10.1111/rssb.12336](https://doi.org/10.1111/rssb.12336)

**\[2\]** Heng, J. and Jacob, P. E. (2019). *Unbiased Hamiltonian Monte
Carlo with couplings*. arXiv preprint arXiv:1909.13339. [arXiv:
1909.13339](https://arxiv.org/abs/1909.13339)

**\[3\]** Roberts, G. O., Gelman, A., and Gilks, W. R. (1997). *Weak
convergence and optimal scaling of random walk Metropolis algorithms*.
Annals of Applied Probability, 7(1), 110-120. [DOI:
10.1214/aoap/1034625254](https://doi.org/10.1214/aoap/1034625254)
