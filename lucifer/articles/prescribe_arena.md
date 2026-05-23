# Algorithm selection and universal benchmarking

## The algorithm selection problem

Bayesian computation has fragmented into a landscape of fundamentally
different inference strategies, each with distinct regimes of efficiency
and failure. Markov chain Monte Carlo methods produce asymptotically
exact samples but scale poorly in high dimensions without geometric
information; variational inference sacrifices accuracy for speed by
optimizing within a restricted family of distributions; Laplace
approximation provides rapid point estimates adequate only when the
posterior is approximately Gaussian near its mode; sequential Monte
Carlo tempers between prior and posterior to handle multimodality; and
simulation-based inference bypasses the likelihood entirely using neural
density estimators. The lucifer package implements over 130 methods
spanning all of these families: 82 MCMC algorithms across 13
subcategories (gradient, adaptive random-walk, componentwise, slice,
ensemble, QMC, geometric, PDMP, multimodal tempering,
constraint-handling, flow-enhanced, surrogate, and augmentation), 8
variational Bayes methods, 18 Laplace optimizers, 3 iterative quadrature
methods, population Monte Carlo, SMC, 4 ABC methods, and 6 SBI methods.

Selecting the right inference strategy for a given problem is a decision
that profoundly affects both the quality of the posterior approximation
and the computational cost of obtaining it. The classical approach is
trial and error, which imposes substantial analyst time and computation.
The literature on algorithm selection in optimization has formalized
this as the Algorithm Selection Problem [\[1\]](#ref1), but no
equivalent systematic framework exists for Bayesian inference beyond
informal rules of thumb.

[`Prescribe()`](https://robustecologies.github.io/lucifer/reference/Prescribe.md)
and
[`Arena()`](https://robustecologies.github.io/lucifer/reference/Arena.md)
address this gap.
[`Prescribe()`](https://robustecologies.github.io/lucifer/reference/Prescribe.md)
is a pre-fit advisor that profiles a Bayesian model without performing
inference, characterizes the problem along several computational axes
(dimensionality, evaluation speed, gradient availability, constraint
structure, posterior conditioning, and multimodality), and produces a
ranked recommendation across all inference families.
[`Arena()`](https://robustecologies.github.io/lucifer/reference/Arena.md)
is a post-fit benchmarking system that accepts heterogeneous fit objects
from any combination of inference methods and computes unified
efficiency, accuracy, and reliability metrics. Together they implement
the full workflow: profile the problem, choose a strategy, fit one or
more methods, and select the winner on empirical evidence.

## Theoretical foundations

### Model profiling and the no free lunch theorem

The no free lunch theorem for optimization [\[2\]](#ref2) establishes
that no algorithm uniformly dominates all others across all possible
problems. In the Bayesian context, this means that NUTS is not always
the best MCMC algorithm, variational inference is not always faster than
MCMC, and Laplace approximation is not always adequate for
low-dimensional problems. The optimal choice depends on problem-specific
properties that can be measured cheaply before committing to a full
inference run.

[`Prescribe()`](https://robustecologies.github.io/lucifer/reference/Prescribe.md)
implements a model profiling phase that evaluates the Model function
approximately \\n\_{\text{profile}} + 2d\\ times (where \\d\\ is the
number of parameters), which is typically negligible compared to the
thousands or millions of evaluations required for inference. The
profiling characterizes the problem along seven axes.

**Dimensionality.** The number of parameters \\d\\ is the single most
important determinant of algorithm performance. Gradient-based MCMC
methods (NUTS, HMC, MALA) exploit the geometry of the posterior and
scale as \\O(d^{5/4})\\ in mixing time [\[3\]](#ref3), vastly superior
to the \\O(d^2)\\ scaling of random-walk methods [\[4\]](#ref4).
Variational inference with mean-field approximations scales as \\O(d)\\
in per-iteration cost but suffers from increasingly severe approximation
error as correlations grow with dimensionality. Laplace approximation
becomes unreliable above approximately \\d = 30\\ because the Gaussian
assumption degrades exponentially in the tails. Iterative quadrature
suffers from the curse of dimensionality and is practically limited to
\\d \leq 10\\.

**Evaluation speed.** The time per model evaluation directly affects
which algorithms are feasible. Componentwise methods (AMWG, MWG, Gibbs,
Slice) evaluate the model \\d\\ times per iteration, making them
impractical when evaluations are slow and \\d\\ is moderate. Ensemble
methods (AIES, DEMC) require \\N\_\text{walkers}\\ evaluations per step.
Gradient-based methods require additional gradient evaluations (or
finite-difference approximations at \\O(d)\\ cost).
[`Prescribe()`](https://robustecologies.github.io/lucifer/reference/Prescribe.md)
measures this via timed evaluations at the initial values, reporting
evals/minute in four regimes: fast (\>1000), medium (100-1000), slow
(10-100), and very slow (\<10).

**Gradient availability.** User-supplied gradients enable the most
efficient MCMC algorithms. NUTS achieves near-optimal scaling by
simulating Hamiltonian dynamics using the gradient of the log-posterior
[\[3\]](#ref3). Without gradients, NUTS falls back to numerical
differentiation at \\O(d)\\ overhead per gradient evaluation, which may
still be worthwhile but changes the computational calculus.
[`Prescribe()`](https://robustecologies.github.io/lucifer/reference/Prescribe.md)
checks for `Data$user.grad` or `Data$gradient` and adjusts scores
accordingly.

**Constraint structure.** Parameters constrained to subsets of
\\\mathbb{R}^d\\ (e.g., positive variances, probability simplices)
affect algorithm choice in two ways. First, componentwise methods
naturally handle constraints by operating within each parameter’s
feasible region. Second, gradient-based methods require
reparameterization (log-transform for positivity, logit for bounded
parameters) to operate in unconstrained space, which introduces Jacobian
corrections and can distort the effective geometry.
[`Prescribe()`](https://robustecologies.github.io/lucifer/reference/Prescribe.md)
uses
[`is.constrained()`](https://robustecologies.github.io/lucifer/reference/is.constrained.md)
to detect the fraction of constrained parameters.

**Posterior conditioning.** The condition number of the local Hessian at
the initial values approximates the posterior’s aspect ratio.
Ill-conditioned posteriors (\\\kappa \> 1000\\) cause random-walk
methods to mix extremely slowly because the proposal must simultaneously
respect the smallest and largest scales. NUTS and other Hamiltonian
methods are more robust to ill-conditioning because they follow the
posterior’s curvature, but numerical integration becomes less stable.
[`Prescribe()`](https://robustecologies.github.io/lucifer/reference/Prescribe.md)
computes the Hessian via
[`numDeriv::hessian()`](https://rdrr.io/pkg/numDeriv/man/hessian.html)
when \\d \leq 100\\.

**Multimodality.** Multimodal posteriors are the bane of single-chain
MCMC, which can become trapped in a local mode for the entire run.
[`Prescribe()`](https://robustecologies.github.io/lucifer/reference/Prescribe.md)
estimates a multimodality score by evaluating the log-posterior at
\\n\_\text{profile}\\ random perturbations of the initial values and
comparing the distribution of LP values to what would be expected under
a unimodal Gaussian. High excess kurtosis or an unusually wide range
relative to the standard deviation suggests multiple modes.

**Posterior correlation structure.** When the Hessian is available (\\d
\leq 100\\),
[`Prescribe()`](https://robustecologies.github.io/lucifer/reference/Prescribe.md)
computes the maximum off-diagonal correlation from the local curvature
matrix. High correlation (\\\> 0.9\\) indicates funnel or banana-shaped
posteriors where mean-field variational approximations (ADVI.mf) fail
because they assume independence across parameters. Full-rank ADVI and
Riemannian methods (RMHMC, LMC) exploit correlated geometry and receive
a scoring bonus when this axis is elevated. The correlation score
appears as the seventh axis in the profile radar chart.

### Subcategory-aware scoring

Within MCMC,
[`Prescribe()`](https://robustecologies.github.io/lucifer/reference/Prescribe.md)
scores algorithms using a unified dispatch across 13 subcategories. Each
subcategory applies characteristic factor profiles. Three registry
metadata fields refine scoring beyond the problem profile. The quality
tier (1 = state-of-the-art, 2 = solid, 3 = niche) applies a
multiplicative bonus: tier-1 algorithms like NUTS, GIST, autoMALA, and
DREAM receive 1.25, while tier-3 algorithms like RWM receive 0.85. The
evaluation cost multiplier penalizes algorithms that require multiple
model evaluations per iteration: geometric methods (RMHMC, LMC, MHMC)
with a Hessian cost of 4 receive a \\1/\sqrt{4} = 0.5\\ penalty, while
PDMP methods (BPS, ZigZag) with amortized cost 0.5 receive a
\\1/\sqrt{0.5} = 1.41\\ bonus. The recommended dimension range applies a
soft penalty when the problem dimension falls outside an algorithm’s
effective operating range, decaying from 1.0 to 0.1 rather than hard
disqualification.

### Multiplicative scoring

The scoring phase assigns each inference strategy a composite score via
a geometric mean of multiplicative factors. Each factor \\f_k\\ measures
how well a particular problem axis matches the algorithm’s capabilities,
with \\f_k = 0\\ disqualifying the method entirely. For a method \\m\\
with factors \\\\f_1, \ldots, f_K\\\\, the composite score is

\\S(m) = \exp\left(\frac{1}{K}\sum\_{k=1}^K \log f_k\right) =
\left(\prod\_{k=1}^K f_k\right)^{1/K},\\

which equals the geometric mean and has the desirable property that a
single zero factor eliminates the method regardless of its performance
on other axes. This is appropriate because some problem-algorithm
mismatches are absolute: running ABC on a problem with an available
likelihood, or running gradient-based MCMC without gradients on a
discrete parameter, cannot be rescued by compensating strengths
elsewhere.

### Universal efficiency metrics

Comparing heterogeneous inference methods requires metrics that are
meaningful across fundamentally different computational strategies.
[`Arena()`](https://robustecologies.github.io/lucifer/reference/Arena.md)
computes three families of metrics.

**Efficiency: ESS per second.** The effective sample size (ESS) measures
the information content of correlated posterior samples by estimating
the equivalent number of independent draws [\[5\]](#ref5). For MCMC, ESS
is computed from the autocorrelation structure of the chain. For
variational inference, ESS comes from the importance weights when
self-normalized importance resampling (SIR) is applied to generate
corrected samples. For Laplace approximation without SIR, ESS equals the
number of generated samples (since they are independently drawn from the
Gaussian approximation). Dividing ESS by wall-clock time in seconds
produces a rate metric that is comparable across methods regardless of
their internal mechanics.

**Accuracy: Wasserstein-1 and KLD.** Given a reference distribution
(analytical truth, consensus, or a high-quality fit),
[`Arena()`](https://robustecologies.github.io/lucifer/reference/Arena.md)
measures how close each method’s marginal posteriors are to the
reference. The Wasserstein-1 distance (earth mover’s distance) between
two univariate distributions with CDFs \\F\\ and \\G\\ is

\\W_1(F, G) = \int_0^1 \|F^{-1}(u) - G^{-1}(u)\| \\ du,\\

which
[`Arena()`](https://robustecologies.github.io/lucifer/reference/Arena.md)
approximates via matched quantiles: \\\hat{W}\_1 = n^{-1}\sum\_{i=1}^n
\|x\_{(i)} - y\_{(i)}\|\\, where \\x\_{(i)}\\ and \\y\_{(i)}\\ are the
order statistics. Unlike the Kullback-Leibler divergence, \\W_1\\ is
always finite even when the supports do not overlap, making it robust
for comparing methods that may explore different regions of the
parameter space. The KL divergence \\D\_\text{KL}(p \\ q) = \int p(x)
\log(p(x)/q(x))\\dx\\ is also computed via kernel density estimates for
users who prefer this information-theoretic measure.

Both metrics are computed marginally for each parameter and averaged
across all \\d\\ parameters to produce a single accuracy score per
method.

**Reliability: convergence and ranking.** Each method is assessed for
convergence using method-specific criteria (stationarity for MCMC, ELBO
convergence for VB, optimizer convergence for Laplace). The composite
ranking combines efficiency, accuracy, and reliability by averaging the
individual ranks, producing a single overall rank per method.

### Pareto frontier

In the time-accuracy plane, methods that are neither faster nor more
accurate than some other method are Pareto-dominated. The Pareto
frontier consists of methods that are non-dominated: each represents the
best accuracy achievable at a given time budget, or the fastest method
at a given accuracy level.
[`Arena()`](https://robustecologies.github.io/lucifer/reference/Arena.md)
identifies the frontier by sorting methods by time and filtering for
monotonically improving accuracy, enabling users to see the tradeoff
curve and choose based on their computational budget.

## Using `Prescribe()`

### A simple low-dimensional model

Consider a normal-normal conjugate model with two parameters (\\\mu\\
and \\\log\sigma\\) estimated from 100 observations. This is the
simplest possible scenario: low-dimensional, fast-evaluating, unimodal,
and well-conditioned.

``` r

library(lucifer)

# Normal-normal model: two parameters
Model <- function(parm, Data) {
  mu <- parm[1]
  sigma <- exp(parm[2])
  LL <- sum(dnorm(Data$y, mu, sigma, log = TRUE))
  LP <- LL + dnormv(mu, 0, 1000, log = TRUE) +
    dhalfcauchy(sigma, 25, log = TRUE)
  yhat <- rep(mu, Data$N)
  return(list(LP = LP, Dev = -2 * LL, Monitor = LP,
    yhat = yhat, parm = parm))
}

set.seed(42)
N <- 100
y <- rnorm(N, mean = 5, sd = 2)
Data <- list(N = N, y = y, mon.names = "LP",
  parm.names = c("mu", "log.sigma"))
IV <- c(0, 0)

# Profile the model and get recommendation
rx <- Prescribe(Model, Data, IV, n.profile = 200)
```

The `print` method shows a one-screen summary:

``` r

print(rx)
```

For a 2-parameter, fast-evaluating model,
[`Prescribe()`](https://robustecologies.github.io/lucifer/reference/Prescribe.md)
should rank iterative quadrature (CAGH) and Laplace approximation highly
because they provide near-exact inference in very low dimensions. NUTS
and other gradient-based MCMC should also score well. Componentwise
methods (AMWG) will be competitive because the per-component overhead is
negligible at \\d = 2\\.

The `summary` method provides the complete scoring table:

``` r

summary(rx)
```

The `plot` method produces a two-panel figure: a horizontal bar chart of
scores on the left, and a radar chart of the model fingerprint on the
right.

``` r

plot(rx)
```

The landscape plot provides an alternative view: a bubble chart where
each algorithm is positioned by its MCMC subcategory (x-axis) and score
(y-axis), sized by evaluation cost multiplier, and colored by
subcategory. This reveals the algorithmic diversity within the MCMC
family at a glance.

``` r

plot(rx, type = "landscape", n.top = 30)
```

### Querying the algorithm registry

The
[`algo_info()`](https://robustecologies.github.io/lucifer/reference/algo_info.md)
function provides programmatic access to the algorithm registry,
returning a data frame of all methods with their metadata. It supports
filtering by category, subcategory, or gradient requirement.

``` r

# All MCMC gradient-based methods
grad_mcmc <- algo_info(category = "MCMC", subcategory = "gradient")
knitr::kable(grad_mcmc[, c("Abbreviation", "Full.Name", "Quality.Tier",
                             "Eval.Cost")],
             caption = "Gradient-based MCMC algorithms")
```

``` r

# All multimodal-specialized methods
multimodal <- algo_info(subcategory = "multimodal")
knitr::kable(multimodal[, c("Abbreviation", "Full.Name",
                              "Multimodal.Affinity", "Quality.Tier")],
             caption = "Multimodal tempering algorithms")
```

### A moderate-dimensional regression with gradients

Now consider a Bayesian logistic regression with \\d = 20\\ parameters
estimated from 500 observations, where user gradients are supplied.

``` r

# Logistic regression model
Model_logistic <- function(parm, Data) {
  beta <- parm
  eta <- Data$X %*% beta
  p <- 1 / (1 + exp(-eta))
  p <- pmax(pmin(p, 1 - 1e-8), 1e-8)
  LL <- sum(Data$y * log(p) + (1 - Data$y) * log(1 - p))
  LP <- LL + sum(dnormv(beta, 0, 100, log = TRUE))
  yhat <- rbinom(Data$N, 1, p)
  return(list(LP = LP, Dev = -2 * LL, Monitor = LP,
    yhat = yhat, parm = parm))
}

set.seed(123)
d <- 20
N <- 500
X <- cbind(1, matrix(rnorm(N * (d - 1)), ncol = d - 1))
beta_true <- rnorm(d, 0, 0.5)
p_true <- 1 / (1 + exp(-X %*% beta_true))
y <- rbinom(N, 1, p_true)

Data_logistic <- list(
  N = N, X = X, y = y,
  mon.names = "LP",
  parm.names = paste0("beta", 0:(d - 1)),
  user.grad = TRUE  # signal gradient availability
)
IV_logistic <- rep(0, d)

rx_logistic <- Prescribe(Model_logistic, Data_logistic, IV_logistic,
                          n.profile = 100)
print(rx_logistic)
```

At \\d = 20\\ with gradients available,
[`Prescribe()`](https://robustecologies.github.io/lucifer/reference/Prescribe.md)
should strongly favor NUTS, which achieves near-optimal scaling by
exploiting Hamiltonian dynamics. Pathfinder (VB) will also score well as
a fast approximate alternative. Random-walk methods should be penalized
for their \\O(d^2)\\ scaling.

### A high-dimensional hierarchical model

For a hierarchical model with \\d = 200\\ parameters, the landscape
changes dramatically. Componentwise methods become impractical (200
model evaluations per iteration), random-walk Metropolis mixes
prohibitively slowly, and only gradient-based MCMC and variational
inference remain competitive.

``` r

# Simulate scenario with d = 200
# (Model code omitted for brevity; use any hierarchical model)
# The key insight is the profile:

# rx_hier <- Prescribe(Model_hier, Data_hier, IV_hier)
# With d = 200, fast evaluation, and no gradients:
#   - VB (Pathfinder, ADVI.mf) will dominate
#   - NUTS will score well if gradients become available
#   - All componentwise and random-walk methods will score near zero
```

This illustrates the core value of
[`Prescribe()`](https://robustecologies.github.io/lucifer/reference/Prescribe.md):
it prevents the user from wasting hours running AMWG on a 200-parameter
model when variational inference would complete in seconds.

## Using `Arena()`

### Mode A: comparing pre-computed fits

The most common use case is comparing fits that have already been
computed. Suppose we fit the normal-normal model with three different
methods:

``` r

# Fit with three methods
fit_nuts <- lucifer(Model, Data, IV,
  Iterations = 10000, Status = 1000, Thinning = 5,
  Algorithm = "No-U-Turn Sampler")

fit_amwg <- lucifer(Model, Data, IV,
  Iterations = 10000, Status = 1000, Thinning = 5,
  Algorithm = "Adaptive Metropolis-within-Gibbs")

fit_vb <- VariationalBayes(Model, Data, parm = IV,
  Iterations = 1000, Method = "Pathfinder", sir = TRUE)

fit_laplace <- LaplaceApproximation(Model, Data, parm = IV,
  Iterations = 500, Method = "SPG", sir = TRUE)

# Compare all four in the Arena
arena <- Arena(list(
  NUTS    = fit_nuts,
  AMWG    = fit_amwg,
  VB      = fit_vb,
  Laplace = fit_laplace
))
```

The `print` method shows a compact metrics table:

``` r

print(arena)
```

The output includes ESS per second (efficiency), Wasserstein distance to
the consensus reference (accuracy), convergence status, and overall
rank.

The `summary` method provides the full comparison:

``` r

summary(arena)
```

This includes per-parameter posterior means and standard deviations
across all methods, the pairwise Wasserstein distance matrix (revealing
which methods agree and which diverge), and the Pareto frontier.

### Analytical reference

When analytical posterior moments are available (conjugate models,
simulation studies), pass them as the reference for accuracy
computation:

``` r

# For the normal-normal model, the analytical posterior is known
# (assuming flat prior on mu, known-sigma case for simplicity)
analytical_ref <- list(
  mean = c(mean(y), log(sd(y))),
  sd   = c(sd(y) / sqrt(N), 1 / sqrt(2 * N))
)

arena_ref <- Arena(list(
  NUTS    = fit_nuts,
  AMWG    = fit_amwg,
  VB      = fit_vb,
  Laplace = fit_laplace
), reference = analytical_ref)

print(arena_ref)
```

### Mode B: automatic fitting

For convenience,
[`Arena()`](https://robustecologies.github.io/lucifer/reference/Arena.md)
can run methods automatically given a specification list:

``` r

arena_auto <- Arena(
  Model = Model, Data = Data, Initial.Values = IV,
  methods = list(
    NUTS = list(type = "MCMC", Algorithm = "No-U-Turn Sampler",
                Iterations = 5000),
    AMWG = list(type = "MCMC", Algorithm = "Adaptive Metropolis-within-Gibbs",
                Iterations = 5000),
    RAM  = list(type = "MCMC", Algorithm = "Robust Adaptive Metropolis",
                Iterations = 5000),
    VB   = list(type = "VB", Method = "Pathfinder"),
    Laplace = list(type = "Laplace", Method = "SPG")
  )
)
```

## Visualization battery

[`Arena()`](https://robustecologies.github.io/lucifer/reference/Arena.md)
produces eight distinct plot types, each revealing a different aspect of
the comparison.

### Efficiency bar chart

The default plot shows methods ordered by ESS per second, the primary
efficiency metric. Color encodes the inference category (MCMC, VB,
Laplace, etc.), and opacity distinguishes converged from non-converged
methods.

``` r

plot(arena, type = "efficiency")
```

### Marginal posterior density overlays

The accuracy plot overlays the marginal posterior densities from each
method, with the reference shown as a thick black dashed line. This is
the most direct visualization of how different methods approximate the
same posterior.

``` r

plot(arena, type = "accuracy")
```

If the model has many parameters, use the `Parms` argument to select a
subset:

``` r

plot(arena, type = "accuracy", Parms = 1:4)
```

### Pareto frontier

The Pareto plot reveals the time-accuracy tradeoff. Methods on the
frontier (filled points, connected by dashed line) represent optimal
choices at different computational budgets. Dominated methods (open
circles) are always inferior to some frontier method.

``` r

plot(arena, type = "pareto")
```

### Convergence traces

For each method, the convergence plot shows the relevant trace:
log-posterior for MCMC and Laplace, ELBO for VB, and ESS history for
SMC. This quickly reveals whether methods have converged or are still in
a transient phase.

``` r

plot(arena, type = "convergence")
```

### Spider/radar chart

The radar chart normalizes four performance axes (efficiency, accuracy,
speed, reliability) to \\\[0,1\]\\ and draws one polygon per method.
This provides an instant gestalt of each method’s strengths and
weaknesses relative to the others.

``` r

plot(arena, type = "radar")
```

### Rank heatmap

The heatmap shows the rank of each method on each metric, with color
intensity proportional to rank (darker = better). This is the most
compact summary when comparing many methods across multiple criteria.

``` r

plot(arena, type = "heatmap")
```

### Forest plot

The forest plot shows posterior means and 95% credible intervals for
each method on selected parameters, with the reference mean as a
vertical dashed line. This is the standard format for meta-analytic
displays and is intuitive for comparing point estimates and uncertainty.

``` r

plot(arena, type = "forest", Parms = c("mu", "log.sigma"))
```

### Pairwise Wasserstein heatmap

The pairwise plot shows the Wasserstein distance between every pair of
methods, averaged across parameters. Methods that agree closely (low
distance) appear light; methods that diverge appear dark. This is
particularly useful for detecting which methods form consensus clusters
and which are outliers.

``` r

plot(arena, type = "pairwise")
```

## Recommended workflow

The intended workflow combines
[`Prescribe()`](https://robustecologies.github.io/lucifer/reference/Prescribe.md),
method fitting, the existing
[`Consort()`](https://robustecologies.github.io/lucifer/reference/Consort.md)
function (for MCMC-specific iterative tuning), and
[`Arena()`](https://robustecologies.github.io/lucifer/reference/Arena.md):

``` r

# Step 1: Profile the problem and get recommendation
rx <- Prescribe(Model, Data, IV)

# Step 2: Fit the recommended method plus alternatives
# (use the ready-to-paste code from rx$Recommend)
fit1 <- lucifer(Model, Data, IV,
  Iterations = 10000, Status = 1000, Thinning = 5,
  Algorithm = "No-U-Turn Sampler")

fit2 <- VariationalBayes(Model, Data, parm = IV,
  Iterations = 1000, Method = "Pathfinder", sir = TRUE)

fit3 <- LaplaceApproximation(Model, Data, parm = IV,
  Iterations = 500, Method = "SPG", sir = TRUE)

# Step 3: For MCMC fits, use Consort for iterative diagnostics
# Consort(fit1)

# Step 4: Compare all methods in the Arena
arena <- Arena(list(NUTS = fit1, Pathfinder = fit2, Laplace = fit3))
print(arena)
plot(arena)
plot(arena, type = "accuracy")
```

## Advanced example: comparing across all inference families

For a problem with both a likelihood-based formulation and a
simulation-based formulation,
[`Arena()`](https://robustecologies.github.io/lucifer/reference/Arena.md)
can compare methods from every available family:

``` r

# Suppose we have a model that supports both likelihood-based
# and simulation-based inference

# MCMC
fit_nuts <- lucifer(Model, Data, IV,
  Iterations = 20000, Status = 2000, Thinning = 10,
  Algorithm = "No-U-Turn Sampler")

# Adaptive random-walk
fit_ram <- lucifer(Model, Data, IV,
  Iterations = 20000, Status = 2000, Thinning = 10,
  Algorithm = "Robust Adaptive Metropolis")

# Ensemble MCMC
fit_aies <- lucifer(Model, Data, IV,
  Iterations = 20000, Status = 2000, Thinning = 10,
  Algorithm = "Affine-Invariant Ensemble Sampler")

# Variational Bayes
fit_pathfinder <- VariationalBayes(Model, Data, parm = IV,
  Iterations = 1000, Method = "Pathfinder", sir = TRUE)

fit_advi <- VariationalBayes(Model, Data, parm = IV,
  Iterations = 1000, Method = "ADVI.mf", sir = TRUE)

# Laplace approximation
fit_laplace <- LaplaceApproximation(Model, Data, parm = IV,
  Iterations = 500, Method = "SPG", sir = TRUE)

# Iterative quadrature (only if d <= 10)
fit_iq <- IterativeQuadrature(Model, Data, parm = IV,
  Iterations = 100, Algorithm = "CAGH", sir = TRUE)

# PMC
fit_pmc <- PMC(Model, Data, IV,
  Iterations = 20, Thinning = 10, N = 1000)

# SMC
fit_smc <- SMC(Model, Data, IV,
  N.particles = 2000, Rejuvenation.steps = 5)

# Full Arena comparison
arena_full <- Arena(list(
  NUTS       = fit_nuts,
  RAM        = fit_ram,
  AIES       = fit_aies,
  Pathfinder = fit_pathfinder,
  ADVI       = fit_advi,
  Laplace    = fit_laplace,
  IQ         = fit_iq,
  PMC        = fit_pmc,
  SMC        = fit_smc
))

# Summary with all pairwise comparisons
summary(arena_full)

# Visualize the full comparison
plot(arena_full, type = "efficiency")
plot(arena_full, type = "pareto")
plot(arena_full, type = "radar")
plot(arena_full, type = "heatmap")
plot(arena_full, type = "accuracy", Parms = 1:4)
```

## Implementation details

### The algorithm registry

Both
[`Prescribe()`](https://robustecologies.github.io/lucifer/reference/Prescribe.md)
and
[`Arena()`](https://robustecologies.github.io/lucifer/reference/Arena.md)
rely on a shared algorithm property registry (`.algo_registry`) that
catalogs every inference method in the package. Each entry stores the
abbreviation, inference category, subcategory, and properties including
gradient requirement, discrete parameter support, componentwise
evaluation, quality tier (1 = state-of-the-art, 2 = solid, 3 = niche),
multimodal affinity, constraint affinity, torch requirement, evaluation
cost multiplier, and recommended dimension range. The MCMC category
alone spans 13 subcategories: gradient (NUTS, GIST, autoMALA, MALT,
Barker, and others), adaptive random-walk (RAM, DRAM, AM), componentwise
(AMWG, Gibbs, Zanella), slice (AFSS, ESS, QSS), ensemble (AIES, DREAM,
Zeus), QMC (8 quasi-Monte Carlo methods), geometric (RMHMC, LMC, MHMC,
exploiting Riemannian metric tensors for ill-conditioned posteriors),
PDMP (BPS, ZigZag, Boomerang, continuous-time piecewise deterministic
samplers), multimodal (SimTemp, NRST, NRPT, WL, tempering and
flat-histogram methods), constraint (ProjLang, ProxMCMC, native
constraint handling via projection), flow (NeuTra, flowMC, neural
transport requiring torch), surrogate (DA, delayed acceptance with cheap
surrogates for expensive likelihoods), and augmentation (PG,
Polya-Gamma). The registry provides a single source of truth for
algorithm metadata used by
[`Consort()`](https://robustecologies.github.io/lucifer/reference/Consort.md),
[`Prescribe()`](https://robustecologies.github.io/lucifer/reference/Prescribe.md),
[`Arena()`](https://robustecologies.github.io/lucifer/reference/Arena.md),
and
[`Crucible()`](https://robustecologies.github.io/lucifer/reference/Crucible.md).
The exported function
[`algo_info()`](https://robustecologies.github.io/lucifer/reference/algo_info.md)
enables users to query this registry, filtering by category,
subcategory, or gradient requirement. Helper functions `.algo_abbrev()`
and `.algo_lookup()` enable bidirectional lookup between full algorithm
names and abbreviations.

### Extraction adapters

[`Arena()`](https://robustecologies.github.io/lucifer/reference/Arena.md)
uses class-specific adapter functions for each S3 class (demonoid, vb,
laplace, iterquad, pmc, smc, abc, sbi) to extract a standardized list
from each fit object. Specialized engines that inherit from demonoid
(sde_fit) are handled automatically through the demonoid adapter with
appropriate category labeling. The standardized list contains posterior
samples, parameter means and standard deviations, ESS, wall-clock time,
convergence status, log marginal likelihood (when available), DIC (when
available), and acceptance rate (when applicable). When posterior
samples are not directly available (e.g., Laplace without SIR), the
adapter generates samples from the Gaussian approximation
\\\mathcal{N}(\hat{\theta}, \hat{\Sigma})\\.

### Consensus reference

When no analytical reference is provided,
[`Arena()`](https://robustecologies.github.io/lucifer/reference/Arena.md)
constructs a consensus reference by computing an ESS-weighted average of
the posterior means and variances from all methods:

\\\mu\_\text{ref} = \sum\_{m=1}^M w_m \hat{\mu}\_m, \qquad
\sigma^2\_\text{ref} = \sum\_{m=1}^M w_m \left(\hat{\sigma}\_m^2 +
(\hat{\mu}\_m - \mu\_\text{ref})^2\right),\\

where \\w_m = \text{ESS}\_m / \sum_j \text{ESS}\_j\\. This weights more
reliable methods more heavily and incorporates between-method variance
in the reference uncertainty. Samples from
\\\mathcal{N}(\mu\_\text{ref}, \sigma^2\_\text{ref})\\ are then used as
the reference distribution for computing divergence metrics.

## References

**\[1\]** Rice, J. R. (1976). “The Algorithm Selection Problem”.
*Advances in Computers*, 15, 65-118.
<doi:%5B10.1016/S0065-2458(08)60520-3>\](<https://doi.org/10.1016/S0065-2458(08)60520-3>)

**\[2\]** Wolpert, D. H. and Macready, W. G. (1997). “No Free Lunch
Theorems for Optimization”. *IEEE Transactions on Evolutionary
Computation*, 1(1), 67-82.
<doi:%5B10.1109/4235.585893>\](<https://doi.org/10.1109/4235.585893>)

**\[3\]** Betancourt, M. (2017). “A Conceptual Introduction to
Hamiltonian Monte Carlo”. *arXiv preprint arXiv:1701.02434*.
<doi:%5B10.48550/arXiv.1701.02434>\](<https://doi.org/10.48550/arXiv.1701.02434>)

**\[4\]** Roberts, G. O. and Rosenthal, J. S. (2001). “Optimal Scaling
for Various Metropolis-Hastings Algorithms”. *Statistical Science*,
16(4), 351-367.
<doi:%5B10.1214/ss/1015346320>\](<https://doi.org/10.1214/ss/1015346320>)

**\[5\]** Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and
Burkner, P.-C. (2021). “Rank-Normalization, Folding, and Localization:
An Improved Rhat for Assessing Convergence of MCMC”. *Bayesian
Analysis*, 16(2), 667-718.
<doi:%5B10.1214/20-BA1221>\](<https://doi.org/10.1214/20-BA1221>)

**\[6\]** Blei, D. M., Kucukelbir, A. and McAuliffe, J. D. (2017).
“Variational Inference: A Review for Statisticians”. *Journal of the
American Statistical Association*, 112(518), 859-877.
<doi:%5B10.1080/01621459.2017.1285773>\](<https://doi.org/10.1080/01621459.2017.1285773>)

**\[7\]** Zhang, L., Carpenter, B., Gelman, A. and Vehtari, A. (2022).
“Pathfinder: Parallel Quasi-Newton Variational Inference”. *Journal of
Machine Learning Research*, 23(306), 1-49.

**\[8\]** Del Moral, P., Doucet, A. and Jasra, A. (2006). “Sequential
Monte Carlo Samplers”. *Journal of the Royal Statistical Society: Series
B*, 68(3), 411-436.
<doi:%5B10.1111/j.1467-9868.2006.00553.x>\](<https://doi.org/10.1111/j.1467-9868.2006.00553.x>)

**\[9\]** Cranmer, K., Brehmer, J. and Louppe, G. (2020). “The Frontier
of Simulation-Based Inference”. *Proceedings of the National Academy of
Sciences*, 117(48), 30055-30062.
<doi:%5B10.1073/pnas.1912789117>\](<https://doi.org/10.1073/pnas.1912789117>)

**\[10\]** Hoffman, M. D. and Gelman, A. (2014). “The No-U-Turn Sampler:
Adaptively Setting Path Lengths in Hamiltonian Monte Carlo”. *Journal of
Machine Learning Research*, 15(47), 1593-1623.

**\[11\]** Goodman, J. and Weare, J. (2010). “Ensemble Samplers with
Affine Invariance”. *Communications in Applied Mathematics and
Computational Science*, 5(1), 65-80.
<doi:%5B10.2140/camcos.2010.5.65>\](<https://doi.org/10.2140/camcos.2010.5.65>)

**\[12\]** Villani, C. (2009). *Optimal Transport: Old and New*.
Springer. ISBN: 978-3-540-71049-3.
<doi:%5B10.1007/978-3-540-71050-9>\](<https://doi.org/10.1007/978-3-540-71050-9>)
