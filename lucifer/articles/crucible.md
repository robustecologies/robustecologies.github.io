# Crucible: automated Bayesian inference

## Motivation

Bayesian inference demands more from the practitioner than specifying a
model and pressing “run.” The analyst must choose an inference engine
from a growing menu of methods, each with distinct failure modes, tuning
parameters, and convergence criteria, then diagnose whether the result
is trustworthy, and, if not, iterate. This workflow is particularly
demanding in the lucifer ecosystem, which offers over 130 inference
methods across eight families: 82 MCMC algorithms spanning 13
subcategories, 8 variational Bayes methods, 18 Laplace optimizers, 3
iterative quadrature schemes, population Monte Carlo, sequential Monte
Carlo, 4 ABC methods, and 6 SBI methods. The combinatorial space of
method-plus-tuning configurations is vast, and navigating it manually
requires both statistical expertise and substantial analyst time.

The standard manual workflow proceeds in four stages. First, the analyst
calls
[`Prescribe()`](https://robustecologies.github.io/lucifer/reference/Prescribe.md)
to profile the model and obtain a ranked recommendation of inference
strategies [\[1\]](#ref1). Second, she fits one or more recommended
methods using
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`VariationalBayes()`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md),
[`LaplaceApproximation()`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
or another fitting function. Third, she inspects convergence with
[`Consort()`](https://robustecologies.github.io/lucifer/reference/Consort.md)
and decides whether the fit is adequate. Fourth, if multiple methods
have been tried, she compares them with
[`Arena()`](https://robustecologies.github.io/lucifer/reference/Arena.md)
to select the most efficient and accurate. Each stage requires manual
intervention, interpretation of console output, and decision-making
about how to proceed when convergence is unsatisfactory.

[`Crucible()`](https://robustecologies.github.io/lucifer/reference/Crucible.md)
automates this entire pipeline. Given a model, data, and initial values,
it profiles the problem, selects the most promising methods, fits them
with parallel chains, diagnoses convergence via a modernized
[`Consort()`](https://robustecologies.github.io/lucifer/reference/Consort.md),
iteratively refines non-converged fits using structured suggestions,
compares all successful fits in
[`Arena()`](https://robustecologies.github.io/lucifer/reference/Arena.md),
and returns a ranked summary with the best method identified. The name
evokes a severe test in which impurities are burned away, fitting the
package’s theological metaphor alongside lucifer (the light-bearer),
Consort (the advisor), and Arena (the arena of combat).

The automation rests on two key innovations. First,
[`Consort()`](https://robustecologies.github.io/lucifer/reference/Consort.md)
has been rewritten from a 1158-line monolith that only handled MCMC
objects and printed to the console into a modular diagnostic engine that
supports all six inference families and returns a structured S3 object
with programmatic suggestions. Second, Crucible’s refinement loop uses
these suggestions to adjust fitting specifications automatically,
implementing the kind of iterative tuning that analysts normally perform
by hand.

## Convergence diagnostics: theoretical background

The reliability of any Bayesian computation depends on whether the
inference algorithm has converged to the target posterior distribution.
The diagnostic criteria differ fundamentally across inference families,
but the underlying principle is the same: we need evidence that the
approximation is close enough to the true posterior to support the
scientific conclusions drawn from it.

### MCMC diagnostics

For MCMC methods, convergence diagnostics assess whether the Markov
chain has mixed well enough to produce reliable estimates. The modern
standard, established by Vehtari et al. [\[2\]](#ref2), combines three
complementary measures.

The split \\\hat{R}\\ statistic compares within-chain and between-chain
variances after splitting each chain in half, folding the draws to
detect problems in the tails, and rank-normalizing to remove sensitivity
to heavy tails. The recommended threshold is \\\hat{R} \< 1.01\\,
stricter than the traditional \\\hat{R} \< 1.1\\ from Gelman and Rubin
[\[3\]](#ref3), which has been shown to be insufficiently conservative
for many practical problems [\[2\]](#ref2). A value near 1.0 indicates
that the chains have mixed well and are sampling from the same
distribution; values above 1.01 suggest that the chains have not yet
converged or are exploring different regions of the parameter space.

The effective sample size (ESS) quantifies how many independent draws
the correlated MCMC output is worth. Bulk ESS measures the reliability
of central tendency estimates (means, medians), while tail ESS measures
the reliability of quantile estimates (credible intervals, tail
probabilities) [\[2\]](#ref2). The recommended thresholds are ESS bulk
\\\geq 400\\ and ESS tail \\\geq 200\\, which ensure that Monte Carlo
standard errors for posterior means and 95% intervals are below 5% of
the posterior standard deviation. These thresholds are substantially
higher than the traditional ESS \\\geq 100\\, reflecting the insight
that reliable inference, particularly for interval estimation, requires
more effective samples than previously appreciated.

The Monte Carlo standard error (MCSE) provides a direct measure of
estimation precision. For a parameter with posterior standard deviation
\\\sigma\\, an MCSE of \\\epsilon\\ means that the posterior mean
estimate has a standard error of \\\epsilon\\. Flegal, Haran, and Jones
[\[4\]](#ref4) demonstrated that MCMC practitioners routinely
underestimate the uncertainty in their estimates and recommended
reporting MCSE alongside posterior summaries.
[`Consort()`](https://robustecologies.github.io/lucifer/reference/Consort.md)
checks that MCSE \\\< 6.27\\\\ of \\\sigma\\ for all parameters, a
threshold derived from requiring that the Monte Carlo error contributes
less than 0.2% additional variance to the posterior variance.

Divergent transitions are a diagnostic specific to Hamiltonian Monte
Carlo and its variants (HMC, NUTS, HMCDA, MALA). A divergence occurs
when the numerical integration of Hamilton’s equations encounters a
region where the step size is too large to accurately track the
posterior geometry, typically near regions of high curvature such as
funnel-shaped posteriors [\[5\]](#ref5). Even a small number of
divergences can bias posterior estimates because they indicate
systematic failure to explore part of the posterior.
[`Consort()`](https://robustecologies.github.io/lucifer/reference/Consort.md)
reports the divergence count and flags any non-zero count as a
diagnostic failure.

For gradient-based algorithms,
[`Consort()`](https://robustecologies.github.io/lucifer/reference/Consort.md)
also checks the energy Bayesian fraction of missing information
(E-BFMI), computed as the variance of the energy transition divided by
the marginal energy variance. A value below 0.3 indicates that the
sampler is failing to explore the energy distribution adequately, often
due to regions where the posterior geometry changes abruptly. This check
applies to all gradient-requiring MCMC methods, not only NUTS.

Beyond the binary pass/fail for each condition,
[`Consort()`](https://robustecologies.github.io/lucifer/reference/Consort.md)
computes a continuous convergence score between 0 and 100 by combining
all diagnostic metrics through saturating transforms. The Rhat component
uses a logistic decay centered at the threshold; the ESS components
saturate at twice the threshold; the MCSE component linearly decays; and
the divergence and E-BFMI components provide sharp penalties. The
composite score does not replace the pass/fail assessment for
determining appeasement, but it enables Arena to rank methods on a finer
scale than binary convergence status alone. The score and the ESS per
second metric are displayed by `print(consort_object)` and can be
visualized with `plot(consort_object, type = "scorecard")`.

All diagnostic thresholds are configurable via the `thresholds` argument
to
[`Consort()`](https://robustecologies.github.io/lucifer/reference/Consort.md),
which accepts a list produced by `.consort_thresholds()`. The defaults
implement the modern recommendations discussed above but can be
overridden for specialized applications.

### Variational inference diagnostics

Variational inference optimizes a tractable approximation \\q(\theta)\\
to the true posterior \\p(\theta \mid y)\\ by maximizing the evidence
lower bound (ELBO),

\\\text{ELBO}(q) = \mathbb{E}\_{q}\left\[\log p(y, \theta) - \log
q(\theta)\right\] \leq \log p(y),\\

which is equivalent to minimizing the Kullback-Leibler divergence
\\\text{KL}(q \\ p)\\ [\[6\]](#ref6). Because variational inference
solves an optimization problem rather than a sampling problem, the
relevant diagnostics differ from MCMC.
[`Consort()`](https://robustecologies.github.io/lucifer/reference/Consort.md)
checks up to four conditions for VB fits: whether the optimizer has
converged (the algorithm’s internal convergence flag), whether the ELBO
trace is stable in its final segment (coefficient of variation below 5%
in the last 10% of iterations, indicating that the optimization has
settled rather than oscillating), whether the Pareto \\k\\ diagnostic
from self-importance sampling is below 0.7 [\[7\]](#ref7), and, when
gradient information is stored, whether the ELBO gradient norm at
convergence is below 1.0 (indicating a true local optimum rather than
premature termination). A Pareto \\k \> 0.7\\ indicates that the
variational approximation is too far from the true posterior for
importance sampling corrections to be reliable, suggesting that the
model may need a richer variational family or a switch to MCMC.

### Laplace approximation diagnostics

The Laplace approximation finds the posterior mode \\\hat{\theta}\\ and
approximates the posterior as Gaussian with covariance equal to the
inverse observed information matrix: \\p(\theta \mid y) \approx
\mathcal{N}(\hat{\theta}, \mathcal{I}(\hat{\theta})^{-1})\\. This is
adequate when the posterior is approximately Gaussian near its mode,
which is guaranteed asymptotically for regular models but may fail
catastrophically for finite samples, multimodal posteriors, or
posteriors with heavy tails or constraints [\[8\]](#ref8).
[`Consort()`](https://robustecologies.github.io/lucifer/reference/Consort.md)
checks four conditions: convergence of the optimizer (the internal
convergence flag), positive-definiteness of the Hessian at the mode
(which ensures the covariance matrix is valid; a non-positive-definite
Hessian indicates a saddle point or ridge), the gradient norm at the
mode being near zero (below \\10^{-3}\\, confirming that the optimizer
has actually found a stationary point rather than stopping prematurely),
and the Hessian condition number being below \\10^6\\ (an
ill-conditioned Hessian produces unreliable standard errors even when
positive-definite).

### Iterative quadrature diagnostics

Iterative quadrature methods approximate the posterior by numerical
integration on an adaptive grid. Their diagnostics focus on convergence
of the iterative scheme and stability of the parameter estimates across
iterations.
[`Consort()`](https://robustecologies.github.io/lucifer/reference/Consort.md)
checks convergence and parameter stability for these methods.

## The Crucible pipeline

Crucible orchestrates the inference workflow through six stages, each
corresponding to a distinct phase of the analyst’s decision process. The
stages run sequentially, with methods fitted one at a time to avoid
nested parallelism (each MCMC method already uses parallel chains
internally).

### Stage 1: profiling

The pipeline begins by calling
[`Prescribe()`](https://robustecologies.github.io/lucifer/reference/Prescribe.md)
to profile the model. This evaluates the Model function approximately
\\n\_\text{profile} + 2d\\ times, measuring evaluation speed, checking
gradient availability, computing the Hessian condition number (for
models with \\d \leq 100\\), estimating a multimodality score, and
classifying the problem along six axes. The cost is typically
negligible: a few seconds for models with fast evaluation, a fraction of
a minute even for slow models. The output is a ranked list of all
inference methods with composite scores reflecting their expected
suitability for the problem.

### Stage 2: fitting

Crucible selects the top `n_methods` methods from Prescribe’s ranking,
filtered by the `families` argument. For MCMC methods, it constructs
sensible defaults: iteration counts scaled by dimensionality (10000 for
\\d \leq 20\\, 20000 for \\d \leq 100\\, 50000 for \\d \> 100\\),
thinning to retain approximately 2000 samples, and algorithm-specific
Specs generated from an internal registry. For mandatory-Specs
algorithms such as AIES and DEMC (ensemble methods that require
specifying the number of walkers), Crucible provides appropriate
defaults (\\N_c = 2d\\ for AIES, \\N_c = 3d\\ for DEMC, following the
recommendations of Goodman and Weare [\[9\]](#ref9) and ter Braak
[\[10\]](#ref10)). Each MCMC method is run with `Chains` parallel chains
using the package’s parallel infrastructure. Non-MCMC methods (VB,
Laplace, IQ) are called with their standard interfaces and with
`sir = TRUE` to enable self-importance resampling for posterior sample
generation.

### Stage 3: diagnosis

Each successful fit is passed to
[`Consort()`](https://robustecologies.github.io/lucifer/reference/Consort.md),
which dispatches to the appropriate family-specific diagnostic engine.
The result is a `consort` S3 object containing the family label, a
logical `appeased` flag, a conditions table (data.frame with condition
name, pass/fail status, and detail text), family-specific diagnostics
(per-parameter MCSE, ESS, and \\\hat{R}\\ for MCMC; ELBO stability
metrics for VB; Hessian condition for Laplace), and a structured
suggestion if the fit is not appeased.

### Stage 4: refinement

For each non-appeased fit, Crucible extracts the structured suggestion
from Consort and re-fits with adjusted specifications. The suggestion
engine is registry-driven rather than hardcoded: it looks up the
algorithm’s properties in the internal registry, assesses which
diagnostic conditions failed, and applies a prioritized decision tree.
The priorities, in decreasing order, are: fix divergent transitions
(reduce step size or switch to a non-gradient method), handle adaptation
failures (switch to a non-adaptive algorithm or a universal
componentwise fallback like AFSS), run the non-adaptive phase (set
`Adaptive = 0`), adjust acceptance rate (modify tuning parameters or
escalate to a more capable algorithm), and increase iterations (multiply
by \\\max(2, \lceil\text{target ESS} / \text{current ESS}\rceil)\\).

MCMC warm-starting is handled via
[`as.initial.values()`](https://robustecologies.github.io/lucifer/reference/as.initial.values.md),
which extracts the last posterior draw as the starting point for the
next round. This avoids re-running the adaptive phase from scratch and
allows the chain to continue exploring from a region of high posterior
density.

The refinement loop repeats for up to `max_rounds` total rounds (initial
fit plus `max_rounds - 1` refinements). Setting `max_rounds = 1`
disables refinement entirely.

### Stage 5: Arena comparison

Once all methods have been fitted and refined, Crucible collects the
final fits (the last round’s fit for each method) and passes them to
[`Arena()`](https://robustecologies.github.io/lucifer/reference/Arena.md)
for cross-method comparison. Arena computes unified metrics across all
inference families: effective samples per second (efficiency),
Wasserstein distance to a reference distribution (accuracy), and
convergence reliability. Methods are ranked by a composite score that
balances efficiency and accuracy [\[1\]](#ref1).

If only one method produced a successful fit, Arena is skipped and the
ranking defaults to appeasement status and ESS.

### Stage 6: assembly

The final stage assembles all results into a `crucible` S3 object
containing the Prescribe prescription, per-method results (fit, Consort
diagnostics, and refinement history for each round), the Arena
comparison, a summary table ranking all methods, and a `Best` field
identifying the winning method with its fit and Consort object ready for
immediate use.

## Example 1: normal-normal model

The simplest possible Bayesian model serves as a clear demonstration of
the pipeline. We observe \\y_1, \ldots, y_N \sim \mathcal{N}(\mu,
\sigma^2)\\ with weakly informative priors \\\mu \sim \mathcal{N}(0,
100^2)\\ and \\\sigma \sim \text{Half-Cauchy}(0, 25)\\. Working on the
log-scale for \\\sigma\\ ensures unconstrained parameterization: we
define \\\theta_2 = \log \sigma\\ and transform back via \\\sigma =
\exp(\theta_2)\\.

This model has a well-conditioned, unimodal posterior with known
analytical properties. For \\N = 50\\ observations drawn from
\\\mathcal{N}(5, 4)\\, the posterior for \\\mu\\ is approximately
\\\mathcal{N}(\bar{y}, \sigma^2/N)\\ and for \\\log\sigma\\ is
approximately normal near \\\log(s)\\ where \\s\\ is the sample standard
deviation. Any competent inference method should converge quickly,
making this an ideal test that the pipeline machinery works correctly.

``` r

# Normal-normal model
Model <- function(parm, Data) {
    mu <- parm[1]
    sigma <- exp(parm[2])
    LL <- sum(dnorm(Data$y, mu, sigma, log = TRUE))
    LP <- LL + dnorm(mu, 0, 100, log = TRUE) +
        dhalfcauchy(sigma, 25, log = TRUE)
    yhat <- rep(mu, Data$N)
    return(list(LP = LP, Dev = -2 * LL, Monitor = LP,
                yhat = yhat, parm = parm))
}

set.seed(42)
Data <- list(N = 50, y = rnorm(50, 5, 2),
             mon.names = "LP",
             parm.names = c("mu", "log.sigma"))
IV <- c(5, log(2))

# Run Crucible with 3 methods, 2 refinement rounds
cr <- Crucible(Model, Data, IV,
               n_methods = 3,
               max_rounds = 2,
               Chains = 2, CPUs = 1,
               families = c("MCMC"))

# Compact dashboard
print(cr)
```

The print method displays a compact dashboard with one row per attempted
method, showing the number of refinement rounds used, whether Consort
was appeased, the minimum bulk ESS across parameters, the elapsed time,
and the Arena rank. A typical output for this model looks like:

    Crucible: automated Bayesian inference
    ──────────────────────────────────────────
    Model: 2 parameters | 50 observations | Profiled in 0.0s

      Method  Category  Rounds Appeased  ESS.min  Minutes Rank
      ───────────────────────────────────────────────────────────
      AMWG    MCMC         1/2      ✔      2132     0.01    1
      DEMC    MCMC         2/2      ✘         2     0.29    2
      AIES    MCMC         2/2      ✘         2     0.79    3

    ✔ Best overall: AMWG (MCMC)
    Total time: 68.5 seconds

The extended summary adds model profile information, per-method
refinement history showing how many Consort conditions passed at each
round, Arena metrics, and the Pareto frontier.

``` r

summary(cr)
```

### Extracting the best fit

The `Best` field provides direct access to the winning method’s fit
object and its Consort diagnostics. From here the analyst can proceed to
posterior inference, prediction, or model comparison exactly as if she
had fitted the model manually.

``` r

# The winning fit
best_fit <- cr$Best$fit
best_consort <- cr$Best$consort

# Posterior samples (combined thinned draws from all chains)
posterior <- best_fit$Posterior1
cat("Posterior dimensions:", nrow(posterior), "draws x",
    ncol(posterior), "parameters\n")

# Point estimates from posterior summary
best_fit$Summary2[, c("Mean", "SD", "ESS")]
```

### Visualization

Crucible provides six plot types accessible via the `type` argument to
[`plot()`](https://rdrr.io/r/graphics/plot.default.html). The default
`"dashboard"` produces a 2x2 grid combining the Prescribe bar chart
(algorithm scores), Arena efficiency bars, ESS convergence trajectory,
and the Pareto frontier. Individual panels can be accessed via
`"prescribe"`, `"arena"`, `"convergence"`, `"accuracy"`, and
`"diagnostics"`.

``` r

# Full dashboard (2x2 grid)
plot(cr)
```

The convergence plot tracks the minimum ESS across parameters at each
refinement round for every method, with a horizontal line at the ESS =
400 threshold. Methods that cross the threshold have converged by
Consort’s standards; those that remain below have not.

``` r

# ESS trajectory across refinement rounds
plot(cr, type = "convergence")
```

The Arena plot compares methods on efficiency (effective samples per
second) and accuracy (Wasserstein distance to a reference distribution).

``` r

# Arena comparison
plot(cr, type = "arena")
```

## Working with Consort directly

While Crucible automates the full pipeline, analysts working
interactively may prefer to call
[`Consort()`](https://robustecologies.github.io/lucifer/reference/Consort.md)
directly on individual fits. The modernized Consort returns a proper S3
object rather than printing to the console, enabling programmatic
inspection and downstream use.

### Fitting and diagnosing MCMC

``` r

# Fit NUTS (gradient-based MCMC)
fit_nuts <- lucifer(Model, Data, IV,
                    Iterations = 10000, Status = 2000, Thinning = 10,
                    Algorithm = "NUTS",
                    Specs = list(A = 500, delta = 0.65,
                                 epsilon = 0.01, Lmax = 20),
                    Chains = 2, CPUs = 1)

# Diagnose
cx <- Consort(fit_nuts)
```

The print method displays a compact verdict with the algorithm name,
timing, and a condition table using colored pass/fail indicators:

    Consort with Lucifer
    ────────────────────
    Algorithm: NUTS (MCMC, gradient)
    Time: 0.45 min | Iterations: 10000 | Parameters: 2

      Condition              Status
      ─────────────────────────────
      Non-adaptive           ✔ pass
      Acceptance rate        ✔ 0.82 in [0.60, 0.90]
      MCSE                   ✔ all < 6.27% SD
      ESS (bulk/tail)        ✔ min bulk = 512 >= 400
      Rhat                   ✔ max = 1.002 < 1.01
      Divergences            ✔ 0
      Diminishing adaptation ✔ pass

    ✔ Lucifer has been appeased.

### Programmatic access

The `consort` object exposes all diagnostic information as structured
data, enabling programmatic decision-making in scripts and pipelines.

``` r

# Overall verdict
cx$appeased        # TRUE or FALSE
cx$family          # "MCMC", "VB", "Laplace", etc.
cx$algorithm       # Full algorithm name

# Condition table
cx$conditions      # data.frame: condition, status, detail

# Per-parameter diagnostics
names(cx$diagnostics)  # list with per-parameter ESS, Rhat, MCSE

# Suggestion (populated when not appeased, NULL otherwise)
if (!cx$appeased && !is.null(cx$suggestion)) {
    cx$suggestion$action     # what to do
    cx$suggestion$rationale  # why this suggestion
}
```

### Summary and per-parameter diagnostics

The summary method extends the print output with a per-parameter
diagnostic table showing MCSE ratio, bulk ESS, tail ESS, and \\\hat{R}\\
for each parameter, plus the full suggested code when the fit is not
appeased.

``` r

summary(cx)
```

### Consort plot types

``` r

# ESS and Rhat dot plots per parameter (default)
plot(cx, type = "diagnostics")

# Pass/fail bar chart of conditions
plot(cx, type = "conditions")

# Deviance trace
plot(cx, type = "trace")
```

The diagnostics plot produces a faceted display with ESS (bulk and tail)
on one panel and \\\hat{R}\\ on another, with vertical threshold lines
at ESS = 400 and \\\hat{R} = 1.01\\. Parameters falling outside the
thresholds are highlighted, making it immediately clear which parameters
are driving convergence failures.

### Diagnosing non-MCMC fits

[`Consort()`](https://robustecologies.github.io/lucifer/reference/Consort.md)
now supports all inference families. For Laplace approximation, the
condition table reports convergence, Hessian positive-definiteness, and
gradient norm at the mode. For variational Bayes, it checks convergence,
ELBO stability, and (where applicable) the Pareto \\k\\ diagnostic. The
suggestion engine provides family-appropriate recommendations: switching
to a different optimizer for Laplace, increasing iterations for VB, or
escalating to MCMC when the approximation method appears inadequate. The
following example demonstrates Consort on a Laplace fit.

``` r

# Fit Laplace approximation
fit_la <- LaplaceApproximation(Model, Data, parm = IV,
                               Iterations = 100,
                               Method = "SPG",
                               sir = TRUE)

# Diagnose
cx_la <- Consort(fit_la)
print(cx_la)

# Check appeasement (works for all fit classes)
is.appeased(fit_la)
```

## Example 2: Neal’s funnel

Neal’s funnel [\[5\]](#ref5) is a canonical stress test for MCMC
algorithms. The model defines a hierarchy where \\v \sim \mathcal{N}(0,
9)\\ controls the scale of \\d - 1\\ lower-level parameters \\\theta_j
\sim \mathcal{N}(0, \exp(v))\\. The posterior has a funnel shape in the
\\(v, \theta_j)\\ plane: when \\v\\ is large, the \\\theta_j\\ are
spread over a wide region; when \\v\\ is small, they are tightly
concentrated near zero. The dramatic change in scale across the funnel
makes this posterior notoriously difficult for fixed-step-size methods,
which either take steps too large in the narrow end (producing divergent
transitions) or too small in the wide end (producing near-zero
acceptance rates).

This model is a reduced version of the hierarchical structures that
arise naturally in Bayesian random-effects models, and its difficulty is
representative of the challenges faced by practitioners fitting
multilevel models. Betancourt and Girolami [\[11\]](#ref11) used it to
demonstrate the superiority of Hamiltonian Monte Carlo over Gibbs
sampling for hierarchical models; Betancourt [\[5\]](#ref5) used it to
explain the diagnostic significance of divergent transitions.

``` r

# Neal's funnel: d-dimensional version
d_funnel <- 3
FunnelModel <- function(parm, Data) {
    v <- parm[1]
    theta <- parm[2:Data$d]

    # Log-likelihood of hierarchical structure
    LL_v <- dnorm(v, 0, 3, log = TRUE)
    LL_theta <- sum(dnorm(theta, 0, exp(v / 2), log = TRUE))
    LL <- LL_v + LL_theta

    # Flat prior (the model IS the posterior)
    LP <- LL
    yhat <- rep(0, Data$N)
    return(list(LP = LP, Dev = -2 * LL, Monitor = LP,
                yhat = yhat, parm = parm))
}

FunnelData <- list(
    N = 1,
    d = d_funnel,
    y = 0,
    mon.names = "LP",
    parm.names = c("v", paste0("theta.", 1:(d_funnel - 1)))
)
FunnelIV <- rep(0, d_funnel)

# Run Crucible: compare gradient-based vs componentwise vs slice.
cr_funnel <- Crucible(FunnelModel, FunnelData, FunnelIV,
                      methods = c("NUTS", "AMWG", "Slice"),
                      max_rounds = 2,
                      Chains = 2, CPUs = 1)

print(cr_funnel)
summary(cr_funnel)
```

On this problem, Crucible reveals how diagnostic-driven refinement works
in practice. The Slice sampler converges at the first round because its
stepping-out procedure adapts naturally to the scale changes across the
funnel. Both NUTS and AMWG fail the diminishing adaptation diagnostic,
triggering Consort’s escalation to AFSS (the universal full-conditional
slice sampler) in round 2. This automatic algorithm switching is one of
Crucible’s key features: when a method struggles, the suggestion engine
identifies a more robust alternative rather than simply increasing
iterations. The escalation map covers all 82 MCMC algorithms with
sensible within-family fallback chains: geometric methods (RMHMC, LMC,
MHMC) escalate to NUTS, PDMP samplers (BPS, ZigZag) fall back through
BPS to NUTS, multimodal methods (SimTemp, NRST, WL, NRPT) escalate
within the tempering family before falling to AFSS, and flow-enhanced
methods (NeuTra, flowMC) revert to NUTS when the neural transport fails.
The convergence plot shows the contrast between Slice’s immediate
success and the other methods’ continued difficulty even after switching
to AFSS. When analytical gradients are available, NUTS would typically
handle funnel geometry more effectively by following the posterior
curvature.

``` r

# Inspect the best method's diagnostics
plot(cr_funnel, type = "diagnostics")

# Compare ESS trajectories across methods
plot(cr_funnel, type = "convergence")
```

The convergence plot for the funnel model is instructive: it shows how
the Slice sampler crosses the ESS \\\geq 400\\ threshold at round 1,
while NUTS and AMWG remain below the threshold even after refinement to
AFSS. The funnel geometry, where the posterior scale varies by orders of
magnitude depending on the value of \\v\\, is fundamentally challenging
for methods that rely on a fixed proposal scale or step size.

The refinement trajectory heatmap reveals which diagnostic conditions
drive each method’s failure and how they evolve across rounds. Each
facet corresponds to a method, each column to a fitting round, and each
row to a diagnostic condition. Green tiles indicate passing conditions,
red tiles indicate failures. This visualization makes immediately
apparent, for instance, that NUTS fails on the ESS condition at round 1
but passes acceptance, while AMWG fails adaptation.

``` r

plot(cr_funnel, type = "refinement")
```

For the best method’s Consort diagnostics, the convergence scorecard
provides a single-number summary (0-100) combining all diagnostic
metrics through saturating transforms. This is useful for quick visual
comparison when examining many fits.

``` r

# Scorecard for the best method's consort diagnostics
plot(cr_funnel$Best$consort, type = "scorecard")
```

## Example 3: multi-family comparison

One of Crucible’s distinctive capabilities is comparing inference
methods across fundamentally different families. The following example
fits the normal-normal model using MCMC, variational Bayes, and Laplace
approximation simultaneously, letting Arena determine which method
provides the best efficiency-accuracy tradeoff.

``` r

# Compare across all families using specific methods
cr_multi <- Crucible(Model, Data, IV,
                     methods = c("NUTS", "AMWG", "Salimans2", "SPG"),
                     max_rounds = 2,
                     Chains = 2, CPUs = 1)

print(cr_multi)

# Arena comparison across families
plot(cr_multi, type = "arena")

# The Pareto frontier identifies methods that are not dominated
# by any other method on both speed and accuracy
plot(cr_multi, type = "dashboard")
```

For the simple normal-normal model, MCMC methods converge immediately at
round 1, while the Laplace and VB methods typically require additional
tuning to achieve appeasement. Although the posterior is genuinely
Gaussian for this model (making Laplace theoretically ideal), the
numerical Hessian computed at the mode can be poorly conditioned for
some optimizers, causing the positive-definiteness check to fail. NUTS
and AMWG produce reliable posterior samples with high ESS. The `methods`
argument lets the analyst explicitly choose which algorithms to compare,
bypassing the automatic Prescribe ranking when cross-family comparison
is desired. For a more complex model like the funnel, Laplace would fail
catastrophically (the Hessian at the mode does not capture the funnel
geometry), VB would struggle with the non-Gaussian tails, and only MCMC
methods would produce reliable results.

## Configuration

### Choosing `n_methods`

The `n_methods` argument controls how many top-ranked methods from
Prescribe are attempted. The default of 5 provides a broad comparison at
the cost of longer total runtime. For quick exploratory work,
`n_methods = 2` suffices; for a thorough benchmarking study,
`n_methods = 8` or more can be used. Each additional method adds one
full fitting cycle (including parallel chains for MCMC) to the total
time, so the tradeoff is linear.

### Choosing `max_rounds`

The `max_rounds` argument controls the maximum number of fitting
attempts per method. The first round uses the default specifications
from Crucible’s internal registry. Subsequent rounds use Consort’s
suggestions, which may increase iteration counts, adjust tuning
parameters, or switch to a more capable algorithm within the same
subcategory. Setting `max_rounds = 1` disables refinement entirely,
which is useful when the analyst wants a quick comparison without
iterative tuning. Setting `max_rounds = 3` (the default) allows two
additional refinement attempts, which is sufficient for most problems.
Very difficult problems (highly multimodal, severely ill-conditioned)
may benefit from `max_rounds = 5`, but if a method has not converged
after 3 rounds, it is unlikely to converge at all without fundamental
changes to the model parameterization.

### Choosing `families`

The `families` argument filters which inference families to consider.
The default includes MCMC, VB, Laplace, and IQ. Setting
`families = "MCMC"` restricts to MCMC algorithms only, which is
appropriate when the analyst specifically needs MCMC samples (e.g., for
complex posterior functionals). Setting `families = c("VB", "Laplace")`
restricts to fast approximation methods, useful for rapid model
iteration during development.

### Diversified method selection

The `diverse` argument (default `TRUE`) controls how Crucible selects
methods from Prescribe’s ranking. When `diverse = TRUE`, the selection
algorithm ensures representation from different MCMC subcategories
rather than picking the top-N by score, which could produce five
gradient-based methods on a gradient-available problem. The
diversification proceeds in three steps: first, it seeds one method per
inference family (MCMC, VB, Laplace, IQ); second, within MCMC, it seeds
the top-scoring method from each subcategory (gradient, ensemble, slice,
geometric, multimodal, etc.); third, it fills any remaining slots from
the global ranking. This portfolio approach increases the probability
that at least one selected method is well-suited to the problem’s true
structure, which may not be fully captured by the profiling phase.
Setting `diverse = FALSE` reverts to the legacy behavior of pure
score-based selection. The effect is visible in the selected methods
table:

``` r

# With diverse = TRUE (default), Crucible selects from different subcategories
cr_diverse <- Crucible(Model, Data, IV,
                       n_methods = 5, diverse = TRUE,
                       max_rounds = 1, Chains = 2, CPUs = 1,
                       verbose = FALSE)

# Inspect which methods were selected and their subcategories
cr_diverse$Methods[, c("Method", "Category", "Subcategory", "Score")]
```

### Parallel chains and CPUs

MCMC methods are always run with `Chains` parallel chains (default 4)
using `min(CPUs, Chains)` workers. Multiple chains are essential for
reliable convergence diagnostics: the split \\\hat{R}\\ statistic
requires at least 2 chains, and 4 chains provide more robust estimates
of between-chain variance [\[2\]](#ref2). Crucible runs methods
sequentially (one method at a time) to avoid nested parallelism, which
can cause resource contention and unpredictable performance on most
operating systems.

## Interpreting results

### The summary table

The `Summary` data frame is the primary output for comparing methods.
Each row represents one method with the following columns:

| Column   | Description                                                  |
|:---------|:-------------------------------------------------------------|
| Method   | Algorithm abbreviation (e.g., NUTS, AMWG, Salimans2)         |
| Category | Inference family (MCMC, VB, Laplace, IQ)                     |
| Rounds   | Rounds used / max_rounds                                     |
| Appeased | Whether Consort was appeased (all diagnostics pass)          |
| ESS.min  | Minimum bulk ESS across parameters (MCMC only)               |
| Minutes  | Total fitting time for this method                           |
| Rank     | Arena rank (1 = best) or ESS-based rank if Arena unavailable |

The Rank column reflects Arena’s composite score when multiple methods
succeed (combining efficiency and accuracy metrics), or a simpler
ordering by appeasement status and ESS when only one method succeeds or
Arena is skipped.

### Extracting results

The `crucible` object provides structured access to all pipeline
outputs.

``` r

# Access the Prescribe prescription
cr$Prescription$Profile$dim_class   # dimensionality class
cr$Prescription$Profile$speed_class # evaluation speed

# Access per-method results (first method as example)
first_method <- names(cr$Results)[1]
cat("First method:", first_method, "\n")
cr$Results[[first_method]]$appeased
cr$Results[[first_method]]$n_rounds
cr$Results[[first_method]]$consort$conditions

# Access Arena metrics (if available)
if (!is.null(cr$Arena))
    cr$Arena$Metrics[, c("Method", "ESS.per.second", "Rank.Overall")]

# Configuration and timing
cr$Config$families
cat("Total pipeline time:", round(cr$Minutes * 60, 1), "seconds\n")
```

### When methods fail

Not all methods succeed on all problems. Crucible wraps each fitting
call in error handling, so individual failures do not abort the
pipeline. Failed methods appear in the Summary table with `Rounds = 0`,
`Appeased = FALSE`, and `ESS.min = NA`. The `Results` list for a failed
method contains the error message in the `error` field. Common failure
modes include convergence failure for Laplace on multimodal posteriors,
numerical overflow for VB on poorly scaled models, and timeout for slow
ensemble methods on high-dimensional problems.

If all methods fail, Crucible returns a `crucible` object with an empty
`Best` field and a warning. This typically indicates a problem with the
model specification or initial values rather than with the inference
methods.

## Consort v2: design and architecture

### From monolith to modular dispatch

The original
[`Consort()`](https://robustecologies.github.io/lucifer/reference/Consort.md)
was a single 1158-line function containing hardcoded `if/else` chains
for each of the 55 MCMC algorithms. Each algorithm’s diagnostic logic
was scattered throughout the function, making it difficult to maintain,
extend, or test. Adding a new algorithm required modifying the function
in multiple places, and the algorithm-specific knowledge was duplicated
between Consort and the main
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
dispatcher.

The modernized Consort uses a two-level dispatch architecture. The
top-level
[`Consort()`](https://robustecologies.github.io/lucifer/reference/Consort.md)
function inspects the class of the fit object and dispatches to a
family-specific diagnostic engine: `.consort_mcmc()` for demonoid
objects, `.consort_vb()` for VB fits, `.consort_laplace()` for Laplace
fits, and so on. Within each family engine, algorithm-specific behavior
is driven by the `algo_registry` rather than hardcoded conditionals. For
example, the acceptance rate check looks up the algorithm’s expected
acceptance range from the registry (`algo_entry$acceptance_range`), and
the diminishing adaptation check consults the `has_adaptive_phase` flag.
This means adding a new MCMC algorithm to the registry automatically
makes it diagnosable by Consort without any changes to the diagnostic
code.

### The suggestion engine

When a fit is not appeased, the suggestion engine constructs a
structured recommendation containing four components: an `action` string
describing what to do (“increase iterations”, “reduce step size”,
“switch to AFSS”), a `next_specs` list that can be programmatically
passed to the fitting function, a human-readable `code` string that the
analyst can paste into a script, and a `rationale` explaining why this
suggestion was chosen.

The MCMC suggestion engine follows a prioritized decision tree.
Divergent transitions take highest priority because they indicate
systematic bias that cannot be fixed by running longer. Adaptation
failures are second, because a failed adaptive phase invalidates all
subsequent diagnostics. The non-adaptive phase check is third: many
adaptive algorithms (AMWG, NUTS, HMCDA) require a non-adaptive phase
after warmup to ensure the stationarity of the chain, and running
without it is a common oversight. Acceptance rate problems are fourth,
typically addressed by adjusting the proposal scale or escalating to a
more capable algorithm. Finally, if all other conditions pass but ESS or
\\\hat{R}\\ are inadequate, the suggestion is simply to increase
iterations, with the multiplier computed from the ratio of target ESS to
current ESS.

Key algorithm escalation paths are preserved from the original Consort:

- ADMG/AMWG \\\rightarrow\\ AFSS: when componentwise adaptive methods
  fail, the universal full-conditional slice sampler is a robust
  fallback that works for any continuously differentiable posterior.
- HMC \\\rightarrow\\ AHMC \\\rightarrow\\ NUTS: gradient-based
  escalation from fixed-trajectory HMC through adaptive HMC to the fully
  adaptive No-U-Turn Sampler.
- AM/AMM/DRAM/RAM \\\rightarrow\\ AMM \\\rightarrow\\ AMWG: adaptive
  random-walk escalation from block proposals to componentwise
  proposals.
- DEMC \\\rightarrow\\ AIES: ensemble method switch when differential
  evolution struggles.

### The `is.appeased()` function

The
[`is.appeased()`](https://robustecologies.github.io/lucifer/reference/is.appeased.md)
function has been rewritten from a fragile implementation based on
[`capture.output()`](https://rdrr.io/r/utils/capture.output.html) plus
[`grep()`](https://rdrr.io/r/base/grep.html) (which parsed the console
output of the old Consort for the string “has been appeased”) to a clean
implementation that calls `Consort(x, verbose = FALSE)$appeased`. This
approach is robust to changes in console formatting, works for all fit
classes (not just demonoid), and avoids the overhead of capturing and
parsing text output.

## Advanced usage

### Combining Crucible with model_spec

For models defined via the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
interface, the `spec` argument provides a convenient alternative to
passing Model, Data, and Initial.Values separately.

``` r

# model_spec uses the lucifer formula DSL to define models declaratively
ms <- model_spec("y ~ normal(mu, sigma); mu ~ normal(0, 100); sigma ~ halfcauchy(25)",
                 data = Data)

# Crucible extracts Model, Data, and Initial.Values from the spec
cr <- Crucible(spec = ms,
               n_methods = 3,
               max_rounds = 2)
```

### Passing extra arguments to Prescribe and Arena

The `prescribe.args` and `arena.args` arguments accept named lists that
are forwarded to
[`Prescribe()`](https://robustecologies.github.io/lucifer/reference/Prescribe.md)
and
[`Arena()`](https://robustecologies.github.io/lucifer/reference/Arena.md)
respectively. This allows fine-tuning the profiling phase (e.g.,
increasing `n.profile` for a more thorough multimodality assessment) or
the comparison phase (e.g., supplying an analytical reference
distribution for accuracy measurement).

``` r

# Increase profiling thoroughness and supply reference
cr <- Crucible(Model, Data, IV,
               methods = c("NUTS", "AMWG", "Slice"),
               max_rounds = 2,
               Chains = 2, CPUs = 1,
               prescribe.args = list(n.profile = 500),
               arena.args = list(verbose = FALSE))
```

### Iterating after Crucible

Crucible’s result is a starting point, not a final answer. The analyst
should always inspect the Best fit critically, especially for
high-stakes analyses. If the best method’s Consort is not appeased
(which can happen if `max_rounds` was insufficient), the analyst can
continue manually:

``` r

# Check if best method is appeased
if (!cr$Best$consort$appeased) {
    # Get the suggestion
    suggestion <- cr$Best$consort$suggestion

    # View the suggested code
    cat(suggestion$code, "\n")

    # Or extract the structured specs for programmatic use
    next_specs <- suggestion$next_specs
}
```

## Diagnostic thresholds summary

The following table summarizes the diagnostic thresholds used by
[`Consort()`](https://robustecologies.github.io/lucifer/reference/Consort.md)
across inference families, with references to the literature
establishing each threshold.

| Family | Condition | Threshold | Reference |
|----|----|----|----|
| MCMC | \\\hat{R}\\ | \\\< 1.01\\ | Vehtari et al. (2021) [\[2\]](#ref2) |
| MCMC | ESS bulk | \\\geq 400\\ | Vehtari et al. (2021) [\[2\]](#ref2) |
| MCMC | ESS tail | \\\geq 200\\ | Vehtari et al. (2021) [\[2\]](#ref2) |
| MCMC | MCSE | \\\< 6.27\\\\ of SD | Flegal et al. (2008) [\[4\]](#ref4) |
| MCMC | Divergences | \\= 0\\ | Betancourt (2017) [\[5\]](#ref5) |
| VB | ELBO stability | CV \\\< 5\\\\ in final 10% | Kucukelbir et al. (2017) [\[6\]](#ref6) |
| VB | Pareto \\k\\ | \\\< 0.7\\ | Yao et al. (2018) [\[7\]](#ref7) |
| Laplace | Gradient norm | \\\< 10^{-3}\\ | Tierney and Kadane (1986) [\[8\]](#ref8) |
| Laplace | Hessian | Positive-definite | Tierney and Kadane (1986) [\[8\]](#ref8) |

## References

**\[1\]** Rice, J.R. (1976). *The algorithm selection problem*. Advances
in Computers, 15, 65-118.
[doi:10.1016/S0065-2458(08)60520-3](https://doi.org/10.1016/S0065-2458(08)60520-3)

**\[2\]** Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and
Burkner, P.-C. (2021). Rank-normalization, folding, and localization: an
improved \\\hat{R}\\ for assessing convergence of MCMC. *Bayesian
Analysis*, 16(2), 667-718.
[doi:10.1214/20-BA1221](https://doi.org/10.1214/20-BA1221)

**\[3\]** Gelman, A. and Rubin, D.B. (1992). Inference from iterative
simulation using multiple sequences. *Statistical Science*, 7(4),
457-472.
[doi:10.1214/ss/1177011136](https://doi.org/10.1214/ss/1177011136)

**\[4\]** Flegal, J.M., Haran, M. and Jones, G.L. (2008). Markov chain
Monte Carlo: can we trust the third significant figure? *Statistical
Science*, 23(2), 250-260.
[doi:10.1214/08-STS257](https://doi.org/10.1214/08-STS257)

**\[5\]** Betancourt, M. (2017). A conceptual introduction to
Hamiltonian Monte Carlo. *arXiv:1701.02434*.
[arXiv:1701.02434](https://arxiv.org/abs/1701.02434)

**\[6\]** Kucukelbir, A., Tran, D., Ranganath, R., Gelman, A. and Blei,
D.M. (2017). Automatic differentiation variational inference. *Journal
of Machine Learning Research*, 18(14), 1-45.
[JMLR](http://jmlr.org/papers/v18/16-107.md)

**\[7\]** Yao, Y., Vehtari, A., Simpson, D. and Gelman, A. (2018). Yes,
but did it work? Evaluating variational inference. *Proceedings of the
35th International Conference on Machine Learning*, PMLR 80, 5581-5590.
[arXiv:1802.02538](https://arxiv.org/abs/1802.02538)

**\[8\]** Tierney, L. and Kadane, J.B. (1986). Accurate approximations
for posterior moments and marginal densities. *Journal of the American
Statistical Association*, 81(393), 82-86.
[doi:10.1080/01621459.1986.10478240](https://doi.org/10.1080/01621459.1986.10478240)

**\[9\]** Goodman, J. and Weare, J. (2010). Ensemble samplers with
affine invariance. *Communications in Applied Mathematics and
Computational Science*, 5(1), 65-80.
[doi:10.2140/camcos.2010.5.65](https://doi.org/10.2140/camcos.2010.5.65)

**\[10\]** ter Braak, C.J.F. (2006). A Markov chain Monte Carlo version
of the genetic algorithm differential evolution: easy Bayesian computing
for real parameter spaces. *Statistics and Computing*, 16(3), 239-249.
[doi:10.1007/s11222-006-8769-1](https://doi.org/10.1007/s11222-006-8769-1)

**\[11\]** Betancourt, M. and Girolami, M. (2015). Hamiltonian Monte
Carlo for hierarchical models. In S.K. Upadhyay, U. Singh, D.K. Dey, and
A. Loganathan (Eds.), *Current Trends in Bayesian Methodology with
Applications*, 79-101. CRC Press.
[arXiv:1312.0906](https://arxiv.org/abs/1312.0906)
