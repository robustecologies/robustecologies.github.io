# Bayesian model selection with the product space method

  

## Introduction

Bayesian model selection requires computing marginal likelihoods \\p(y
\mid \mathcal{M}\_k) = \int p(y \mid \theta_k, \mathcal{M}\_k) \\
p(\theta_k \mid \mathcal{M}\_k) \\ d\theta_k\\, high-dimensional
integrals that are notoriously difficult to evaluate. The ratio of
marginal likelihoods for two models gives the Bayes factor \\B\_{ab} =
p(y \mid \mathcal{M}\_a) / p(y \mid \mathcal{M}\_b)\\, which quantifies
the relative evidence the data provide for one model over another,
independent of prior model probabilities. Standard approaches to
estimating Bayes factors include the harmonic mean estimator (unstable
and unbounded variance), the Laplace-Metropolis approximation (requires
unimodal posteriors), bridge sampling (requires careful design of the
bridge distribution), and sequential Monte Carlo (reliable but
computationally expensive for many models).

The product space method, introduced by Carlin and Chib [\[1\]](#ref1),
offers an alternative that is conceptually straightforward and operates
entirely within the standard MCMC framework. Rather than estimating
marginal likelihoods separately for each model, it constructs a single
“supermodel” that encompasses all \\K\\ candidate models simultaneously.
A discrete model indicator \\M \in \\1, \ldots, K\\\\ selects the active
model at each MCMC iteration, and the posterior proportion of time spent
in each model directly estimates the posterior model probabilities. The
[`ProductSpace()`](https://robustecologies.github.io/lucifer/reference/ProductSpace.md)
function in **lucifer** implements this method with automatic
pseudoprior estimation from pilot MCMC fits, bisection calibration for
balanced model activation [\[2\]](#ref2), posterior model space
complexity metrics, and a hybrid product-space sampler that combines
Robust Adaptive Metropolis (RAM, [\[6\]](#ref6)) for active model
parameters with exact enumeration for the model indicator and direct
sampling from pseudopriors for inactive parameters.

  

## Theory

### The product space construction

Suppose we wish to compare \\K\\ models \\\mathcal{M}\_1, \ldots,
\mathcal{M}\_K\\, each with its own parameter vector \\\theta_k\\ of
dimension \\d_k\\. The product space method, introduced by Carlin and
Chib [\[1\]](#ref1) and subsequently refined by Godsill [\[7\]](#ref7)
and Dellaportas, Forster, and Ntzoufras [\[8\]](#ref8), augments the
parameter space by concatenating all model-specific parameter vectors
with a discrete model indicator:

\\ \boldsymbol{\theta} = (\theta_1, \theta_2, \ldots, \theta_K, M),
\qquad M \in \\1, \ldots, K\\. \\

The dimension of this augmented space is \\d = \sum\_{k=1}^K d_k + 1\\,
which grows linearly in the total number of parameters across all
models. The construction rests on a simple but powerful idea: rather
than integrating out each model’s parameters separately to obtain
marginal likelihoods, we define a joint distribution over the entire
concatenated parameter space and let the MCMC sampler explore it. The
proportion of time the chain spends in each model directly estimates the
posterior model probability, circumventing the need for explicit
marginal likelihood computation.

The augmented joint posterior is defined so that when model \\M = k\\ is
active, the likelihood and prior of model \\k\\ determine the data fit,
while the parameters of all other models \\j \neq k\\ are drawn from
pseudoprior distributions \\q(\theta_j \mid M = k)\\. The target density
factorizes as:

\\ \pi(\theta_1, \ldots, \theta_K, M{=}k \mid y) \propto \underbrace{p(y
\mid \theta_k, \mathcal{M}\_k) \\ p(\theta_k \mid
\mathcal{M}\_k)}\_{\text{active model}} \times \underbrace{\prod\_{j
\neq k} q(\theta_j \mid M{=}k)}\_{\text{pseudopriors (known)}} \times
p(M{=}k) \\

This factorization reveals the essential structure of the product space.
The likelihood involves only the active model’s parameters \\\theta_k\\;
the inactive parameters \\\theta_j\\ for \\j \neq k\\ appear exclusively
in their pseudoprior terms, which are known densities chosen by the
analyst. The inactive parameters are therefore **conditionally
independent** of the data and of the active parameters given the model
indicator \\M = k\\. Their full conditional distribution under the
product space target is exactly the pseudoprior: \\\theta_j \mid
\text{rest}, M{=}k \sim q(\theta_j \mid M{=}k)\\ for all \\j \neq k\\.
This conditional independence is not an approximation but an exact
consequence of the product space construction, and it has profound
implications for sampling efficiency.

Three consequences follow directly from this structure. First, the
inactive parameters can be **sampled directly** from their pseudoprior
distributions at each iteration, producing exact, independent draws with
zero autocorrelation and zero model evaluations. There is no need to run
an MCMC kernel for these parameters, since their conditional
distribution is a known Gaussian from which we can draw immediately.
Second, the active model’s parameters \\\theta_k\\ can be updated with
**any valid MCMC kernel** targeting the model-specific posterior
\\p(\theta_k \mid y, \mathcal{M}\_k)\\. In a Metropolis-Hastings update
for \\\theta_k\\, the pseudoprior terms and the model prior \\p(M{=}k)\\
cancel in the acceptance ratio because they depend only on the
(unchanged) inactive parameters and the model indicator, both of which
are held fixed during the active-parameter update. This opens the door
to efficient block proposals, adaptive methods, or gradient-based
samplers for the active model. Third, the model indicator \\M\\ must be
sampled via discrete enumeration (a Gibbs step), which requires
evaluating the log-posterior at each of the \\K\\ candidate models. This
step is \\O(K)\\ and unavoidable, but for moderate \\K\\ it is typically
the cheap part of the computation. The
[`ProductSpace()`](https://robustecologies.github.io/lucifer/reference/ProductSpace.md)
default hybrid sampler exploits all three of these consequences; the
section “The hybrid sampler” below describes the algorithm in detail.

The **cost reduction** from exploiting conditional independence is
substantial. A naive Gibbs sampler that updates every continuous
parameter via univariate slice sampling (Neal, [\[4\]](#ref4)) performs
\\\sim 10 \sum d_k + K\\ model evaluations per iteration (the factor of
10 accounts for the stepping-out and shrinkage phases of slice
sampling). The hybrid sampler reduces this to \\K + 1\\: \\K\\
evaluations for the discrete enumeration of \\M\\, and 1 evaluation for
the block RAM proposal on the active model. For the seven-model time
series example in this vignette (\\\sum d_k = 22\\), this is a \\\sim 30
\times\\ reduction in model evaluations per iteration.

The pseudoprior distributions are proper densities that integrate to 1.
Crucially, they do not affect the marginal likelihood under any model,
because they factor out of the integral over the inactive parameters:

\\ p(y \mid \mathcal{M}\_k) = \int p(y \mid \theta_k, \mathcal{M}\_k) \\
p(\theta_k \mid \mathcal{M}\_k) \\ d\theta_k \times
\underbrace{\prod\_{j \neq k} \int q(\theta_j \mid M{=}k) \\
d\theta_j}\_{= 1}. \\

Because the pseudopriors integrate to unity, they vanish from the
marginal \\p(y \mid \mathcal{M}\_k)\\, and the posterior model odds
computed from the product space are exactly the true Bayes factor times
the prior model odds:

\\ \frac{p(M{=}a \mid y)}{p(M{=}b \mid y)} = B\_{ab} \times
\frac{p(M{=}a)}{p(M{=}b)}. \\

This identity holds regardless of the choice of pseudoprior densities,
so long as they are proper. The pseudopriors affect only the mixing
efficiency of the MCMC sampler, not the validity of the Bayes factor
estimates (see the next subsection). Carlin and Chib [\[1\]](#ref1)
originally proved this result for two competing models; the extension to
arbitrary \\K\\ is straightforward and follows from the same
marginalization argument. Godsill [\[7\]](#ref7) provided an alternative
derivation showing that the product space method arises naturally as a
special case of a general framework for composite model spaces, where
each model contributes a “component” to a mixture target and the model
indicator selects among components.

### Pseudopriors and mixing

Although the choice of pseudoprior distributions has no effect on the
theoretical Bayes factor, it is the single most important determinant of
MCMC mixing in the product space. The mechanism is straightforward: when
the chain proposes a transition from model \\k\\ to model \\j\\, the
parameters \\\theta_j\\ must take on values that are plausible under
model \\j\\’s posterior, since the proposed state will be evaluated
against the likelihood \\p(y \mid \theta_j, \mathcal{M}\_j)\\. If the
pseudoprior \\q(\theta_j \mid M{=}k)\\ concentrates mass in a region far
from the posterior \\p(\theta_j \mid y, \mathcal{M}\_j)\\, the
log-posterior at the proposed state will be very low, and the model
switch will be rejected. This causes the model indicator to “stick” at
one value for long stretches, producing poor Bayes factor estimates. In
the extreme case where the pseudoprior and the posterior have negligible
overlap, the chain will never switch models, making the Bayes factor
estimate meaningless regardless of how long the chain runs.

The optimal pseudoprior for mixing would be exactly the model-specific
posterior, \\q^\*(\theta_j \mid M{=}k) = p(\theta_j \mid y,
\mathcal{M}\_j)\\, but this is of course unknown; if we knew the
posteriors exactly, we would not need MCMC. The standard recommendation,
established by Carlin and Chib [\[1\]](#ref1) and elaborated by
Lodewyckx et al. [\[2\]](#ref2), is to use moment-matched Gaussian
approximations: run each model separately in a pilot MCMC phase, compute
the posterior mean \\\hat{\mu}\_j\\ and standard deviation
\\\hat{\sigma}\_j\\ for each parameter, and set \\q(\theta_j \mid M{=}k)
= \mathcal{N}(\hat{\mu}\_j, \hat{\sigma}\_j^2)\\ independently for each
component of \\\theta_j\\. The independence assumption (diagonal
covariance) is a simplification, but it suffices for moment matching
because the primary goal is to place the pseudoprior’s mass in the right
region of parameter space, not to capture the exact posterior geometry.

A mild inflation of the pseudoprior variance is standard practice.
[`ProductSpace()`](https://robustecologies.github.io/lucifer/reference/ProductSpace.md)
multiplies the posterior standard deviation by a factor of 1.5
(configurable via `Pseudoprior.Inflation`), producing pseudopriors that
are slightly overdispersed relative to the true posteriors. This
asymmetry is deliberate: a pseudoprior that is too narrow misses the
tails of the posterior and causes rejected model switches, while a
pseudoprior that is too wide wastes some proposals in low-density
regions but does not prevent switching entirely. Dellaportas, Forster,
and Ntzoufras [\[8\]](#ref8) discuss alternative pseudoprior
constructions, including kernel density estimates of the pilot posterior
and mixtures of normals for multimodal posteriors. In practice, the
simple moment-matched Gaussian with mild inflation performs well for the
unimodal posteriors typical of well-identified statistical models.

[`ProductSpace()`](https://robustecologies.github.io/lucifer/reference/ProductSpace.md)
automates the entire pseudoprior estimation process: it runs \\K\\
independent pilot MCMC fits (one per model) using the algorithm
specified by `Pilot.Algorithm` (default `"CHARM"`), extracts posterior
summaries from the second half of each chain, and constructs the
pseudopriors automatically. The pilot fits run in parallel via
[`callr::r_bg()`](https://callr.r-lib.org/reference/r_bg.html) when
`CPUs > 1`. Users who wish to supply their own pseudopriors, for example
from a previous analysis or from an analytical approximation, can do so
via the `Pseudopriors` argument, bypassing the pilot phase entirely.

### Bisection calibration

Even with well-constructed pseudopriors, the raw posterior model
probabilities can be extremely unbalanced when one model is strongly
preferred by the data. Consider a case where \\\hat{p}(M{=}1 \mid y) =
0.999\\ and \\\hat{p}(M{=}2 \mid y) = 0.001\\: the chain visits model 2
for only 0.1% of its iterations, producing a Bayes factor estimate based
on a very small effective sample size. The Bayes factor itself may be
accurate in sign (model 1 is preferred), but the magnitude \\B\_{12}\\
will have high Monte Carlo variance because it depends on the ratio of
small estimated probabilities.

Lodewyckx et al. [\[2\]](#ref2) addressed this problem with a bisection
algorithm that adjusts the prior model probabilities \\\mathbf{p} =
(p(M{=}1), \ldots, p(M{=}K))\\ to equalize posterior model activation.
The idea is to inflate the prior probability of disfavored models and
deflate that of favored models, so that the MCMC chain spends
approximately equal time in each model. After the production MCMC run,
the corrected posterior probabilities are recovered analytically by
dividing out the calibrated priors:

\\ \tilde{p}(M{=}k \mid y) = \frac{\hat{p}(M{=}k \mid y) /
p(M{=}k)}{\sum\_{j=1}^K \hat{p}(M{=}j \mid y) / p(M{=}j)}. \\

This correction is exact and does not introduce bias, because it is
simply Bayes’ rule applied in reverse: the model probabilities produced
by the chain under modified priors are converted back to what they would
have been under the original (usually uniform) priors. The bisection
algorithm iterates: run a short MCMC chain, measure the posterior model
activation frequencies, adjust the prior probabilities to push
underrepresented models toward equal activation, and repeat. Convergence
is monitored by checking whether all models are activated within
`Bisection.Tolerance` (default 0.10) of the uniform \\1/K\\. The
algorithm typically converges in 3–5 rounds for moderate evidence
imbalances.

[`ProductSpace()`](https://robustecologies.github.io/lucifer/reference/ProductSpace.md)
implements bisection calibration by default (`Bisection = TRUE`), with
configurable parameters: `Bisection.Iterations` controls the length of
each calibration chain, `Bisection.Max.Rounds` caps the number of
bisection iterations, and `Bisection.Tolerance` sets the convergence
criterion. When the evidence is so overwhelming that no amount of prior
rebalancing can force the chain to visit the disfavored model (e.g.,
comparing white noise against a strongly autocorrelated AR(1) process),
the algorithm detects the failure, defined as zero model switching for 3
consecutive rounds, terminates bisection early, and issues a warning. In
such cases, the resulting Bayes factors should be interpreted as lower
bounds on the true evidence.

### Posterior model space complexity

Given the corrected posterior model probabilities \\\tilde{p}\_1,
\ldots, \tilde{p}\_K\\,
[`ProductSpace()`](https://robustecologies.github.io/lucifer/reference/ProductSpace.md)
computes five information-theoretic metrics that characterize the
complexity of the posterior model space. These metrics quantify how
concentrated the posterior is on a single model versus spread across
many, providing a scalar summary of model selection uncertainty.

**Shannon entropy.** The Shannon entropy \\H = -\sum\_{k=1}^K
\tilde{p}\_k \ln \tilde{p}\_k\\ measures the total uncertainty in the
posterior model distribution. It ranges from 0 (all probability on one
model) to \\\ln K\\ (uniform across all models). The normalized version
\\H / \ln K\\ rescales this to \\\[0, 1\]\\ for comparability across
analyses with different numbers of candidate models.

**Effective number of models.** The quantity \\K\_{\mathrm{eff}} =
\exp(H)\\ gives the Hill number of order 1, interpretable as the number
of equiprobable models that would produce the same Shannon entropy. It
ranges from 1 (complete dominance) to \\K\\ (uniform), and provides an
intuitive count-like summary.

**Simpson diversity.** The index \\D = 1 - \sum\_{k=1}^K
\tilde{p}\_k^2\\ gives the probability that two independently drawn MCMC
iterations come from different models. It ranges from 0 to \\(K-1)/K\\
and is less sensitive to rare models than the Shannon entropy.

**Gini effective number.** The reciprocal of the Herfindahl-Hirschman
index, \\1 / \sum\_{k=1}^K \tilde{p}\_k^2\\, is the Hill number of order
2. It provides a concentration-weighted count of effective models, more
influenced by the dominant model than the Shannon-based
\\K\_{\mathrm{eff}}\\.

### Relationship to other methods

The product space method belongs to a family of transdimensional MCMC
techniques that enable posterior inference over collections of models
with different parameter dimensions. The most prominent member of this
family is reversible-jump MCMC (RJMCMC), introduced by Green
[\[5\]](#ref5), which constructs dimension-matching bijections between
parameter spaces of different models and uses Metropolis-Hastings
proposals that jump between model-dimension pairs. Godsill
[\[7\]](#ref7) showed that the product space method and RJMCMC are
special cases of the same general framework for composite model spaces;
the difference lies in how the inactive parameters are handled. In the
product space, inactive parameters are explicitly represented and
sampled from pseudopriors, while in RJMCMC they are marginalized out via
the dimension-matching transformation. RJMCMC avoids sampling
unnecessary parameters, which can be more efficient for problems where
\\\sum d_k\\ is large, but it requires designing bijective mappings and
computing Jacobian determinants for each pair of models, a task that
becomes increasingly difficult as the models grow more dissimilar in
structure. The product space method requires no such mappings; the only
design choice is the pseudoprior, which can be estimated automatically
from pilot fits.

Dellaportas, Forster, and Ntzoufras [\[8\]](#ref8) provided a systematic
comparison of the product space method, RJMCMC, and Gibbs variable
selection (the indicator-based approach of George and McCulloch, 1993)
in the context of Bayesian variable selection for linear models. They
found that the product space method performed comparably to RJMCMC when
pseudopriors were well-calibrated, and noted that its main advantage was
ease of implementation within standard MCMC software. The hybrid sampler
implemented in
[`ProductSpace()`](https://robustecologies.github.io/lucifer/reference/ProductSpace.md)
mitigates the historical cost of sampling inactive parameters (the
traditional weakness of the product space approach) by direct-sampling
them from pseudopriors rather than running MCMC kernels, recovering much
of RJMCMC’s efficiency advantage while retaining the product space
method’s simplicity.

Bridge sampling [\[9\]](#ref9) and its generalizations provide another
class of Bayes factor estimators that work by constructing an optimal
importance sampling bridge between the posterior distributions of two
models. Bridge sampling can be very efficient when the bridge
distribution is well-chosen, but it estimates only pairwise Bayes
factors and requires separate posterior samples from each model. For
\\K\\ models, this gives \\K(K-1)/2\\ pairwise comparisons, each
requiring its own bridge construction. The product space method, by
contrast, produces all pairwise Bayes factors simultaneously from a
single MCMC run. Sequential Monte Carlo (SMC) methods, available via
[`SMC()`](https://robustecologies.github.io/lucifer/reference/SMC.md) in
lucifer, provide marginal likelihood estimates as a byproduct of the
tempering schedule and can handle model comparison naturally. SMC avoids
the pseudoprior design problem entirely, since it works by tempering
from the prior to the posterior within each model separately. For
problems with many models or where pseudoprior calibration is difficult,
SMC may be more reliable, though at higher computational cost per model.
Han and Carlin [\[10\]](#ref10) provide a comprehensive comparative
review of MCMC-based methods for computing Bayes factors, including the
product space method, RJMCMC, bridge sampling, and several importance
sampling variants. The
[`BayesFactor()`](https://robustecologies.github.io/lucifer/reference/BayesFactor.md)
function in lucifer computes Bayes factors from log-marginal likelihood
estimates produced by
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation()`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
or
[`VariationalBayes()`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md),
offering yet another route.

Information criteria such as
[`WAIC()`](https://robustecologies.github.io/lucifer/reference/WAIC.md)
and
[`LOO()`](https://robustecologies.github.io/lucifer/reference/LOO.md)
provide asymptotic approximations to model comparison that avoid both
the pseudoprior design problem and the explicit computation of marginal
likelihoods. These criteria are computationally cheap and asymptotically
consistent under regularity conditions, but they approximate the
expected log predictive density rather than the marginal likelihood, and
their finite-sample behavior can diverge from Bayes factor rankings for
small datasets or models with informative priors. Kass and Raftery
[\[11\]](#ref11) provide a comprehensive review of Bayes factors, their
computation, and their relationship to other model comparison criteria.

### Interpreting Bayes factors

The Bayes factor \\B\_{ab}\\ quantifies the evidence the data provide
for model \\\mathcal{M}\_a\\ over model \\\mathcal{M}\_b\\ on a
continuous scale, and it does so in a way that is independent of the
prior model probabilities. Jeffreys [\[12\]](#ref12) proposed a
widely-used calibration scale for interpreting the magnitude of Bayes
factors: values between 1 and 3 are considered “barely worth
mentioning,” between 3 and 10 as “substantial,” between 10 and 30 as
“strong,” between 30 and 100 as “very strong,” and above 100 as
“decisive.” These thresholds are conventions rather than formal decision
boundaries and should be interpreted in context; a Bayes factor of 5 may
be compelling in a screening study but insufficient for a confirmatory
analysis. Kass and Raftery [\[11\]](#ref11) modified Jeffreys’ original
scale and recommended reporting Bayes factors on the logarithmic scale
(\\2 \ln B\_{ab}\\), which is approximately on the same scale as
likelihood ratio test statistics and facilitates comparison with
frequentist model selection.
[`ProductSpace()`](https://robustecologies.github.io/lucifer/reference/ProductSpace.md)
reports both \\\log\_{10} B\_{ab}\\ and the Jeffreys interpretation in
its summary output, and the `plot(type = "bayes_factors")` method
annotates the forest plot with color-coded evidence bands.

  

## Algorithm and implementation

### Workflow

[`ProductSpace()`](https://robustecologies.github.io/lucifer/reference/ProductSpace.md)
proceeds in six stages:

1.  **Input validation.** Checks that all Model functions follow the
    lucifer contract (`function(parm, Data)` returning
    `list(LP, Dev, Monitor, yhat, parm)`), that Data lists contain
    `parm.names` and `mon.names`, and that dimensions are consistent.

2.  **Pilot fits.** Unless user-supplied pseudopriors are provided, the
    function runs \\K\\ independent pilot MCMC fits (one per model) in
    parallel using
    [`callr::r_bg()`](https://callr.r-lib.org/reference/r_bg.html). Each
    pilot uses the algorithm specified by `Pilot.Algorithm` (default
    `"CHARM"`) for `Pilot.Iterations` iterations. Posterior means and
    (inflated) standard deviations from these pilots become the
    pseudoprior parameters.

3.  **Supermodel construction.** The super-parameter vector
    \\\boldsymbol{\theta} = (\theta_1, \ldots, \theta_K, M)\\ is
    assembled by concatenating all model parameters and appending the
    model indicator. A super-Model function evaluates the active model’s
    log-posterior plus the pseudoprior contributions of all inactive
    models plus the log prior on \\M\\.

4.  **Bisection calibration (optional).** Short MCMC chains are run with
    iteratively adjusted prior model probabilities until posterior model
    activation is approximately balanced (within `Bisection.Tolerance`).
    The algorithm stops early if any model is never visited for 3
    consecutive rounds, indicating overwhelming evidence.

5.  **Production MCMC.** The `Production.Algorithm` argument controls
    which sampler is used:

    - `"RAM"` (default): a hybrid sampler that exploits the conditional
      independence structure of the product space. At each iteration,
      the model indicator \\M\\ is updated by exact enumeration over
      \\\\1, \ldots, K\\\\ (\\K\\ model evaluations), the active model’s
      parameters are updated with a block RAM proposal (1 model
      evaluation), and inactive parameters are drawn directly from their
      pseudopriors (0 model evaluations). Total cost: \\K + 1\\ model
      evaluations per iteration.
    - `"Gibbs"`: the original auto-FC Gibbs sampler with univariate
      stepping-out slice sampling for continuous parameters and exact
      enumeration for \\M\\. Cost: \\\sim 10 \sum d_k + K\\ model
      evaluations per iteration.

6.  **Post-processing.** The function extracts the model indicator
    trace, computes raw and corrected posterior model probabilities,
    pairwise Bayes factor matrices, per-model posterior summaries,
    transition matrix, switching rate, model indicator ESS, and
    complexity metrics.

### The hybrid sampler: hybrid product space sampler for scalable Bayesian model selection

The default `Production.Algorithm = "RAM"` implements a three-component
hybrid sampler that exploits the factorization of the product space
target. When \\M = k\\, the inactive parameters \\\theta_j\\ (\\j \neq
k\\) are conditionally independent of the data: their full conditional
is the pseudoprior \\q(\theta_j \mid M{=}k) = \mathcal{N}(\mu_j,
\sigma_j^2)\\. This means they can be drawn exactly, with zero
autocorrelation and zero model evaluations. The active model’s
parameters \\\theta_k\\ are updated with a block Robust Adaptive
Metropolis (RAM) proposal [\[6\]](#ref6), which adapts a full-rank
proposal covariance to target acceptance rate 0.234 (optimal for
Gaussian targets). Each model maintains its own adaptation state,
initialized from the pilot posterior covariance scaled by \\2.38^2 /
d_k\\. When the model indicator switches from model \\a\\ to model
\\b\\, the sampler resumes model \\b\\’s adaptation state where it left
off.

The cost reduction is substantial. The original Gibbs sampler performs
\\\sim 10\\ model evaluations per continuous parameter (stepping-out +
shrinkage in slice sampling) plus \\K\\ for the model indicator, giving
\\\sim 10 \sum d_k + K\\ per iteration. The hybrid sampler performs
\\K + 1\\: \\K\\ for enumerating \\M\\ and 1 for the RAM proposal on the
active model. For the seven-model time series example in this vignette
(\\\sum d_k = 22\\), this is a \\\sim 30\\x reduction in model
evaluations per iteration.

### Return value

[`ProductSpace()`](https://robustecologies.github.io/lucifer/reference/ProductSpace.md)
returns an S3 object of class `product_space` with the following
components:

- `Posterior.Model.Probabilities`: corrected posterior probabilities.
- `Bayes.Factors`, `Log.Bayes.Factors`: pairwise Bayes factor matrices.
- `Per.Model.Posteriors`, `Per.Model.Summaries`: model-conditional
  posteriors.
- `Model.Indicator`: integer trace of \\M\\ across iterations.
- `Transition.Matrix`: \\K \times K\\ model-switch probabilities.
- `Switching.Rate`, `M.ESS`: mixing diagnostics.
- `Complexity`: list with Shannon entropy, effective K, Simpson
  diversity.
- `Pilot.Fits`, `Pseudopriors`, `Bisection.History`: intermediate
  objects.
- `MCMC.Fit`: the raw demonoid object from the production MCMC.

### S3 methods

- [`print()`](https://rdrr.io/r/base/print.html): one-screen summary
  with probabilities, Bayes factors, complexity metrics, and mixing
  diagnostics.
- [`summary()`](https://rdrr.io/r/base/summary.html): extended output
  including full Bayes factor matrix with Jeffreys interpretation,
  transition matrix, per-model posterior summaries, and bisection
  history.
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html): seven plot
  types via the `type` argument:
  - `"all"` (default): combined 5-panel layout.
  - `"probabilities"`: bar chart with complexity annotation.
  - `"trace"`: model indicator trace colored by model.
  - `"transition"`: heatmap of transition matrix.
  - `"bayes_factors"`: forest plot with Jeffreys scale legend.
  - `"posteriors"`: overlapping density plots faceted by parameter.
  - `"complexity"`: bullet chart of model space complexity metrics.

  

## Example: time series model selection

### Data generation

We simulate a time series from a known deterministic process
contaminated by additive stochastic noise. The generative model is an
AR(1) process with intercept \\\alpha = 2\\, autoregressive coefficient
\\\phi_1 = 0.7\\, and innovation standard deviation \\\sigma = 1\\:

\\ y_t = \alpha + \phi_1 (y\_{t-1} - \alpha) + \varepsilon_t, \qquad
\varepsilon_t \sim \mathcal{N}(0, \sigma^2). \\

We will fit seven competing models, including the true AR(1), and use
[`ProductSpace()`](https://robustecologies.github.io/lucifer/reference/ProductSpace.md)
to identify which model the data support.

``` r

library(lucifer)
set.seed(666)

## True parameters
n <- 200
alpha_true <- 2
phi1_true  <- 0.7
sigma_true <- 0.1

## Simulate AR(1)
y <- numeric(n)
y[1] <- alpha_true + rnorm(1, 0, sigma_true / sqrt(1 - phi1_true^2))
for (t in 2:n) {
    y[t] <- alpha_true + phi1_true * (y[t - 1] - alpha_true) +
            rnorm(1, 0, sigma_true)
}

plot(y, type = "l", col = "#4e79a7", lwd = 1.2,
     xlab = "Time", ylab = expression(y[t]),
     main = "Simulated AR(1) time series")
abline(h = alpha_true, lty = 2, col = "grey50")
```

### Candidate models

We define seven models spanning different levels of structural
complexity. All models follow the lucifer contract and use weakly
informative priors.

**Model 1: White noise (AR(0)).** The simplest model assumes independent
observations with a common mean and variance. This model ignores all
temporal dependence.

``` r

Model_WN <- function(parm, Data) {
    mu    <- parm[1]
    sigma <- interval(parm[2], 1e-100, Inf)
    parm[2] <- sigma
    LL <- sum(dnorm(Data$y, mu, sigma, log = TRUE))
    LP <- LL + dnorm(mu, 0, 100, log = TRUE) +
          dhalfcauchy(sigma, 10, log = TRUE)
    list(LP = LP, Dev = -2 * LL, Monitor = LP,
         yhat = rnorm(Data$n, mu, sigma), parm = parm)
}
Data_WN <- list(y = y, n = n,
                parm.names = c("mu", "sigma"), mon.names = "LP")
```

**Model 2: AR(1).** The true generative model. Estimates an intercept,
one autoregressive coefficient, and the innovation variance.

``` r

Model_AR1 <- function(parm, Data) {
    mu   <- parm[1]
    phi1 <- interval(parm[2], -0.999, 0.999)
    parm[2] <- phi1
    sigma <- interval(parm[3], 1e-100, Inf)
    parm[3] <- sigma
    resid <- Data$y[-1] - mu - phi1 * (Data$y[-Data$n] - mu)
    LL <- sum(dnorm(resid, 0, sigma, log = TRUE))
    LP <- LL + dnorm(mu, 0, 100, log = TRUE) +
          dnorm(phi1, 0, 1, log = TRUE) +
          dhalfcauchy(sigma, 10, log = TRUE)
    list(LP = LP, Dev = -2 * LL, Monitor = LP,
         yhat = rnorm(Data$n, mu, sigma), parm = parm)
}
Data_AR1 <- list(y = y, n = n,
                 parm.names = c("mu", "phi1", "sigma"), mon.names = "LP")
```

**Model 3: AR(2).** Extends AR(1) with a second-order lag. Since the
true \\\phi_2 = 0\\, this model is overparameterized but nested with the
truth.

``` r

Model_AR2 <- function(parm, Data) {
    mu   <- parm[1]
    phi1 <- interval(parm[2], -0.999, 0.999)
    parm[2] <- phi1
    phi2 <- interval(parm[3], -0.999, 0.999)
    parm[3] <- phi2
    sigma <- interval(parm[4], 1e-100, Inf)
    parm[4] <- sigma
    nn <- Data$n
    resid <- Data$y[3:nn] - mu -
             phi1 * (Data$y[2:(nn - 1)] - mu) -
             phi2 * (Data$y[1:(nn - 2)] - mu)
    LL <- sum(dnorm(resid, 0, sigma, log = TRUE))
    LP <- LL + dnorm(mu, 0, 100, log = TRUE) +
          dnorm(phi1, 0, 1, log = TRUE) +
          dnorm(phi2, 0, 1, log = TRUE) +
          dhalfcauchy(sigma, 10, log = TRUE)
    list(LP = LP, Dev = -2 * LL, Monitor = LP,
         yhat = rnorm(nn, mu, sigma), parm = parm)
}
Data_AR2 <- list(y = y, n = n,
                 parm.names = c("mu", "phi1", "phi2", "sigma"),
                 mon.names = "LP")
```

**Model 4: AR(3).** A further overparameterized extension with three
lags.

``` r

Model_AR3 <- function(parm, Data) {
    mu   <- parm[1]
    phi1 <- interval(parm[2], -0.999, 0.999); parm[2] <- phi1
    phi2 <- interval(parm[3], -0.999, 0.999); parm[3] <- phi2
    phi3 <- interval(parm[4], -0.999, 0.999); parm[4] <- phi3
    sigma <- interval(parm[5], 1e-100, Inf);   parm[5] <- sigma
    nn <- Data$n
    resid <- Data$y[4:nn] - mu -
             phi1 * (Data$y[3:(nn - 1)] - mu) -
             phi2 * (Data$y[2:(nn - 2)] - mu) -
             phi3 * (Data$y[1:(nn - 3)] - mu)
    LL <- sum(dnorm(resid, 0, sigma, log = TRUE))
    LP <- LL + dnorm(mu, 0, 100, log = TRUE) +
          dnorm(phi1, 0, 1, log = TRUE) +
          dnorm(phi2, 0, 1, log = TRUE) +
          dnorm(phi3, 0, 1, log = TRUE) +
          dhalfcauchy(sigma, 10, log = TRUE)
    list(LP = LP, Dev = -2 * LL, Monitor = LP,
         yhat = rnorm(nn, mu, sigma), parm = parm)
}
Data_AR3 <- list(y = y, n = n,
                 parm.names = c("mu", "phi1", "phi2", "phi3", "sigma"),
                 mon.names = "LP")
```

**Model 5: MA(1).** A moving-average model that captures serial
dependence through lagged innovations rather than lagged observations.

``` r

Model_MA1 <- function(parm, Data) {
    mu    <- parm[1]
    theta <- interval(parm[2], -0.999, 0.999)
    parm[2] <- theta
    sigma <- interval(parm[3], 1e-100, Inf)
    parm[3] <- sigma
    nn <- Data$n
    ## Compute residuals recursively
    eps <- numeric(nn)
    eps[1] <- Data$y[1] - mu
    for (t in 2:nn) {
        eps[t] <- Data$y[t] - mu - theta * eps[t - 1]
    }
    LL <- sum(dnorm(eps, 0, sigma, log = TRUE))
    LP <- LL + dnorm(mu, 0, 100, log = TRUE) +
          dnorm(theta, 0, 1, log = TRUE) +
          dhalfcauchy(sigma, 10, log = TRUE)
    list(LP = LP, Dev = -2 * LL, Monitor = LP,
         yhat = rnorm(nn, mu, sigma), parm = parm)
}
Data_MA1 <- list(y = y, n = n,
                 parm.names = c("mu", "theta1", "sigma"),
                 mon.names = "LP")
```

**Model 6: ARMA(1,1).** Combines one autoregressive and one
moving-average term. This model is more flexible than either AR(1) or
MA(1) alone but overparameterized for the true AR(1) process.

``` r

Model_ARMA11 <- function(parm, Data) {
    mu    <- parm[1]
    phi1  <- interval(parm[2], -0.999, 0.999); parm[2] <- phi1
    theta <- interval(parm[3], -0.999, 0.999); parm[3] <- theta
    sigma <- interval(parm[4], 1e-100, Inf);   parm[4] <- sigma
    nn <- Data$n
    eps <- numeric(nn)
    eps[1] <- Data$y[1] - mu
    for (t in 2:nn) {
        eps[t] <- Data$y[t] - mu - phi1 * (Data$y[t - 1] - mu) -
                  theta * eps[t - 1]
    }
    LL <- sum(dnorm(eps, 0, sigma, log = TRUE))
    LP <- LL + dnorm(mu, 0, 100, log = TRUE) +
          dnorm(phi1, 0, 1, log = TRUE) +
          dnorm(theta, 0, 1, log = TRUE) +
          dhalfcauchy(sigma, 10, log = TRUE)
    list(LP = LP, Dev = -2 * LL, Monitor = LP,
         yhat = rnorm(nn, mu, sigma), parm = parm)
}
Data_ARMA11 <- list(y = y, n = n,
                    parm.names = c("mu", "phi1", "theta1", "sigma"),
                    mon.names = "LP")
```

**Model 7: Random walk with drift.** A unit-root process \\y_t =
\delta + y\_{t-1} + \varepsilon_t\\. This is non-stationary and
fundamentally different from the true AR(1), which is stationary.
Including it tests whether the method can discriminate between
stationary and non-stationary dynamics.

``` r

Model_RW <- function(parm, Data) {
    delta <- parm[1]
    sigma <- interval(parm[2], 1e-100, Inf)
    parm[2] <- sigma
    ## Random walk: y_t = delta + y_{t-1} + eps_t
    resid <- diff(Data$y) - delta
    LL <- sum(dnorm(resid, 0, sigma, log = TRUE))
    LP <- LL + dnorm(delta, 0, 10, log = TRUE) +
          dhalfcauchy(sigma, 10, log = TRUE)
    list(LP = LP, Dev = -2 * LL, Monitor = LP,
         yhat = rnorm(Data$n, delta, sigma), parm = parm)
}
Data_RW <- list(y = y, n = n,
                parm.names = c("delta", "sigma"), mon.names = "LP")
```

### Running the product space analysis

We assemble all seven models and their data lists into the
[`ProductSpace()`](https://robustecologies.github.io/lucifer/reference/ProductSpace.md)
call. With seven models and a total of 22 parameters plus the model
indicator, the supermodel has dimension 23. The auto-FC Gibbs sampler
handles this by slice-sampling each continuous parameter and
exact-enumerating the discrete model indicator.

``` r

Models <- list(
    WN      = Model_WN,
    AR1     = Model_AR1,
    AR2     = Model_AR2,
    AR3     = Model_AR3,
    MA1     = Model_MA1,
    ARMA11  = Model_ARMA11,
    RW      = Model_RW
)

Data_list <- list(
    WN      = Data_WN,
    AR1     = Data_AR1,
    AR2     = Data_AR2,
    AR3     = Data_AR3,
    MA1     = Data_MA1,
    ARMA11  = Data_ARMA11,
    RW      = Data_RW
)

ps <- ProductSpace(
    Models = Models,
    Data   = Data_list,
    Pilot.Iterations = 5000,
    Pilot.Algorithm  = "CHARM",
    Iterations       = 50000,
    Thinning         = 10,
    Chains           = 3,
    CPUs             = 3,
    Bisection        = TRUE,
    Bisection.Tolerance  = 0.10,
    Bisection.Max.Rounds = 10,
    Bisection.Iterations = 5000
)
```

### Results

#### Print and summary

The [`print()`](https://rdrr.io/r/base/print.html) method provides a
compact one-screen summary. The corrected posterior model probabilities
account for any non-uniform prior model probabilities introduced by
bisection calibration. The Bayes factors are reported relative to the
best model, with Jeffreys scale interpretations. The model space
complexity metrics give a scalar summary of how concentrated the
posterior model distribution is.

``` r

print(ps)
```

The [`summary()`](https://rdrr.io/r/base/summary.html) method provides
the full Bayes factor matrix, transition matrix, per-model posterior
parameter summaries, and bisection calibration history.

``` r

summary(ps)
```

#### Posterior model space complexity

The complexity metrics are accessible as `ps$Complexity`. For this
example, where the AR(1) model should be strongly preferred, we expect a
low normalized Shannon entropy (near 0, indicating clear preference), an
effective number of models close to 1, and low Simpson diversity.

``` r

cx <- ps$Complexity
cat(sprintf("Shannon entropy (normalized): %.3f\n", cx$Shannon.Normalized))
cat(sprintf("Effective number of models:   %.2f / %d\n", cx$Effective.K, ps$K))
cat(sprintf("Simpson diversity:            %.3f\n", cx$Simpson.Diversity))
cat(sprintf("Gini effective K:             %.2f\n", cx$Gini.Effective.K))
```

#### Visualization

The default [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
call produces a combined five-panel layout with the model probabilities
bar chart, transition matrix heatmap, Bayes factor forest plot,
per-model posterior densities, and the full-width model indicator trace.

``` r

plot(ps)
```

Individual plot types can be accessed separately:

``` r

plot(ps, type = "probabilities")
```

The bar chart annotates each bar with the posterior probability,
highlights the preferred model in full saturation, and shows the Shannon
entropy and effective number of models in the subtitle. The dashed
horizontal line marks the uniform baseline \\1/K\\.

``` r

plot(ps, type = "trace")
```

The trace plot colors each iteration by the active model, making
switching patterns immediately visible. Frequent, irregular color
changes indicate good mixing.

``` r

plot(ps, type = "bayes_factors")
```

The Bayes factor forest plot shows log Bayes factors relative to the
best model. Background shading encodes the Jeffreys scale: red for
evidence against, yellow for anecdotal, teal for substantial, and green
for strong to decisive evidence. The caption at the bottom provides a
text key.

``` r

plot(ps, type = "transition")
```

The transition matrix shows how often the chain moves between each pair
of models. High off-diagonal values indicate that the sampler can freely
switch between models.

``` r

plot(ps, type = "posteriors")
```

The density plots show parameter posteriors conditional on each model
being active. For the AR(1) model, the posterior on \\\phi_1\\ should
concentrate near the true value 0.7.

``` r

plot(ps, type = "complexity")
```

The complexity bullet chart summarizes four information-theoretic
metrics of the posterior model probability distribution, each normalized
to \\\[0, 1\]\\ so they share a common axis. A value near zero indicates
that a single model dominates the posterior (decisive selection), while
a value near one indicates that posterior mass is spread uniformly
across all \\K\\ candidates (maximum ambiguity). Three background zones,
decisive (green, \\\[0, 1/3\]\\), moderate (yellow, \\(1/3, 2/3\]\\),
and ambiguous (red, \\(2/3, 1\]\\), provide immediate visual
classification. The raw metric value and its theoretical range \\\[\min,
\max\]\\ are annotated to the right of each bar; for instance, the
effective number of models ranges from 1 (all mass on a single model) to
\\K\\ (uniform), while the Simpson diversity ranges from 0 to
\\(K-1)/K\\. The subtitle reports an overall regime assessment derived
from the average of the four normalized scores.

#### Extracting per-model posteriors

The per-model posterior samples are stored as matrices in
`ps$Per.Model.Posteriors`. These are the MCMC draws for each model’s
parameters conditional on that model being the active model, and can be
used for any further analysis (credible intervals, posterior predictive
checks, etc.).

``` r

## AR(1) posterior summary
ps$Per.Model.Summaries$AR1

## Number of samples from each model
sapply(ps$Per.Model.Posteriors, nrow)
```

  

## Benchmark: model recovery under increasing noise

A natural question is how well the product space method recovers the
true generative model as the signal-to-noise ratio decreases. We
simulate AR(1) data at five noise levels \\\sigma \in \\0.5, 1.0, 2.0,
3.0, 5.0\\\\ while keeping the deterministic structure fixed (\\\alpha =
2\\, \\\phi_1 = 0.7\\). To isolate the effect of noise amplitude from
random-seed variability, we generate a single standard-normal innovation
sequence and rescale it by \\\sigma\\ at each level, so the same
underlying random pattern drives all five series. For each noise level,
we run
[`ProductSpace()`](https://robustecologies.github.io/lucifer/reference/ProductSpace.md)
with the same seven candidate models and record the posterior
probability assigned to the true AR(1) model, the effective number of
models, and the Bayes factor of AR(1) versus the best alternative.

``` r

library(lucifer)

sigma_levels <- c(0.5, 1.0, 2.0, 3.0, 5.0)
n <- 200
alpha_true <- 2.0
phi1_true  <- 0.7

## Common innovation sequence: same standard-normal draws scaled by sigma,
## so only the noise amplitude varies across levels (no seed-to-seed noise).
set.seed(42)
eps_common <- rnorm(n)

results <- data.frame(
    sigma     = numeric(0),
    p_AR1     = numeric(0),
    eff_K     = numeric(0),
    BF_AR1    = numeric(0),
    H_norm    = numeric(0),
    best      = character(0),
    stringsAsFactors = FALSE
)
```

``` r

for (sig in sigma_levels) {

    ## Simulate AR(1) using common innovations scaled by sigma
    y <- numeric(n)
    y[1] <- alpha_true + eps_common[1] * sig / sqrt(1 - phi1_true^2)
    for (t in 2:n) {
        y[t] <- alpha_true + phi1_true * (y[t - 1] - alpha_true) +
                eps_common[t] * sig
    }

    ## Rebuild Data lists with new y
    mk_data <- function(pnames) list(y = y, n = n,
        parm.names = pnames, mon.names = "LP")

    DL <- list(
        WN     = mk_data(c("mu", "sigma")),
        AR1    = mk_data(c("mu", "phi1", "sigma")),
        AR2    = mk_data(c("mu", "phi1", "phi2", "sigma")),
        AR3    = mk_data(c("mu", "phi1", "phi2", "phi3", "sigma")),
        MA1    = mk_data(c("mu", "theta1", "sigma")),
        ARMA11 = mk_data(c("mu", "phi1", "theta1", "sigma")),
        RW     = mk_data(c("delta", "sigma"))
    )

    ps_i <- ProductSpace(
        Models           = Models,
        Data             = DL,
        Pilot.Iterations = 5000,
        Iterations       = 50000,
        Thinning         = 10,
        Chains           = 3,
        CPUs             = 3,
        Bisection        = TRUE,
        Bisection.Max.Rounds = 8,
        Bisection.Iterations = 5000,
        Status           = 100000,
        Verbose          = FALSE
    )

    p_ar1 <- ps_i$Posterior.Model.Probabilities["AR1"]
    best_model <- names(which.max(ps_i$Posterior.Model.Probabilities))

    ## BF of AR1 vs. best alternative
    if (best_model == "AR1") {
        others <- setdiff(names(ps_i$Posterior.Model.Probabilities), "AR1")
        bf_best_alt <- min(ps_i$Bayes.Factors["AR1", others])
    } else {
        bf_best_alt <- ps_i$Bayes.Factors["AR1", best_model]
    }

    results <- rbind(results, data.frame(
        sigma  = sig,
        p_AR1  = p_ar1,
        eff_K  = ps_i$Complexity$Effective.K,
        BF_AR1 = bf_best_alt,
        H_norm = ps_i$Complexity$Shannon.Normalized,
        best   = best_model,
        stringsAsFactors = FALSE
    ))

    cat(sprintf("sigma = %.1f: p(AR1) = %.3f, eff.K = %.2f, best = %s\n",
                sig, p_ar1, ps_i$Complexity$Effective.K, best_model))
}
```

``` r

print(results)
```

### Visualizing the benchmark

``` r

library(ggplot2)

p1 <- ggplot(results, aes(x = sigma, y = p_AR1)) +
    geom_line(color = "#4e79a7", linewidth = 1) +
    geom_point(color = "#4e79a7", size = 3) +
    geom_hline(yintercept = 1/7, linetype = "dashed", color = "grey55") +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = expression(sigma ~ "(noise level)"),
         y = "p(AR1 | y)",
         title = "Posterior probability of the true model",
         subtitle = "Dashed line: uniform baseline (1/7)") 

p2 <- ggplot(results, aes(x = sigma, y = eff_K)) +
    geom_line(color = "#e15759", linewidth = 1) +
    geom_point(color = "#e15759", size = 3) +
    geom_hline(yintercept = 1, linetype = "dotted", color = "grey55") +
    geom_hline(yintercept = 7, linetype = "dotted", color = "grey55") +
    scale_y_continuous(limits = c(0.8, 7.2)) +
    labs(x = expression(sigma ~ "(noise level)"),
         y = expression(K[eff]),
         title = "Effective number of models",
         subtitle = "Dotted lines: 1 (full certainty) and 7 (maximum uncertainty)")

gridExtra::grid.arrange(p1, p2, ncol = 2)
```

Because every series is driven by the same innovation sequence, the
trend across \\\sigma\\ is clean and monotonic. At low noise (\\\sigma =
0.5\\) the autoregressive structure is unmistakable, and the posterior
probability of the true AR(1) model is near 1. As noise increases, the
AR(1) signal is progressively buried, the posterior probability of AR(1)
decreases, and the effective number of models rises, reflecting growing
uncertainty over which model generated the data. At very high noise
(\\\sigma = 5\\) the series resembles white noise and simpler models
(WN) may become competitive or preferred, demonstrating the Bayesian
Occam’s razor in action.

  

## Diagnostics and troubleshooting

### Model indicator mixing

The most important diagnostic is the trace plot of the model indicator
\\M\\. A well-mixing chain switches frequently between models, producing
a trace that looks like random noise across the \\K\\ levels. A poorly
mixing chain gets “stuck” at one model value for long stretches,
indicating that the pseudopriors are inadequate or the model
probabilities need calibration.

The switching rate, reported by
[`print()`](https://rdrr.io/r/base/print.html) and accessible as
`ps$Switching.Rate`, quantifies the proportion of consecutive iterations
where \\M\\ changed. Values above 0.2 indicate good mixing; values below
0.05 suggest problems. The transition matrix
(`plot(ps, type = "transition")`) provides more granular information:
ideally, all off-diagonal entries should be substantially positive.

### When bisection is needed

If the raw posterior model probabilities are very unbalanced (e.g., 0.99
vs. 0.01), the chain spends almost all its time in the preferred model,
leaving very few samples to estimate the Bayes factor. The bisection
algorithm adjusts the prior model probabilities to force approximately
equal posterior activation, then corrects the final probabilities
analytically. This is enabled by default (`Bisection = TRUE`) and is
generally recommended. Setting `Bisection = FALSE` is appropriate when
the models are expected to have similar posterior probabilities, or when
computation time is a concern and uncalibrated estimates suffice.

When the evidence is overwhelming (e.g., comparing a white noise model
against AR(1) for strongly autocorrelated data), the bisection algorithm
cannot force the chain to visit the disfavored model regardless of how
extreme the prior model probability is set. In this case,
[`ProductSpace()`](https://robustecologies.github.io/lucifer/reference/ProductSpace.md)
detects the failure (no model switching for 3 consecutive rounds),
terminates bisection early, and issues a warning. The resulting Bayes
factors should be interpreted as lower bounds on the true evidence.

### Pseudoprior quality

The automatic pseudoprior estimation relies on pilot MCMC fits. If the
pilot fits have poor convergence (too few iterations, bad starting
values), the pseudopriors will be inaccurate and mixing will suffer.
Increasing `Pilot.Iterations` or using a more robust `Pilot.Algorithm`
(e.g., `"HARM"` or `"twalk"`) can help. User-supplied pseudopriors via
the `Pseudopriors` argument allow full control when the automatic
estimates are unsatisfactory.

### Number of models and scalability

The product space method scales linearly in \\K\\ for the discrete
enumeration step. The hybrid sampler (`Production.Algorithm = "RAM"`,
default) eliminates the dependence on \\\sum d_k\\ in the per-iteration
cost: inactive parameters are direct-sampled regardless of their
dimension, and the active model gets a single block proposal. This makes
the method practical for problems with many models or high-dimensional
parameter spaces that were previously bottlenecked by univariate slice
sampling. For very large model sets (\\K \> 20\\), the enumeration step
itself can become expensive; in such cases, consider pre-screening with
information criteria
([`WAIC()`](https://robustecologies.github.io/lucifer/reference/WAIC.md),
[`LOO()`](https://robustecologies.github.io/lucifer/reference/LOO.md))
and applying
[`ProductSpace()`](https://robustecologies.github.io/lucifer/reference/ProductSpace.md)
to a shortlist.

The original Gibbs sampler is available via
`Production.Algorithm = "Gibbs"` and may be preferred for models with
very low-dimensional parameter spaces where univariate slice sampling
converges quickly.

  

## References

**\[1\]** Carlin, B.P. and Chib, S. (1995). “Bayesian model choice via
Markov chain Monte Carlo methods.” *Journal of the Royal Statistical
Society: Series B*, 57(3), p. 473-484.

**\[2\]** Lodewyckx, T., Kim, W., Lee, M.D., Tuerlinckx, F., Kuppens,
P., and Wagenmakers, E.-J. (2011). “A tutorial on Bayes factor
estimation with the product space method.” *Journal of Mathematical
Psychology*, 55(5), p. 331-347.
[doi:10.1016/j.jmp.2011.06.001](https://doi.org/10.1016/j.jmp.2011.06.001)

**\[3\]** Tenan, S., O’Hara, R.B., Hendriks, I., and Tavecchia, G.
(2014). “Bayesian model selection: the steepest mountain to climb.”
*Ecological Modelling*, 283, p. 62-69.
[doi:10.1016/j.ecolmodel.2014.03.017](https://doi.org/10.1016/j.ecolmodel.2014.03.017)

**\[4\]** Neal, R.M. (2003). “Slice sampling.” *Annals of Statistics*,
31(3), p. 705-767.
[doi:10.1214/aos/1056562461](https://doi.org/10.1214/aos/1056562461)

**\[5\]** Green, P.J. (1995). “Reversible jump Markov chain Monte Carlo
computation and Bayesian model determination.” *Biometrika*, 82(4),
p. 711-732.
[doi:10.1093/biomet/82.4.711](https://doi.org/10.1093/biomet/82.4.711)

**\[6\]** Vihola, M. (2012). “Robust adaptive Metropolis algorithm with
coerced acceptance rate.” *Statistics and Computing*, 22(5),
p. 997-1008.
[doi:10.1007/s11222-011-9269-5](https://doi.org/10.1007/s11222-011-9269-5)

**\[7\]** Godsill, S.J. (2001). “On the relationship between Markov
chain Monte Carlo methods for model uncertainty.” *Journal of
Computational and Graphical Statistics*, 10(2), p. 230-248.
[doi:10.1198/10618600152627924](https://doi.org/10.1198/10618600152627924)

**\[8\]** Dellaportas, P., Forster, J.J., and Ntzoufras, I. (2002). “On
Bayesian model and variable selection using MCMC.” *Statistics and
Computing*, 12(1), p. 27-36.
[doi:10.1023/A:1013164120801](https://doi.org/10.1023/A:1013164120801)

**\[9\]** Meng, X.-L. and Wong, W.H. (1996). “Simulating ratios of
normalizing constants via a simple identity: a theoretical exploration.”
*Statistica Sinica*, 6(4), p. 831-860.

**\[10\]** Han, C. and Carlin, B.P. (2001). “Markov chain Monte Carlo
methods for computing Bayes factors: a comparative review.” *Journal of
the American Statistical Association*, 96(455), p. 1122-1132.
[doi:10.1198/016214501753208780](https://doi.org/10.1198/016214501753208780)

**\[11\]** Kass, R.E. and Raftery, A.E. (1995). “Bayes factors.”
*Journal of the American Statistical Association*, 90(430), p. 773-795.
[doi:10.1080/01621459.1995.10476572](https://doi.org/10.1080/01621459.1995.10476572)

**\[12\]** Jeffreys, H. (1961). *Theory of probability* (3rd ed.).
Oxford University Press.
