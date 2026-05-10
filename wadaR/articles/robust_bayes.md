# Robust Bayesian sensitivity analysis

## Motivation: why a single posterior is never enough

Bayesian inference produces a posterior distribution \\p(\theta \| y)\\
that synthesizes prior knowledge and observed data. In practice, every
posterior is conditioned on at least three layers of modeling choices,
each of which can introduce fragility that standard posterior summaries
leave invisible.

The first layer is the prior \\\pi(\theta)\\. When data are sparse
relative to parameter dimensionality, the prior can dominate the
posterior in ways that are difficult to detect by inspecting summaries
alone. A posterior mean that appears precise may shift substantially
under modest prior perturbation, revealing that the apparent precision
is an artifact of the prior rather than a consequence of the data. The
problem is compounded in hierarchical models, where hyperpriors
propagate influence through multiple levels.

The second layer is the likelihood \\L(y\|\theta)\\. Misspecification of
the observation model, whether through incorrect distributional
assumptions, omitted covariates, or structural simplifications, can
produce posteriors that are internally consistent but externally
invalid. Two analysts fitting different but equally plausible likelihood
specifications to the same data may reach qualitatively different
conclusions, and neither posterior alone reveals the extent of this
disagreement.

The third layer is the data itself. Individual observations can exert
disproportionate leverage on the posterior, particularly in small
samples or when the model lacks robust error distributions. An outlier
that shifts the posterior mean by two standard deviations effectively
determines the conclusion, yet this influence is invisible in the
marginal posterior.

[`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md)
provides a unified interface for diagnosing all three sources of
sensitivity through five complementary analyses. Rather than requiring
separate tools for each diagnostic, it accepts any lucifer fit object
(or a family of fits) and returns a single `robust_bayes` S3 object with
print, summary, and plot methods. The function exploits a key property
of the lucifer Model interface: every Model evaluation returns both `LP`
(log-posterior) and `Dev` (\\-2 \times\\ log-likelihood), enabling
automatic decomposition into log-prior and log-likelihood without
requiring the user to write separate functions for these components.

## The LP/Dev decomposition

At the core of several
[`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md)
modules lies a decomposition that the lucifer Model interface provides
for free. Every lucifer Model function returns a list with `LP` (the
log-posterior, equal to log-prior plus log-likelihood) and `Dev` (the
deviance, equal to \\-2\\ times the log-likelihood). From these two
quantities, we recover the components:

\\\ell(\theta) = -\frac{\text{Dev}(\theta)}{2}\\

\\\log \pi(\theta) = \text{LP}(\theta) - \ell(\theta) =
\text{LP}(\theta) + \frac{\text{Dev}(\theta)}{2}\\

This decomposition is exact regardless of the model’s internal
structure. It does not require the user to code the prior and likelihood
as separate functions, and it works identically for simple conjugate
models and for complex hierarchical specifications where the prior
itself involves integration over latent variables. The computational
cost is one Model evaluation per posterior sample, which
[`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md)
performs once and caches for reuse across the power-scaling, conflict,
and influence modules.

## Mathematical foundations

### Module 1: power-scaling sensitivity

#### The power-scaling framework

Power-scaling sensitivity analysis, formalized by Kallioinen et
al. [\[1\]](#ref1), evaluates how the posterior changes when the prior
or likelihood is raised to a power \\\alpha \> 0\\. The power-scaled
posterior for prior perturbation is

\\\tilde{p}(\theta \| y, \alpha) \propto \pi(\theta)^\alpha \\
L(y\|\theta)\\

and for likelihood perturbation is

\\\tilde{p}(\theta \| y, \alpha) \propto \pi(\theta) \\
L(y\|\theta)^\alpha\\

When \\\alpha = 1\\ both expressions reduce to the original posterior.
Values \\\alpha \< 1\\ weaken the component (making the prior flatter or
the likelihood less informative), while \\\alpha \> 1\\ strengthen it
(making the prior more concentrated or the likelihood more dominant).
The key insight is that the power-scaled posterior can be approximated
without refitting by importance reweighting from the original posterior
samples.

#### Importance weights and PSIS stabilization

For a set of posterior samples \\\\\theta_s\\\_{s=1}^S\\ drawn from
\\p(\theta\|y)\\, the importance weights for approximating the
prior-perturbed posterior are

\\\log w_s^{(\alpha)} = (\alpha - 1) \cdot \log \pi(\theta_s)\\

and for the likelihood-perturbed posterior:

\\\log w_s^{(\alpha)} = (\alpha - 1) \cdot \ell(\theta_s)\\

These raw importance weights can be highly variable, particularly when
\\\alpha\\ is far from 1, because a few posterior samples may have
extreme prior or likelihood values. Pareto smoothed importance sampling
(PSIS) [\[4\]](#ref4) addresses this by fitting a generalized Pareto
distribution (GPD) to the upper tail of the log-weight distribution and
replacing the most extreme weights with smoothed quantiles from the
fitted GPD. The Zhang and Stephens (2009) estimator with the weakly
informative prior adjustment of Vehtari et al. (2024) is used for GPD
fitting.

The Pareto shape parameter \\\hat{k}\\ from the GPD fit serves as a
diagnostic for the reliability of the importance sampling approximation.
When \\\hat{k} \< 0.5\\, the approximation is reliable; when \\0.5 \leq
\hat{k} \< 0.7\\, it is acceptable; when \\\hat{k} \geq 0.7\\, the
effective sample size has dropped too far for the approximation to be
trusted, and the sensitivity estimate at that \\\alpha\\ value should be
interpreted cautiously.

#### Jensen-Shannon divergence as a sensitivity measure

For each value of \\\alpha\\ in the grid,
[`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md)
quantifies the distance between the original and power-scaled posteriors
using the Jensen-Shannon (JS) divergence on marginal posteriors. For
parameter \\d\\, the JS divergence between the original marginal \\p_d\\
and the reweighted marginal \\\tilde{p}\_d^{(\alpha)}\\ is

\\\text{JS}(p_d \\ \tilde{p}\_d^{(\alpha)}) = \frac{1}{2}\text{KL}(p_d
\\ m_d) + \frac{1}{2}\text{KL}(\tilde{p}\_d^{(\alpha)} \\ m_d)\\

where \\m_d = \frac{1}{2}(p_d + \tilde{p}\_d^{(\alpha)})\\ is the
mixture distribution. The JS divergence is symmetric, bounded between 0
and \\\log 2\\, and well-defined even when the supports of the two
distributions do not overlap perfectly. It is estimated via histograms
with \\B = 200\\ bins, where the original posterior uses uniform weights
\\1/S\\ and the reweighted posterior uses the PSIS-smoothed weights.

#### Cumulative JS distance and scalar sensitivity scores

A single JS divergence at one value of \\\alpha\\ provides limited
information. What matters is how rapidly the posterior changes as
\\\alpha\\ moves away from 1. The cumulative Jensen-Shannon distance
integrates the instantaneous divergence across the perturbation path:

\\D\_{\text{CJS}}^{(\alpha)} = \int_1^{\alpha}
\sqrt{\text{JS}(\tilde{p}^{(a)}, p)} \\ da\\

This integral is approximated by trapezoidal quadrature over the
\\\alpha\\ grid. The scalar sensitivity score for each parameter is the
absolute derivative of \\D\_{\text{CJS}}\\ at \\\alpha = 1\\, estimated
by central differences:

\\s_d = \left\|\frac{dD\_{\text{CJS},d}}{d\alpha}\right\|\_{\alpha=1}
\approx \frac{\|D\_{\text{CJS},d}^{(\alpha\_+)} -
D\_{\text{CJS},d}^{(\alpha\_-)}\|}{\alpha\_+ - \alpha\_-}\\

where \\\alpha\_+\\ and \\\alpha\_-\\ are the grid points immediately
above and below 1. A large sensitivity score indicates that even a small
perturbation to the prior or likelihood materially changes the marginal
posterior for that parameter.

#### Diagnostic classification

Following Kallioinen et al. [\[1\]](#ref1),
[`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md)
classifies each parameter into one of four diagnostic categories based
on whether its prior sensitivity and likelihood sensitivity exceed a
user-specified threshold (default 0.05):

| Prior sensitivity | Likelihood sensitivity | Diagnosis |
|:---|:---|:---|
| Low | Low | Posterior robust |
| High | Low | Prior-dominated (weak data) |
| Low | High | Likelihood-dominated (strong data) |
| High | High | Prior-data conflict |

Diagnostic classification from power-scaling sensitivity. {.table}

A “robust” posterior is one whose conclusions survive reasonable
perturbations to both the prior and the likelihood; this is the ideal
case. “Prior-dominated” parameters indicate that the data carry
insufficient information to override the prior specification; the
analyst should consider whether the prior is scientifically justified or
merely a default. “Likelihood-dominated” parameters are well-informed by
the data and insensitive to prior choice; this is typical for
well-identified parameters with moderate to large samples. “Conflict”
parameters are sensitive to perturbations from both directions, often
indicating that the prior and the data pull in different directions.

### Module 2: cross-model posterior comparison

#### The need for divergence measures

When multiple model specifications are scientifically plausible, an
analyst should understand the extent to which posterior inferences
depend on the choice of model. Two models may yield qualitatively
similar posteriors for some parameters but diverge sharply for others.
[`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md)
computes pairwise divergence matrices across all models, quantifying
both marginal (per-parameter) and joint (multivariate) disagreement.

#### Jensen-Shannon divergence

The JS divergence between two marginal posteriors from models \\k\\ and
\\k'\\ for parameter \\d\\ is estimated via histograms on the union
support range:

\\\text{JS}(P\_{kd} \\ P\_{k'd}) = \frac{1}{2}\text{KL}(P\_{kd} \\
M_d) + \frac{1}{2}\text{KL}(P\_{k'd} \\ M_d)\\

where \\P\_{kd}\\ and \\P\_{k'd}\\ are empirical distributions from the
posterior samples and \\M_d = \frac{1}{2}(P\_{kd} + P\_{k'd})\\. Both
distributions are histogrammed into \\B = 200\\ bins on their joint
range, and the KL divergences are computed as sums over bin
probabilities. The JS divergence is zero when the posteriors are
identical and maximal (\\\log 2 \approx 0.693\\) when they have disjoint
supports.

#### Wasserstein-2 distance

The Wasserstein-2 distance provides a metric that accounts for the
geometry of the parameter space, unlike the JS divergence which depends
only on density overlap:

\\W_2(P, Q) = \left(\int_0^1 (F_P^{-1}(t) - F_Q^{-1}(t))^2 \\
dt\right)^{1/2}\\

For empirical distributions this has an exact closed-form solution: sort
both samples into order statistics \\x\_{(1)} \leq \cdots \leq
x\_{(S)}\\ and \\y\_{(1)} \leq \cdots \leq y\_{(S)}\\, then \\W_2 =
\left(\frac{1}{S}\sum\_{s=1}^S (x\_{(s)} - y\_{(s)})^2\right)^{1/2}\\.
When sample sizes differ,
[`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md)
linearly interpolates the smaller sample’s quantile function onto a
common grid of \\\max(S_k, S\_{k'})\\ quantile levels. The Wasserstein
distance has the advantage of being interpretable in the units of the
parameter; a Wasserstein distance of 0.5 for a regression coefficient
means the posteriors disagree by half a unit on average.

#### Overlap coefficient

The overlap coefficient provides the most directly interpretable measure
of posterior agreement:

\\\text{OVL}(P, Q) = \int \min(f_P(x), f_Q(x)) \\ dx\\

This quantity equals 1 when the distributions are identical and
approaches 0 when they have disjoint supports. It is estimated from
histogram bins as the sum of bin-wise minima of the normalized
histograms. An overlap of 0.85 means that 85% of the probability mass is
shared between the two posteriors; the remaining 15% represents the
region where they disagree.

#### Energy distance for multivariate comparison

Marginal divergences can miss important joint structure. Two models
might agree on the marginal distributions of each parameter individually
while disagreeing on the correlation structure. The energy distance
captures multivariate disagreement without requiring density estimation:

\\E(P, Q) = 2\mathbb{E}\\X - Y\\ - \mathbb{E}\\X - X'\\ -
\mathbb{E}\\Y - Y'\\\\

where \\X, X'\\ are independent draws from \\P\\ and \\Y, Y'\\ are
independent draws from \\Q\\, and \\\\\cdot\\\\ denotes the Euclidean
norm. This is estimated by random subsampling of 1000 pairs from each
posterior for computational tractability. The energy distance is zero if
and only if \\P = Q\\, and it metrizes convergence in distribution,
making it a proper metric on the space of probability measures.

#### Parameter sensitivity ranking

For each parameter \\d\\, the maximum divergence across all
\\\binom{K}{2}\\ model pairs identifies which parameters are most
sensitive to model choice. This ranking helps focus attention on the
parameters whose inferences are least robust, regardless of which
specific pair of models drives the disagreement.

### Module 3: prior-data conflict diagnostics

#### The concept of prior-data conflict

Prior-data conflict occurs when the prior distribution and the
likelihood provide contradictory information about the parameters. This
situation is qualitatively different from having an “uninformative” or
“diffuse” prior; a diffuse prior simply adds little information, while a
conflicting prior actively pulls the posterior away from where the data
suggest it should be. The consequences can be severe: the posterior may
be biased toward the prior, have artificially narrow credible intervals,
or exhibit poor predictive performance.

Two complementary diagnostics are implemented: the Box prior predictive
p-value, which tests for global prior-data conflict, and the
Evans-Moshonov calibration, which tests each parameter individually.

#### Box (1980) prior predictive p-value

Box’s approach [\[2\]](#ref2) uses the prior predictive distribution to
assess whether the observed data are “surprising” relative to what the
prior expects. The algorithm proceeds as follows:

1.  Draw \\R = 1000\\ parameter vectors from the prior distribution via
    `Data$PGF(Data)`.
2.  For each prior draw \\\theta^{(r)}\\, evaluate the Model function
    and extract the deviance \\T^{(r)} = -2\ell(\theta^{(r)})\\.
3.  Compute the observed deviance \\T\_{\text{obs}}\\ at the posterior
    mean \\\hat\theta = \bar\theta\_{\text{post}}\\.
4.  The p-value is the proportion of prior-predictive deviances that
    exceed the observed deviance: \\p\_{\text{Box}} =
    \frac{1}{R}\sum\_{r=1}^R \mathbf{1}\[T^{(r)} \geq
    T\_{\text{obs}}\]\\.

Interpretation requires care. A p-value near 0 indicates that the
observed data fit better than nearly all data generated under the prior,
which can happen when the prior places mass far from the data-generating
region, making typical prior-predictive datasets very poorly-fitting
relative to the actual data. A p-value near 1 indicates the opposite
situation. Both extremes (below 0.05 or above 0.95) suggest prior-data
conflict; the former is more common in practice.

The Box p-value has the virtue of simplicity and generality: it works
for any model that can be evaluated at prior draws, requires no special
distributional assumptions, and provides a single scalar summary of
global conflict. Its limitation is that it conflates different sources
of conflict; a small p-value does not tell you which parameters or
aspects of the prior are problematic.

#### Evans and Moshonov (2006) per-parameter calibration

The Evans-Moshonov calibration [\[3\]](#ref3) addresses the limitation
of the Box p-value by providing parameter-level diagnostics. For each
parameter \\d\\:

1.  Compute the posterior mean \\\hat\theta_d\\ from the posterior
    samples.
2.  Estimate the marginal prior density at the posterior mean,
    \\\hat{f}\_{\pi,d}(\hat\theta_d)\\, using a kernel density estimate
    on the prior draws.
3.  Compute the proportion of prior draws with density less than or
    equal to the density at the posterior mean: \\C_d =
    \Pr\_\pi\[\hat{f}\_{\pi,d}(\theta_d) \leq
    \hat{f}\_{\pi,d}(\hat\theta_d)\]\\.

The calibration value \\C_d\\ measures how extreme the posterior mean is
relative to the prior. If \\C_d\\ is small (e.g., below 0.05), it means
the posterior has concentrated in a region of low prior density, the
tails of the prior distribution, indicating that the data are pulling
the posterior away from where the prior thinks the parameter should be.
This is the hallmark of prior-data conflict for that parameter.

[`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md)
classifies parameters as “none” (\\C_d \geq 0.05\\), “mild” (\\0.01 \leq
C_d \< 0.05\\), or “severe” (\\C_d \< 0.01\\) conflict. The
classification uses kernel density estimation with Silverman’s
rule-of-thumb bandwidth, providing a nonparametric assessment that does
not assume any particular prior form.

### Module 4: observation influence analysis

#### Why observation-level diagnostics matter

In any finite dataset, some observations contribute more to the
posterior than others. This is expected and not inherently problematic:
observations that are highly informative about the parameters should
carry more weight. The concern arises when a small number of
observations dominate the posterior to the extent that removing them
would substantially change the conclusions. Such influential
observations may be outliers, data entry errors, or genuine but rare
events that the model handles poorly.

#### PSIS-LOO importance weights

The leave-one-out (LOO) posterior for observation \\i\\ is

\\p(\theta \| y\_{-i}) \propto p(\theta \| y) \cdot \frac{1}{p(y_i \|
\theta)}\\

which can be approximated by reweighting the full posterior with
importance weights

\\\log w_s^{(-i)} = -\log p(y_i \| \theta_s)\\

These weights are then stabilized via PSIS, producing smoothed,
normalized weights and a Pareto-\\\hat{k}\\ diagnostic for each
observation. This requires the \\N \times S\\ matrix of pointwise
log-likelihoods, the same input used by
[`LOO()`](https://robustecologies.github.io/lucifer/reference/LOO.md)
and
[`WAIC()`](https://robustecologies.github.io/lucifer/reference/WAIC.md).
The user can supply this matrix directly via the `log_lik` argument, or
[`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md)
will attempt to reconstruct it from Model evaluations.

The Pareto-\\\hat{k}\\ values have a direct interpretation as measures
of observation leverage. When \\\hat{k}\_i \> 0.7\\, the importance
weights for observation \\i\\ are so variable that the LOO approximation
is unreliable, which itself indicates that observation \\i\\ has strong
influence on the posterior. Values above 1.0 indicate that the
observation is so influential that no finite-variance importance
sampling scheme can approximate its leave-one-out effect.

#### Influence on posterior moments

Beyond the Pareto-\\\hat{k}\\ flag,
[`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md)
quantifies the direction and magnitude of each observation’s influence.
For each observation \\i\\ and parameter \\d\\, the shift in posterior
mean is

\\\Delta\bar\theta_d^{(-i)} = \sum\_{s=1}^S w_s^{(-i)} \theta\_{s,d} -
\bar\theta_d\\

and the shift in posterior standard deviation is computed analogously
from the weighted variance. These quantities reveal not just which
observations are influential, but how they influence each parameter: an
outlier might pull the mean of one parameter while inflating the
variance of another.

The maximum absolute mean shift across parameters, \\\max_d
\|\Delta\bar\theta_d^{(-i)}\|\\, provides a single scalar summary of
each observation’s influence, useful for identifying the most impactful
data points in a single pass.

### Module 5: robust posterior aggregation

#### The M-open perspective

In the “M-open” view of model uncertainty [\[6\]](#ref6), the true
data-generating process is not assumed to be among the candidate models.
Under this perspective, model selection (choosing a single model)
discards valuable information, and model averaging that accounts for
model uncertainty can produce better calibrated predictions. Three
complementary weighting strategies are implemented, each reflecting a
different philosophical and practical stance on how models should be
combined.

#### Bayesian model averaging (BMA) weights

The classical BMA weights are proportional to the marginal likelihood of
each model:

\\w_k^{\text{BMA}} = \frac{p(y \| \mathcal{M}\_k)}{\sum\_{j=1}^K p(y \|
\mathcal{M}\_j)} = \frac{\exp(\text{LML}\_k)}{\sum\_{j=1}^K
\exp(\text{LML}\_j)}\\

where \\\text{LML}\_k\\ is the log marginal likelihood (available as
`$LML` in every lucifer fit object). BMA weights have a clear Bayesian
interpretation as posterior model probabilities under flat model priors,
and they are computationally cheap since the LML is already computed
during fitting. Their limitation, well-documented by Yao et
al. [\[5\]](#ref5), is that they optimize for identifying the true
model, not for predictive accuracy. In the M-open setting where no
candidate model is true, BMA tends to concentrate weight on a single
model, producing overconfident predictions that do not adequately hedge
across model uncertainty.

[`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md)
extracts the LML from each fit object and computes BMA weights in
log-space for numerical stability. Models without a valid LML (e.g.,
fits that did not converge) are excluded with a warning.

#### Stacking weights

Stacking [\[5\]](#ref5) directly optimizes the leave-one-out predictive
performance of the model mixture:

\\\hat{w} = \arg\max\_{w \in \mathcal{S}} \sum\_{i=1}^N
\log\left(\sum\_{k=1}^K w_k \\ \hat{p}(y_i \| y\_{-i},
\mathcal{M}\_k)\right)\\

where \\\mathcal{S} = \\w \in \mathbb{R}^K : w_k \geq 0, \sum_k w_k =
1\\\\ is the probability simplex and \\\hat{p}(y_i \| y\_{-i},
\mathcal{M}\_k)\\ is the LOO predictive density for observation \\i\\
under model \\k\\, obtained from PSIS-LOO.

The optimization is solved via
[`stats::constrOptim()`](https://rdrr.io/r/stats/constrOptim.html) with
BFGS updates, using analytical gradients for efficiency. The gradient of
the negative stacking criterion with respect to \\w_k\\ is

\\\frac{\partial}{\partial w_k}\left\[-\sum\_{i=1}^N \log\left(\sum_j
w_j \hat{p}\_{ij}\right)\right\] = -\sum\_{i=1}^N
\frac{\hat{p}\_{ik}}{\sum_j w_j \hat{p}\_{ij}}\\

Stacking weights have several advantages over BMA: they optimize for
prediction rather than model identification; they naturally handle the
M-open setting; and they tend to distribute weight more evenly across
models that contribute complementary predictive information, rather than
concentrating on a single “best” model.

#### Pseudo-BMA+ weights

Pseudo-BMA+ [\[5\]](#ref5) provides a Bayesian bootstrap correction that
accounts for the estimation uncertainty in the LOO-PSIS pointwise log
predictive densities. The algorithm is:

1.  For \\b = 1, \ldots, B\\ bootstrap replicates:
    1.  Draw \\\alpha^{(b)} \sim \text{Dirichlet}(1, \ldots, 1)\\ by
        normalizing \\K\\ independent \\\text{Gamma}(1,1)\\ draws.
    2.  Compute the bootstrap-weighted elpd for each model:
        \\\text{elpd}\_k^{(b)} = \sum\_{i=1}^N \alpha_i^{(b)} \cdot
        \text{elpd}\_{ik}\\.
    3.  Compute softmax weights for this replicate: \\w_k^{(b)} \propto
        \exp(\text{elpd}\_k^{(b)} - \max_j \text{elpd}\_j^{(b)})\\.
2.  Average across replicates: \\w_k^{\text{pBMA+}} =
    \frac{1}{B}\sum\_{b=1}^B w_k^{(b)}\\.

The Dirichlet resampling perturbs the relative importance of
observations, which propagates into uncertainty about the model elpd
rankings. Models that are consistently preferred across bootstrap
replicates receive high weight; models whose advantage depends on a few
observations receive lower weight. This correction is particularly
important when the elpd differences between models are small relative to
their estimation uncertainty.

#### Aggregated posterior

The final step constructs an aggregated posterior by mixture sampling.
For each of \\S\\ desired posterior draws, a model \\k\\ is selected
with probability \\w_k\\ (using the stacking weights by default, or BMA
if stacking is unavailable), and then a sample is drawn uniformly from
that model’s posterior. The result is an \\S \times D\\ matrix of
posterior samples from the model-averaged distribution, which can be
used for any downstream inference task (credible intervals, predictions,
derived quantities) while automatically accounting for model
uncertainty.

## Implementation details

### C++ backend

The computationally intensive operations, histogram-based divergence
estimation and power-scaling across the \\\alpha\\ grid, are implemented
in C++ with OpenMP parallelization. The `power_scale_divergence_cpp()`
function parallelizes over \\\alpha\\ values: for each \\\alpha\\, it
computes the raw log-importance-ratios, applies PSIS smoothing via the
shared `psis_smooth()` function from `psis_internal.h`, normalizes the
weights, and computes the JS divergence per parameter against the
uniform-weighted original posterior. The `pairwise_divergence_cpp()`
function parallelizes over the \\\binom{K}{2} \times D\\ parameter-pair
combinations, with each work unit computing the chosen divergence
measure (JS, Wasserstein-2, or overlap) for one parameter and one pair
of models.

### Posterior extraction across fit classes

[`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md)
extracts posterior samples from all nine lucifer fit classes through a
unified internal function. Demonoid objects use `$Posterior2` (the
post-burnin, thinned posterior) when available, falling back to
`$Posterior1`. PMC objects similarly use `$Posterior2`. Variational
Bayes, Laplace approximation, iterative quadrature, and Bayesian
quadrature objects use `$Posterior` (the posterior sample matrix
generated by sampling importance resampling or direct quadrature). SMC
and ABC objects use `$Posterior` directly. SBI objects check for
`$Posterior` first, then fall back to `$mcmc_fit$Posterior2` for
MCMC-refined SBI fits.

### Module selection and graceful degradation

Not all modules can run on every input. Single-model analyses cannot
perform cross-model comparison or aggregation; fits without
`Model`/`Data` cannot run power-scaling, conflict, or influence
analysis; models without `Data$PGF` cannot run conflict diagnostics.
Rather than failing,
[`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md)
skips unavailable modules with informative messages and returns `NULL`
for those components. This design allows the function to extract maximum
information from whatever inputs are provided.

## Worked example 1: power-scaling sensitivity

This example demonstrates the interplay between power-scaling
sensitivity and prior-data conflict detection in a setting where the
prior and data disagree. We specify a normal likelihood with unknown
mean and variance, and an informative prior on the mean centered at zero
with standard deviation 2. The true data-generating mean is 5, so the
prior and data pull in different directions. Crucially, with \\n = 50\\
observations from N(5, 1), the data dominate the posterior so thoroughly
that perturbing the prior (the power-scaling question) barely moves the
marginal posteriors. Power-scaling therefore classifies `mu` as
“likelihood-dominated” or “robust”, which is correct: the posterior is
genuinely insensitive to prior perturbation because the data overwhelm
any prior specification. The prior-data *conflict* itself is detected by
the Evans-Moshonov calibration (Module 3), not by power-scaling. Running
both modules on the same fit illustrates this important distinction.

``` r

library(lucifer)

# Simulate data: true mean = 5, true sigma = 1
n <- 50
set.seed(42)
y <- rnorm(n, mean = 5, sd = 1)
```

``` r

# Model: informative prior on mu centered at 0
Model <- function(parm, Data) {
    mu <- parm[1]
    sigma <- interval(parm[2], 1e-100, Inf)
    parm[2] <- sigma
    ### Log-prior
    # N(0, 2) prior on mu: centered at 0, but data is at 5
    mu.prior <- dnorm(mu, 0, 2, log = TRUE)
    sigma.prior <- dhalfcauchy(sigma, 5, log = TRUE)
    ### Log-likelihood
    LL <- sum(dnorm(Data$y, mu, sigma, log = TRUE))
    ### Log-posterior
    LP <- LL + mu.prior + sigma.prior
    Modelout <- list(LP = LP, Dev = -2 * LL, Monitor = LP,
                     yhat = rnorm(Data$n, mu, sigma), parm = parm)
    return(Modelout)
}

Data <- list(
    y = y, n = n,
    parm.names = c("mu", "sigma"),
    mon.names = "LP",
    PGF = function(Data) c(rnorm(1, 0, 2), abs(rcauchy(1, 0, 5)))
)

# Fit with MCMC
Initial.Values <- c(mean(y), sd(y))
Fit <- lucifer(Model, Data, Initial.Values,
               Iterations = 20000, Status = 5000,
               Thinning = 10, Algorithm = "NUTS")
```

``` r

# Run power-scaling and conflict modules
rb <- RobustBayes(Fit, Model = Model, Data = Data,
                  modules = c("power", "conflict"))

# One-screen summary
print(rb)

# Full diagnostic tables
summary(rb)

# Power-scaling heatmap:
# x-axis is alpha, y-axis is parameter, fill is sqrt(JS)
# Facets separate prior perturbation from likelihood perturbation
plot(rb, type = "power")
plot(rb, type = "power_density")
```

In this example, the `mu` parameter shows *low* prior sensitivity
despite the N(0, 2) prior disagreeing with the data centered at 5. This
apparent paradox illustrates a key conceptual point: power-scaling
sensitivity measures how much the posterior *moves* when the prior is
perturbed, not whether the prior and data *disagree*. With \\n = 50\\
observations providing strong likelihood information, the posterior for
`mu` sits near \\\bar{y} \approx 5\\ regardless of whether we raise the
prior to \\\pi(\theta)^{0.5}\\ or \\\pi(\theta)^{1.5}\\. The data simply
overpower the prior at any reasonable perturbation strength. The `sigma`
parameter is similarly robust because the half-Cauchy(5) prior is
diffuse and compatible with the true standard deviation of 1. If the
conflict module is run alongside (as configured above), the
Evans-Moshonov calibration correctly flags `mu` as having severe
prior-data conflict: the posterior mean lies deep in the tail of the
N(0, 2) prior, confirming that the prior and data disagree even though
the posterior is insensitive to prior perturbation.

The power-scaling heatmap reveals the sensitivity pattern across the
full \\\alpha\\ grid. Parameters with warm-colored (high divergence)
bands at \\\alpha\\ values away from 1 are sensitive; parameters with
uniformly cool colors are robust. The vertical dashed line at \\\alpha =
1\\ marks the original posterior.

The density overlay (`type = "power_density"`) provides a complementary
view inspired by priorsense (Kallioinen et al. 2024). Each curve
represents the posterior density reweighted to a different \\\alpha\\
value, colored on a gradient from blue (\\\alpha \< 1\\, component
weakened) through dark navy (\\\alpha = 1\\, base posterior) to orange
(\\\alpha \> 1\\, component strengthened). When all curves overlap
tightly the posterior is insensitive to that perturbation; when they fan
apart the posterior depends materially on the component’s specification.
Dashed lines (when present) flag \\\alpha\\ values where the PSIS Pareto
\\\hat{k}\\ exceeds 0.7, indicating that the importance sampling
approximation may be unreliable at that perturbation strength.

## Worked example 2: cross-model comparison

When multiple model specifications are scientifically plausible,
[`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md)
quantifies how posteriors diverge across the model family. This is
particularly useful in applied settings where the analyst must choose
between, say, different link functions, error distributions, or
covariate sets.

``` r

# Suppose we have three models for the same data:
# M1: Normal errors with N(0, 10) prior on mu
# M2: Normal errors with N(0, 1) prior on mu (tighter)
# M3: Normal errors with N(3, 2) prior on mu (shifted)

Model1 <- function(parm, Data) {
    mu <- parm[1]; sigma <- interval(parm[2], 1e-100, Inf)
    parm[2] <- sigma
    LL <- sum(dnorm(Data$y, mu, sigma, log = TRUE))
    LP <- LL + dnorm(mu, 0, 10, log = TRUE) + dhalfcauchy(sigma, 5, log = TRUE)
    list(LP = LP, Dev = -2 * LL, Monitor = LP,
         yhat = rnorm(Data$n, mu, sigma), parm = parm)
}

Model2 <- function(parm, Data) {
    mu <- parm[1]; sigma <- interval(parm[2], 1e-100, Inf)
    parm[2] <- sigma
    LL <- sum(dnorm(Data$y, mu, sigma, log = TRUE))
    LP <- LL + dnorm(mu, 0, 1, log = TRUE) + dhalfcauchy(sigma, 5, log = TRUE)
    list(LP = LP, Dev = -2 * LL, Monitor = LP,
         yhat = rnorm(Data$n, mu, sigma), parm = parm)
}

Model3 <- function(parm, Data) {
    mu <- parm[1]; sigma <- interval(parm[2], 1e-100, Inf)
    parm[2] <- sigma
    LL <- sum(dnorm(Data$y, mu, sigma, log = TRUE))
    LP <- LL + dnorm(mu, 3, 2, log = TRUE) + dhalfcauchy(sigma, 5, log = TRUE)
    list(LP = LP, Dev = -2 * LL, Monitor = LP,
         yhat = rnorm(Data$n, mu, sigma), parm = parm)
}

# Fit all three
Fit1 <- lucifer(Model1, Data, Initial.Values,
                Iterations = 20000, Status = 5000,
                Thinning = 10, Algorithm = "NUTS")
Fit2 <- lucifer(Model2, Data, Initial.Values,
                Iterations = 20000, Status = 5000,
                Thinning = 10, Algorithm = "NUTS")
Fit3 <- lucifer(Model3, Data, Initial.Values,
                Iterations = 20000, Status = 5000,
                Thinning = 10, Algorithm = "NUTS")
```

``` r

fits <- list(M1 = Fit1, M2 = Fit2, M3 = Fit3)

# Cross-model comparison using Jensen-Shannon divergence
rb <- RobustBayes(fits, modules = "compare", divergence = "js")
summary(rb)

# Tile heatmap of pairwise energy distances
plot(rb, type = "divergence")
```

The divergence tile plot displays the \\K \times K\\ matrix of pairwise
energy distances, with numeric labels overlaid on each cell. Models that
produce similar posteriors will have small values (cool colors); models
that disagree substantially will have large values (warm colors). The
parameter sensitivity ranking in the summary output identifies which
parameters drive the disagreement, helping the analyst focus on the
aspects of the model that matter most.

With \\n = 50\\ observations from N(5, 1), Models M1 (diffuse prior) and
M3 (prior near the data) should produce similar posteriors because the
data dominate; Model M2 (tight prior at 0) should produce a posterior
for `mu` that is pulled toward zero, generating substantial divergence
from the other two models.

## Worked example 3: prior-data conflict

This example constructs a deliberately misspecified prior to demonstrate
the conflict detection diagnostics.

``` r

# Prior: mu ~ N(0, 1) -- tight and far from data
# Data: y ~ N(5, 1) -- posterior mean is near 5
# The prior concentrates at 0 while the data demand mu near 5
Model.conflict <- function(parm, Data) {
    mu <- parm[1]
    sigma <- interval(parm[2], 1e-100, Inf)
    parm[2] <- sigma
    LL <- sum(dnorm(Data$y, mu, sigma, log = TRUE))
    LP <- LL + dnorm(mu, 0, 1, log = TRUE) + dhalfcauchy(sigma, 5, log = TRUE)
    list(LP = LP, Dev = -2 * LL, Monitor = LP,
         yhat = rnorm(Data$n, mu, sigma), parm = parm)
}

Data.conflict <- Data
Data.conflict$PGF <- function(Data) c(rnorm(1, 0, 1), abs(rcauchy(1, 0, 5)))

Fit.conflict <- lucifer(Model.conflict, Data.conflict, Initial.Values,
                        Iterations = 20000, Status = 5000,
                        Thinning = 10, Algorithm = "NUTS")
```

``` r

rb <- RobustBayes(Fit.conflict, Model = Model.conflict, Data = Data.conflict,
                  modules = "conflict")
summary(rb)
plot(rb, type = "conflict")
```

The conflict forest plot shows, for each parameter, the 95% prior
interval (green) alongside the posterior location. Parameters flagged
with conflict (red X markers) indicate where the posterior has
concentrated in the tails of the prior. In this example, `mu` should be
flagged because the posterior mean (near 5) is far in the tail of the
N(0, 1) prior, while `sigma` should show no conflict because the
half-Cauchy(5) prior easily accommodates \\\sigma \approx 1\\.

The Box prior predictive p-value provides the global assessment. Drawing
1000 parameter vectors from the prior and evaluating the deviance at
each, the observed deviance (at the posterior mean near 5) is compared
against this prior-predictive reference distribution. The extreme
p-value confirms the global prior-data conflict.

## Worked example 4: observation influence

Identifying influential observations is essential for assessing
posterior robustness to data perturbations. This example shows how to
detect and visualize outlier influence.

``` r

# Add an outlier to the data
y.outlier <- c(y, 50)  # one extreme observation
n.out <- length(y.outlier)

Data.outlier <- list(
    y = y.outlier, n = n.out,
    parm.names = c("mu", "sigma"),
    mon.names = "LP",
    PGF = function(Data) c(rnorm(1, 0, 10), abs(rcauchy(1, 0, 5)))
)

Model.outlier <- function(parm, Data) {
    mu <- parm[1]
    sigma <- interval(parm[2], 1e-100, Inf)
    parm[2] <- sigma
    LL <- sum(dnorm(Data$y, mu, sigma, log = TRUE))
    LP <- LL + dnorm(mu, 0, 10, log = TRUE) + dhalfcauchy(sigma, 5, log = TRUE)
    list(LP = LP, Dev = -2 * LL, Monitor = LP,
         yhat = rnorm(Data$n, mu, sigma), parm = parm)
}

Fit.outlier <- lucifer(Model.outlier, Data.outlier,
                       c(mean(y.outlier), sd(y.outlier)),
                       Iterations = 20000, Status = 5000,
                       Thinning = 10, Algorithm = "NUTS")
```

To run the influence module, we need the pointwise log-likelihood
matrix. For models where the user monitors per-observation
log-likelihoods, this matrix is directly available. Alternatively,
[`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md)
can reconstruct it from Model evaluations if `Model` and `Data` are
supplied.

``` r

# Option A: supply log_lik matrix directly if available
# rb <- RobustBayes(Fit.outlier, modules = "influence",
#                   log_lik = log_lik_matrix)

# Option B: let RobustBayes compute it from Model/Data
rb <- RobustBayes(Fit.outlier, Model = Model.outlier, Data = Data.outlier,
                  modules = "influence")

summary(rb)
plot(rb, type = "influence")
```

The influence plot has two panels. The top panel shows
Pareto-\\\hat{k}\\ values for each observation, color-coded by
reliability category (blue for good, yellow for acceptable, orange for
problematic, red for very bad). The horizontal dashed line at \\k =
0.7\\ marks the reliability threshold. The bottom panel shows the
maximum absolute mean shift per observation, with flagged observations
highlighted in red. In this example, observation 51 (the outlier at 50)
should stand out dramatically in both panels: its Pareto-\\\hat{k}\\
will be high (indicating strong influence), and its removal would shift
the posterior mean substantially.

## Worked example 5: model aggregation

When no single model dominates, stacking produces an optimally weighted
predictive mixture that hedges across model uncertainty.

``` r

fits <- list(M1 = Fit1, M2 = Fit2, M3 = Fit3)

# Compute all three weighting strategies
rb <- RobustBayes(fits, modules = "aggregate",
                  aggregate.method = "stacking")

# Compare weights side by side
summary(rb)

# Grouped bar chart: BMA vs stacking vs pseudo-BMA+
plot(rb, type = "weights")
```

The weights bar chart displays the three weighting strategies side by
side for each model. A common pattern is for BMA to concentrate weight
on one model (the one with highest marginal likelihood) while stacking
distributes weight more evenly across models that contribute
complementary predictive information. When stacking and BMA agree, it
suggests that one model genuinely dominates; when they disagree,
stacking’s weights are typically preferable for prediction.

The summary output also reports the pointwise elpd-LOO for each model,
which provides the raw material for understanding why the weights
differ. A model with high average elpd but high variance across
observations may receive less stacking weight than a model with moderate
but consistent elpd.

## Full dashboard

Running all modules on a single analysis produces a comprehensive
diagnostic dashboard. Each module’s results are displayed in a single
panel of a multi-panel grid.

``` r

rb <- RobustBayes(fits, Model = Model1, Data = Data)

# Default dashboard: one panel per available module
plot(rb)

# Individual plot types for deeper investigation
plot(rb, type = "power")
plot(rb, type = "power_density")
plot(rb, type = "divergence")
plot(rb, type = "conflict")
plot(rb, type = "influence")
plot(rb, type = "weights")

# Subset parameters with regex
plot(rb, type = "power", Parms = "mu")
```

## Interpreting results: a practical guide

### When to worry about power-scaling sensitivity

A sensitivity score above the threshold (default 0.05) does not
automatically invalidate the analysis. Prior-dominated parameters are
expected when the sample size is small relative to the parameter
dimensionality; the important question is whether the prior
specification is scientifically defensible. If the informative prior
encodes genuine domain knowledge (e.g., a physical constraint), prior
dominance is appropriate. If the prior is a convenient default, prior
sensitivity is a warning that the conclusions may not survive under
alternative reasonable priors.

Likelihood sensitivity is less concerning in most applications, since it
indicates that the data strongly inform the parameter. However, high
likelihood sensitivity combined with model misspecification can produce
posteriors that are precise but wrong.

### When to worry about cross-model divergence

Small marginal divergences (JS \\\< 0.01\\, overlap \\\> 0.95\\)
indicate that the conclusions are robust to the model specifications
considered. Large divergences warrant investigation: which parameters
diverge, and do the differences affect substantive conclusions? If the
scientific question depends on a parameter whose posteriors diverge
across models, model averaging (via stacking) provides a principled way
to account for this uncertainty.

### When to worry about prior-data conflict

Severe conflict (\\C_d \< 0.01\\) indicates that the prior and data are
fundamentally incompatible for that parameter. The appropriate response
depends on the context. If the prior was intended to be weakly
informative, the conflict suggests a coding error or an unintended
consequence of the prior parameterization. If the prior encodes genuine
prior knowledge, the conflict may indicate that the data are
inconsistent with previous findings, which is itself scientifically
interesting.

### When to worry about influential observations

Observations with Pareto \\\hat{k} \> 0.7\\ deserve individual
investigation. They may be legitimate extreme values, in which case a
more robust error model (e.g., Student-t errors) may be appropriate.
They may be data entry errors, in which case correction or removal is
warranted. Or they may be observations from a different population,
suggesting the need for a mixture model or a more flexible likelihood
specification.

## Computational considerations

The most expensive component of
[`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md)
is the LP/Dev decomposition, which requires \\S\\ Model evaluations (one
per posterior sample). For models with expensive likelihood computation,
this can dominate the total runtime. The `cache` argument to the
internal decomposition function avoids redundant computation when
multiple modules share the same decomposition.

The C++ backend for histogram-based divergence computation is
OpenMP-parallelized and handles the \\A \times K\\ grid of power-scaling
divergences and the \\\binom{K}{2} \times D\\ pairwise divergences
efficiently. For typical analyses with \\S \sim 1000\\, \\K \sim 3\\
models, \\D \sim 10\\ parameters, and \\A \sim 20\\ alpha values, the
C++ computation adds negligible time relative to the Model evaluations.

Stacking optimization via
[`constrOptim()`](https://rdrr.io/r/stats/constrOptim.html) is fast for
small numbers of models (\\K \< 20\\) but may converge slowly for larger
model families. The pseudo-BMA+ bootstrap (\\B = 1000\\ replicates) is
embarrassingly parallel in principle, though the current implementation
uses a sequential loop; for most applications with \\K \< 10\\ models
and \\N \< 10000\\ observations, this completes in under a second.

## References

**\[1\]** Kallioinen, N., Paananen, T., Buerkner, P.-C., and Vehtari, A.
(2024). Detecting and diagnosing prior and likelihood sensitivity with
power-scaling. *Statistics and Computing*, 34, 57. [DOI:
10.1007/s11222-023-10366-5](https://doi.org/10.1007/s11222-023-10366-5)

**\[2\]** Box, G.E.P. (1980). Sampling and Bayes’ inference in
scientific modelling and robustness. *JRSS-A*, 143(4), 383-430.

**\[3\]** Evans, M. and Moshonov, H. (2006). Checking for prior-data
conflict. *Bayesian Analysis*, 1(4), 893-914. [DOI:
10.1214/06-BA129](https://doi.org/10.1214/06-BA129)

**\[4\]** Vehtari, A., Simpson, D., Gelman, A., Yao, Y., and Gabry, J.
(2024). Pareto smoothed importance sampling. *JMLR*, 25(72), 1-58.

**\[5\]** Yao, Y., Vehtari, A., Simpson, D., and Gelman, A. (2018).
Using stacking to average Bayesian predictive distributions. *Bayesian
Analysis*, 13(3), 917-1007. [DOI:
10.1214/17-BA1091](https://doi.org/10.1214/17-BA1091)

**\[6\]** Berger, J.O. (1994). An overview of robust Bayesian analysis.
*Test*, 3, 5-124.

**\[7\]** Zhang, J. and Stephens, M.A. (2009). A new and efficient
estimation method for the generalized Pareto distribution.
*Technometrics*, 51(3), 316-325.

**\[8\]** Bernardo, J.M. and Smith, A.F.M. (1994). *Bayesian Theory*.
John Wiley & Sons. ISBN: 978-0-471-49464-5.

**\[9\]** Gelman, A., Hwang, J., and Vehtari, A. (2014). Understanding
predictive information criteria for Bayesian models. *Statistics and
Computing*, 24(6), 997-1016. [DOI:
10.1007/s11222-013-9416-2](https://doi.org/10.1007/s11222-013-9416-2)

**\[10\]** Sivula, T., Magnusson, M., Matamoros, A.A., and Vehtari, A.
(2022). Uncertainty in Bayesian leave-one-out cross-validation based
model comparison. *arXiv:2008.10296v3*.
