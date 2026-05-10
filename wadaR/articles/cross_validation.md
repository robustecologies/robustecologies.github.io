# Cross-validation for Bayesian models

## Introduction

Model comparison and predictive assessment are central tasks in Bayesian
statistics. Given a set of candidate models \\\\M_1, \ldots, M_K\\\\
fitted to observed data \\y = (y_1, \ldots, y_N)\\, the fundamental
question is: which model will predict future observations best? The
expected log pointwise predictive density (elpd) provides a principled
answer, measuring the out-of-sample predictive accuracy in the natural
logarithmic score. Cross-validation estimates this quantity by
systematically evaluating predictions on held-out data.

The `lucifer` package implements a unified cross-validation framework
that covers standard leave-one-out (LOO) approximation via Pareto
smoothed importance sampling (PSIS), the widely applicable information
criterion (WAIC), exact K-fold cross-validation, leave-future-out
cross-validation (LFO) for time series, and Bayesian model stacking for
combining models. All methods return objects with compatible `pointwise`
data frames, enabling seamless comparison through
[`loo_compare()`](https://robustecologies.github.io/lucifer/reference/loo_compare.md)
and combination through
[`stacking_weights()`](https://robustecologies.github.io/lucifer/reference/stacking_weights.md).

## Theoretical foundations

### The expected log pointwise predictive density

The gold standard for evaluating predictive accuracy is the expected log
pointwise predictive density [\[1\]](#ref1), defined as

\\\text{elpd} = \sum\_{i=1}^{N} \int \log p(\tilde{y}\_i \mid y) \\
f(\tilde{y}\_i) \\ d\tilde{y}\_i\\

where \\f(\tilde{y}\_i)\\ is the true data-generating distribution for a
new observation and

\\p(\tilde{y}\_i \mid y) = \int p(\tilde{y}\_i \mid \theta) \\ p(\theta
\mid y) \\ d\theta\\

is the posterior predictive density. Since \\f\\ is unknown, we must
estimate the elpd from the observed data. Cross-validation provides a
nearly unbiased estimate by using portions of the data as surrogate
“future” observations.

### Leave-one-out cross-validation

The LOO-CV estimate of elpd is

\\\widehat{\text{elpd}}\_{\text{loo}} = \sum\_{i=1}^{N} \log p(y_i \mid
y\_{-i})\\

where \\p(y_i \mid y\_{-i}) = \int p(y_i \mid \theta) \\ p(\theta \mid
y\_{-i}) \\ d\theta\\ is the leave-one-out predictive density obtained
by fitting the model to all data except observation \\i\\. Direct
computation requires \\N\\ model refits, but PSIS provides an efficient
single-fit approximation [\[2\]](#ref2) [\[3\]](#ref3).

The key insight is that the LOO predictive density can be rewritten as
an importance sampling estimate:

\\p(y_i \mid y\_{-i}) = \frac{\int p(y_i \mid \theta) \\ p(\theta \mid
y\_{-i}) \\ d\theta}{1} = E\_{p(\theta \mid y)} \left\[ \frac{1}{p(y_i
\mid \theta)} \right\]^{-1}\\

The importance ratios \\r_i(\theta) = 1 / p(y_i \mid \theta)\\ have
potentially infinite variance, so PSIS stabilizes the upper tail by
fitting a generalized Pareto distribution (GPD) to the largest ratios,
following Zhang and Stephens [\[4\]](#ref4), and replacing extreme
values with smoothed quantiles. The shape parameter \\\hat{k}\\ of the
fitted GPD diagnoses reliability:

- \\\hat{k} \< 0.5\\: the estimate is reliable; the central limit
  theorem for importance ratios holds and the convergence rate is close
  to \\1/\sqrt{S}\\.
- \\0.5 \leq \hat{k} \< 0.7\\: the estimate is acceptable but Monte
  Carlo error may be non-negligible; the convergence rate is slower than
  the parametric rate.
- \\0.7 \leq \hat{k} \< 1.0\\: the estimate is problematic; the
  importance weights are highly variable and the finite-sample bias may
  be substantial.
- \\\hat{k} \geq 1.0\\: the estimate is unreliable; the moment of the
  importance ratios required for the estimate to be consistent does not
  exist, meaning the importance sampling approximation does not converge
  at any rate.

### WAIC

The Widely Applicable Information Criterion (WAIC) of Watanabe
[\[5\]](#ref5) provides an alternative estimate of the elpd that is
asymptotically equivalent to LOO-CV. It is computed as

\\\widehat{\text{elpd}}\_{\text{waic}} = \text{lppd} -
p\_{\text{waic}}\\

where the log pointwise predictive density is \\\text{lppd} =
\sum\_{i=1}^{N} \log \left( \frac{1}{S} \sum\_{s=1}^{S} p(y_i \mid
\theta^s) \right)\\ and the effective number of parameters is
\\p\_{\text{waic}} = \sum\_{i=1}^{N} \text{Var}\_s \left\[ \log p(y_i
\mid \theta^s) \right\]\\, where the variance is computed across the
\\S\\ posterior samples. The latter formula (\\p\_{\text{waic2}}\\) is
preferred by Gelman et al. [\[1\]](#ref1) over the alternative
\\p\_{\text{waic1}} = 2 \sum_i \left\[ \log \left( \frac{1}{S} \sum_s
p(y_i \mid \theta^s) \right) - \frac{1}{S} \sum_s \log p(y_i \mid
\theta^s) \right\]\\ because it more closely matches LOO-CV in practice.

Unlike LOO-PSIS, WAIC provides no per-observation diagnostic. When the
posterior is dominated by a few observations (high leverage), WAIC can
silently produce unreliable estimates. For this reason, Vehtari et
al. [\[2\]](#ref2) recommend LOO-PSIS with Pareto \\\hat{k}\\
diagnostics as the default method, resorting to WAIC only when
computational constraints prevent PSIS computation.

## LOO and WAIC in lucifer

Both
[`LOO()`](https://robustecologies.github.io/lucifer/reference/LOO.md)
and
[`WAIC()`](https://robustecologies.github.io/lucifer/reference/WAIC.md)
accept an \\N \times S\\ matrix of pointwise log-likelihoods, where each
row corresponds to an observation and each column to a posterior sample.
The
[`LOO()`](https://robustecologies.github.io/lucifer/reference/LOO.md)
function calls the C++ PSIS backend with OpenMP parallelization over
observations, while
[`WAIC()`](https://robustecologies.github.io/lucifer/reference/WAIC.md)
is pure R. Both return objects with a `pointwise` data frame suitable
for
[`loo_compare()`](https://robustecologies.github.io/lucifer/reference/loo_compare.md).

``` r

set.seed(42)
N <- 50; S <- 1000
log_lik <- matrix(rnorm(N * S, mean = -2, sd = 0.5), nrow = N, ncol = S)

# LOO-PSIS
loo_result <- LOO(log_lik)
print(loo_result)
summary(loo_result)
plot(loo_result)

# WAIC
waic_result <- WAIC(log_lik)
print(waic_result)
plot(waic_result)
```

The `summary.loo()` method flags observations with Pareto \\\hat{k} \geq
0.7\\, reports the proportion of reliable estimates per diagnostic
category, and provides actionable guidance on whether to trust the
approximation or switch to an alternative method. The `plot.loo()`
method visualizes the Pareto \\\hat{k}\\ values per observation, colored
by diagnostic category, with threshold lines at 0.5, 0.7, and 1.0.

## When LOO fails: the case for structured cross-validation

Standard LOO-CV assumes that observations are conditionally independent
given the parameters, so that removing a single observation produces a
valid posterior that does not leak information about the held-out value.
This assumption breaks down in two important classes of models.

### Time series and temporal dependence

For time series models, the likelihood factors as \\p(y \mid \theta) =
p(y_1 \mid \theta) \prod\_{t=2}^{N} p(y_t \mid y\_{1:t-1}, \theta)\\,
reflecting the sequential dependence. Leaving out observation \\y_t\\
from the likelihood while retaining \\y\_{t+1}, y\_{t+2}, \ldots\\
creates a paradox: the future observations carry information about
\\y_t\\ through the temporal structure, so the “leave-one-out” posterior
\\p(\theta \mid y\_{-t})\\ is not truly ignorant of \\y_t\\. The
resulting elpd estimate is optimistically biased because the held-out
prediction benefits from information leakage through future data
[\[6\]](#ref6).

Burkner, Gabry and Vehtari [\[6\]](#ref6) demonstrated this bias on
autoregressive models and proposed leave-future-out cross-validation as
the principled alternative. Their key result is that for a correctly
specified AR(1) model, the bias of LOO-CV relative to exact
one-step-ahead prediction increases with the autocorrelation coefficient
and can reach several elpd units for strongly autocorrelated series.

### Hierarchical and spatial models

In hierarchical models with shared group-level parameters, leaving out a
single observation within a group still allows the remaining
within-group observations to inform the group effect, which in turn
informs the prediction for the held-out observation. For spatial models,
nearby observations transmit information through the spatial correlation
structure. In both cases, LOO-CV overestimates predictive performance,
and the magnitude of the bias depends on the strength of the
within-group or spatial correlation [\[7\]](#ref7).

K-fold cross-validation with appropriately structured folds (e.g.,
leaving out entire groups, spatial blocks, or temporal blocks) provides
a principled solution by ensuring that the training and test sets are
independent given the model structure.

## Leave-future-out cross-validation for time series

### Theory

Leave-future-out cross-validation (LFO-CV) evaluates predictive
performance by expanding a training window forward through time, making
predictions only about future observations that were not used for
inference [\[6\]](#ref6). The LFO-CV estimate of elpd for
\\M\\-step-ahead prediction is

\\\widehat{\text{elpd}}\_{\text{lfo}} = \sum\_{i=L}^{N-M} \log
p(y\_{i+1:i+M} \mid y\_{1:i})\\

where \\L\\ is the minimum training window and \\p(y\_{i+1:i+M} \mid
y\_{1:i}) = \int p(y\_{i+1:i+M} \mid y\_{1:i}, \theta) \\ p(\theta \mid
y\_{1:i}) \\ d\theta\\ is the posterior predictive density for the
future block given the data up to time \\i\\.

Exact LFO-CV requires refitting the model at each time step, which is
prohibitively expensive. The approximate algorithm of Burkner et
al. [\[6\]](#ref6) reduces the number of refits by using PSIS to
approximate the successive posteriors from a single fit, refitting only
when the importance sampling approximation becomes unreliable.

### Algorithm

The algorithm proceeds as follows. Starting from a model fitted on
\\y\_{1:L}\\, we accumulate importance weights as the training window
expands:

1.  Fit the model on \\y\_{1:L}\\ to obtain posterior samples
    \\\theta^{(1)}, \ldots, \theta^{(S)} \sim p(\theta \mid y\_{1:L})\\.
2.  For \\i = L, L+1, \ldots, N-M\\:
    1.  Compute cumulative log-importance-ratios since the last refit at
        \\i^\*\\: \\\log r_i(\theta^{(s)}) = \sum\_{j=i^\*+1}^{i} \log
        p(y_j \mid y\_{1:j-1}, \theta^{(s)})\\
    2.  Apply PSIS to the importance ratios to obtain normalized weights
        \\w_i^{(s)}\\ and the Pareto \\\hat{k}\\ diagnostic.
    3.  If \\\hat{k} \> k\_{\text{threshold}}\\ (default 0.7): refit the
        model on \\y\_{1:i}\\, reset the importance ratios, and set
        \\i^\* = i\\.
    4.  Compute the weighted predictive log-density: \\\log
        p(y\_{i+1:i+M} \mid y\_{1:i}) \approx \log \sum\_{s=1}^{S}
        w_i^{(s)} \\ p(y\_{i+1:i+M} \mid y\_{1:i}, \theta^{(s)})\\

The number of refits depends on how quickly the posterior evolves as new
data arrives. For well-specified models with stable parameters, the
posterior changes slowly and few refits are needed (typically \< 5% of
the total time steps). For misspecified models or models with
time-varying parameters, more frequent refits indicate that the
importance sampling approximation degrades quickly, which is itself
informative about model adequacy.

### Predictive log-likelihoods via deviance differencing

Computing \\\log p(y\_{i+1} \mid y\_{1:i}, \theta)\\ requires access to
the conditional predictive density, which many model specifications do
not expose directly. The
[`LFO()`](https://robustecologies.github.io/lucifer/reference/LFO.md)
function uses deviance differencing as the default approach: for each
posterior sample \\\theta^{(s)}\\, the model is evaluated on the data
\\y\_{1:i+1}\\ and \\y\_{1:i}\\, and the pointwise log-likelihood is
obtained as

\\\log p(y\_{i+1} \mid y\_{1:i}, \theta^{(s)}) =
-\frac{\text{Dev}(y\_{1:i+1}, \theta^{(s)})}{2} +
\frac{\text{Dev}(y\_{1:i}, \theta^{(s)})}{2}\\

This is exact for models with additive log-likelihoods, which covers the
vast majority of practical cases. For models with non-decomposable
likelihoods (joint distributions, copulas), users can supply a custom
pointwise log-likelihood function via `Data$ll_fun`.

### Usage

The following example fits an AR(1) model and evaluates one-step-ahead
predictive performance using LFO-CV with an initial training window of
10 observations.

``` r

# AR(1) model specification
Model <- function(parm, Data) {
    phi <- parm[1]
    sigma <- exp(parm[2])
    mu <- phi * Data$y[-Data$n]
    LL <- sum(dnorm(Data$y[-1], mu, sigma, log = TRUE))
    LP <- LL + dnorm(phi, 0, 1, log = TRUE) +
        dgamma(sigma, 1, 1, log = TRUE)
    yhat <- c(Data$y[1], mu)
    list(LP = LP, Dev = -2 * LL, Monitor = LP, yhat = yhat, parm = parm)
}

set.seed(123)
y <- arima.sim(list(ar = 0.7), n = 30)
Data <- list(y = as.numeric(y), n = 30,
             parm.names = c("phi", "log.sigma"),
             mon.names = "LP")

# Fit the model
fit <- lucifer(Model, Data, c(0, 0), Iterations = 2000,
               Status = .Machine$integer.max,
               Algorithm = "HARM", Specs = NULL, Chains = 1)

# One-step-ahead LFO-CV with initial window of 10 observations
lfo_result <- LFO(fit, Model, Data, L = 10, M = 1, verbose = FALSE)
```

``` r

print(lfo_result)
summary(lfo_result)
plot(lfo_result)
```

The two-panel plot shows cumulative elpd over time (top panel) with
refit points marked as vertical dashed lines, and Pareto \\\hat{k}\\
diagnostics (bottom panel) colored by reliability category. Refit
locations cluster where the data is most informative about parameter
changes, providing a diagnostic view of model stability over time.

The prediction horizon `M` controls multi-step-ahead evaluation. Setting
`M = 1` gives one-step-ahead predictions (the default and most common
choice); `M = 4` evaluates four-step-ahead predictive performance, which
is more stringent and directly relevant for forecasting applications.
Note that with \\M \> 1\\, the predictive block \\y\_{i+1:i+M}\\ is
evaluated jointly, so that the total number of evaluated predictions is
\\N - L - M + 1\\.

The `L` parameter sets the minimum training window. It defaults to
\\\max(10, 2d)\\ where \\d\\ is the number of parameters, ensuring that
the initial posterior has seen enough data relative to the model
complexity. Too small a value of \\L\\ yields unreliable initial fits;
too large a value wastes data that could be used for evaluation. A
reasonable heuristic is to set \\L\\ to the smallest sample size at
which the model can produce stable posterior estimates.

### PSIS-only mode

When a log-likelihood matrix is available but refitting is impractical,
[`LFO()`](https://robustecologies.github.io/lucifer/reference/LFO.md)
can operate in PSIS-only mode by accepting a matrix directly. In this
mode, importance weights are accumulated without the possibility of
refitting, so large Pareto \\\hat{k}\\ values cannot be corrected. This
mode is useful for rapid exploratory comparisons but should be
complemented by full LFO-CV with refitting for final model selection.

``` r

# PSIS-only mode using the log-likelihood matrix from the LOO example
lfo_approx <- LFO(log_lik, L = 20, M = 1, verbose = FALSE)
print(lfo_approx)
```

## K-fold cross-validation

### Theory

K-fold cross-validation partitions the \\N\\ observations into \\K\\
non-overlapping folds \\F_1, \ldots, F_K\\, refits the model \\K\\ times
(each time leaving out one fold as the test set), and evaluates the
predictive log density for each held-out observation. The K-fold
estimate of elpd is

\\\widehat{\text{elpd}}\_{\text{kfold}} = \sum\_{k=1}^{K} \sum\_{i \in
F_k} \log p(y_i \mid y\_{-F_k})\\

where \\p(y_i \mid y\_{-F_k}) = \int p(y_i \mid \theta) \\ p(\theta \mid
y\_{-F_k}) \\ d\theta\\ is the predictive density for observation \\i\\
under the posterior obtained from all data except fold \\k\\.

Unlike LOO-PSIS, K-fold CV is exact in the sense that it does not rely
on importance sampling approximations. It works even when Pareto
\\\hat{k}\\ diagnostics indicate that PSIS-LOO is unreliable. The
tradeoff is computational cost: K model refits are required, making
K-fold CV \\K\\ times more expensive than a single fit (but \\N/K\\
times cheaper than exact LOO via refitting).

The choice of \\K\\ balances bias and variance. Larger \\K\\ reduces
bias (each training set is closer to the full data) but increases
variance (the \\K\\ training sets overlap more, making the fold-level
estimates correlated). \\K = 10\\ is a common default that provides a
reasonable tradeoff; \\K = N\\ gives exact LOO via refitting; \\K = 2\\
gives split-half validation with maximum diversity between folds
[\[8\]](#ref8).

### Fold assignment strategies

The default random fold assignment works well for models with
independent observations. For structured data, the fold assignment
should respect the dependence structure:

- **Blocked folds for time series.** Assigning contiguous time blocks to
  each fold prevents temporal leakage between training and test sets.
  For instance, with \\N = 100\\ and \\K = 5\\, observations 1–20 form
  fold 1, 21–40 form fold 2, and so on. This is less powerful than
  LFO-CV (which uses an expanding window) but simpler to implement for
  non-sequential models.

- **Group-level folds for hierarchical models.** Assigning all
  observations within the same group to the same fold prevents
  within-group information leakage. This yields leave-one-group-out CV
  when \\K\\ equals the number of groups.

- **Spatial blocking.** For spatial models, observations within the same
  spatial block should be assigned to the same fold, where blocks are
  defined to minimize spatial correlation between training and test
  sets.

### Usage

``` r

# Normal model
Model <- function(parm, Data) {
    mu <- parm[1]
    sigma <- exp(parm[2])
    LL <- sum(dnorm(Data$y, mu, sigma, log = TRUE))
    LP <- LL + dnorm(mu, 0, 100, log = TRUE) +
        dgamma(sigma, 1, 1, log = TRUE)
    yhat <- rep(mu, Data$n)
    list(LP = LP, Dev = -2 * LL, Monitor = LP, yhat = yhat, parm = parm)
}

set.seed(456)
Data <- list(y = rnorm(30, 5, 2), n = 30,
             parm.names = c("mu", "log.sigma"),
             mon.names = "LP")

fit <- lucifer(Model, Data, c(0, 0), Iterations = 2000,
               Status = .Machine$integer.max,
               Algorithm = "HARM", Specs = NULL, Chains = 1)

# 5-fold CV
kfold_result <- Kfold(Model, Data, fit = fit, K = 5, verbose = FALSE)
```

``` r

print(kfold_result)
summary(kfold_result)
plot(kfold_result)
plot(kfold_result, type = "folds")
```

The [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method with
`type = "pointwise"` (default) shows pointwise elpd values colored by
fold assignment. The `type = "folds"` variant shows the total elpd
contribution from each fold as a bar chart, which is useful for
identifying folds that are substantially harder to predict than others,
potentially indicating heterogeneity or outlier clusters in the data.

Custom fold assignments can be provided to implement blocked or
group-level designs:

``` r

# Blocked folds for time series (contiguous blocks)
N <- 100; K <- 5
folds <- rep(seq_len(K), each = N / K)

# Group-level folds (leave-one-group-out)
group_membership <- Data$group  # integer vector of group IDs
folds <- as.integer(as.factor(group_membership))
```

### Data subsetting contract

The
[`Kfold()`](https://robustecologies.github.io/lucifer/reference/Kfold.md)
and
[`LFO()`](https://robustecologies.github.io/lucifer/reference/LFO.md)
functions rely on an internal data subsetting engine
(`.cv_subset_data()`) that creates training-set data lists by indexing
into the original Data. The engine follows a simple heuristic: numeric
vectors and matrix rows of length/nrow equal to \\N\\ are subsetted;
scalars, parameter names, and metadata fields are preserved. For models
with complex data structures that this heuristic cannot handle, users
can supply a custom subsetting function via `Data$subset_fn`:

``` r

Data$subset_fn <- function(Data, indices) {
    d <- Data
    d$y <- Data$y[indices]
    d$X <- Data$X[indices, , drop = FALSE]
    d$weights <- Data$weights[indices]
    d$n <- length(indices)
    d
}
```

## Model comparison with loo_compare

### Theory

Comparing models requires not only ranking them by their elpd estimates
but also quantifying the uncertainty in that ranking. Given two models
with pointwise elpd estimates \\\widehat{\text{elpd}}\_{1,i}\\ and
\\\widehat{\text{elpd}}\_{2,i}\\ (\\i = 1, \ldots, N\\), the difference
in total elpd is \\\Delta = \sum_i (\widehat{\text{elpd}}\_{1,i} -
\widehat{\text{elpd}}\_{2,i})\\.

The standard error of this difference must account for the positive
correlation between the two sets of pointwise estimates, which are
evaluated on the same observations. The naive standard error
\\\sqrt{\text{se}\_1^2 + \text{se}\_2^2}\\ treats the estimates as
independent and is always an overestimate. The correct formula is

\\\text{se}\_{\text{diff}} = \sqrt{N \cdot
\text{Var}\left(\widehat{\text{elpd}}\_{1,i} -
\widehat{\text{elpd}}\_{2,i}\right)}\\

which is smaller than the naive estimate because \\\text{Var}(A - B) =
\text{Var}(A) + \text{Var}(B) - 2\text{Cov}(A, B)\\ and the covariance
is typically positive [\[2\]](#ref2). This correlated standard error
makes the comparison more powerful: differences that would be
insignificant under naive standard errors may be clearly significant
under the correct analysis.

### Usage

The
[`loo_compare()`](https://robustecologies.github.io/lucifer/reference/loo_compare.md)
function accepts any combination of `loo`, `waic`, `kfold`, and `lfo`
objects. Models are ranked by elpd (best first), and differences are
computed relative to the best model.

``` r

set.seed(789)
N <- 50; S <- 1000

# Three models with complementary predictive strengths.
# The shared `obs` term creates observations where model 1 excels and
# others where model 2 is better, producing heterogeneous pointwise
# elpd patterns and realistic correlated standard errors.
obs <- rnorm(N, 0, 0.3)
mu1 <- -1.5 + obs + rnorm(N, 0, 0.3)
mu2 <- -1.7 - obs + rnorm(N, 0, 0.3)
mu3 <- -2.8 + rnorm(N, 0, 0.5)

log_lik1 <- matrix(rnorm(N * S, mu1, 0.5), nrow = N)
log_lik2 <- matrix(rnorm(N * S, mu2, 0.5), nrow = N)
log_lik3 <- matrix(rnorm(N * S, mu3, 0.8), nrow = N)

loo1 <- LOO(log_lik1)
loo2 <- LOO(log_lik2)
loo3 <- LOO(log_lik3)

comp <- loo_compare(best = loo1, mid = loo2, worst = loo3)
print(comp)
plot(comp)
```

The forest plot shows elpd differences relative to the best model with
\\\pm 2\\ standard error intervals. A model whose interval includes zero
cannot be reliably distinguished from the best model at the
approximately 95% level. When the difference is large relative to its
standard error (say, \\\|\Delta\| / \text{se}\_{\text{diff}} \> 2\\),
there is strong evidence that the models differ in predictive
performance.

Comparison across different criteria types is supported: one can compare
a LOO estimate for one model against a K-fold estimate for another,
provided both were evaluated on the same observations. This flexibility
is useful when some models have LOO-PSIS problems requiring K-fold,
while others have reliable LOO estimates.

``` r

# Mixed comparison: LOO for model 1, K-fold for model 2
comp <- loo_compare(m1 = loo_result, m2 = kfold_result)
```

## Model averaging with Bayesian stacking

### Theory

Model selection picks a single model and discards the rest, wasting
information about the predictive value of non-selected models. Model
averaging addresses this by combining predictions from multiple models
with data-determined weights. The classical Bayesian model averaging
(BMA) approach weights models by their marginal likelihoods, but this is
problematic when all candidate models are misspecified (as they always
are in practice), because marginal likelihoods favor models that are
closest to the data-generating process in a particular KL sense that may
not align with predictive performance [\[9\]](#ref9).

Bayesian stacking [\[9\]](#ref9) directly optimizes predictive
performance by finding weights \\w = (w_1, \ldots, w_K)\\ that maximize
the leave-one-out predictive log density of the combined model:

\\\hat{w} = \arg\max\_{w} \sum\_{i=1}^{N} \log \left( \sum\_{k=1}^{K}
w_k \exp\left(\widehat{\text{elpd}}\_{k,i}\right) \right) \quad
\text{subject to} \quad w_k \geq 0, \quad \sum_k w_k = 1\\

This is a convex optimization problem on the simplex. The
[`stacking_weights()`](https://robustecologies.github.io/lucifer/reference/stacking_weights.md)
function solves it via unconstrained optimization on the softmax
parameterization \\w_k = \exp(\alpha_k) / \sum_j \exp(\alpha_j)\\, which
automatically satisfies the simplex constraints.

Stacking accounts for model complementarity: a model with lower overall
elpd may still receive positive weight if it captures aspects of the
data that the other models miss. This makes stacking fundamentally
different from model selection, which ignores the possibility that
different models may be best at predicting different subsets of the
data.

### Pseudo-BMA and pseudo-BMA+

Two simpler alternatives are available for rapid computation:

**Pseudo-BMA** assigns weights proportional to
\\\exp(\widehat{\text{elpd}}\_k)\\ where \\\widehat{\text{elpd}}\_k =
\sum_i \widehat{\text{elpd}}\_{k,i}\\ is the total elpd of model \\k\\.
This is equivalent to a Bayesian model averaging approximation using the
LOO log-likelihood in place of the marginal log-likelihood, and tends to
assign weight 1 to the model with highest elpd and 0 to all others when
the elpd differences are large.

**Pseudo-BMA+** adds Bayesian bootstrap uncertainty to the pseudo-BMA
weights [\[9\]](#ref9). For each bootstrap replicate, the pointwise elpd
values are reweighted by Dirichlet(1,…,1) random weights, and the best
model under the reweighted criterion is recorded. The final weights are
the proportion of bootstrap replicates in which each model was best.
This produces smoother weights that are less sensitive to sampling
variability in the pointwise elpd estimates.

### Usage

``` r

# Stacking weights (recommended)
w_stack <- stacking_weights(m1 = loo1, m2 = loo2, m3 = loo3,
                            method = "stacking")
print(w_stack)

# Pseudo-BMA (faster, less optimal)
w_pbma <- stacking_weights(m1 = loo1, m2 = loo2, m3 = loo3,
                           method = "pseudobma")
print(w_pbma)

# Pseudo-BMA+ with Bayesian bootstrap
set.seed(999)
w_pbma_plus <- stacking_weights(m1 = loo1, m2 = loo2, m3 = loo3,
                                method = "pseudobma+")
print(w_pbma_plus)
```

The stacking weights can be used to form a combined predictive
distribution \\p\_{\text{stack}}(\tilde{y} \mid y) = \sum_k \hat{w}\_k
\\ p_k(\tilde{y} \mid y)\\, which generally outperforms any individual
model in terms of log predictive density. In practice, use the weights
to average posterior predictive samples: draw from each model’s
posterior predictive proportionally to its stacking weight.

## Practical workflow

### Step 1: start with LOO-PSIS

For any model comparison task, begin by computing LOO-PSIS for each
candidate model. This is computationally free (it uses the existing
posterior samples) and provides diagnostic information about the
reliability of the approximation.

``` r

loo_m1 <- LOO(log_lik_model1)
loo_m2 <- LOO(log_lik_model2)
summary(loo_m1)
summary(loo_m2)
```

### Step 2: check diagnostics

If all Pareto \\\hat{k} \< 0.7\\, the LOO estimates are reliable and you
can proceed to model comparison. If some observations have problematic
diagnostics, the appropriate next step depends on the data structure.

### Step 3: choose the right CV method

For independent observations with problematic LOO diagnostics, use
K-fold CV:

``` r

kfold_m1 <- Kfold(Model1, Data, fit = fit_m1, K = 10)
```

For time series data, use LFO-CV regardless of LOO diagnostics:

``` r

lfo_m1 <- LFO(fit_m1, Model1, Data, L = 20, M = 1)
```

### Step 4: compare and combine

``` r

# Compare models
comp <- loo_compare(m1 = loo_m1, m2 = loo_m2)
print(comp)

# If no single model dominates, compute stacking weights
w <- stacking_weights(m1 = loo_m1, m2 = loo_m2)
print(w)
```

## Decision table

The following table summarizes which cross-validation method to use in
different scenarios.

| Scenario | Method | Function |
|:---|:---|---:|
| Independent obs, all k \< 0.7 | LOO-PSIS | LOO() |
| Independent obs, some k \> 0.7 | LOO-PSIS (interpret cautiously) | LOO() |
| Independent obs, some k \> 1.0 | K-fold CV (K = 10) | Kfold() |
| Time series (any k) | LFO-CV | LFO() |
| Hierarchical / grouped data | K-fold CV (group-level folds) | Kfold(folds = group_folds) |
| Spatial data | K-fold CV (spatial block folds) | Kfold(folds = spatial_folds) |
| Quick comparison, any data | WAIC | WAIC() |

Cross-validation method selection guide {.table}

## Computational considerations

The computational cost of each method differs substantially. LOO-PSIS
and WAIC are essentially free, requiring only a pass through the
existing posterior samples with no model refitting. K-fold CV requires
\\K\\ model refits, each on approximately \\(1 - 1/K) \cdot N\\
observations. LFO-CV requires a variable number of refits, typically
much fewer than \\N\\, depending on how quickly the posterior evolves.

For models that are expensive to fit (e.g., Bayesian neural networks,
large spatial models, models with C++ particle filter likelihoods), the
refitting cost can dominate. Several strategies can mitigate this:

- **Warm-starting.** Both
  [`Kfold()`](https://robustecologies.github.io/lucifer/reference/Kfold.md)
  and
  [`LFO()`](https://robustecologies.github.io/lucifer/reference/LFO.md)
  use the most recent posterior as the starting point for each refit,
  reducing the number of iterations needed for convergence. Pass reduced
  `Iterations` via `fit_args` or `refit_args`.

- **Parallel folds.**
  [`Kfold()`](https://robustecologies.github.io/lucifer/reference/Kfold.md)
  supports parallel execution of folds via `CPUs > 1` on Unix systems.

- **LFO threshold tuning.** Increasing `k_threshold` from 0.7 to, say,
  0.9 reduces the number of LFO refits at the cost of potentially less
  accurate importance weights. Monitor the resulting Pareto \\\hat{k}\\
  values to assess whether the tradeoff is acceptable.

- **PSIS-only LFO.** When refitting is impractical,
  [`LFO()`](https://robustecologies.github.io/lucifer/reference/LFO.md)
  can operate in PSIS-only mode by accepting a log-likelihood matrix
  directly. The Pareto \\\hat{k}\\ diagnostics indicate where the
  approximation is unreliable, even if refitting cannot be done.

## References

**\[1\]** Gelman, A., Hwang, J., and Vehtari, A. (2014). “Understanding
predictive information criteria for Bayesian models.” *Statistics and
Computing*, 24(6), 997–1016.
[doi:10.1007/s11222-013-9416-2](https://doi.org/10.1007/s11222-013-9416-2)

**\[2\]** Vehtari, A., Gelman, A., and Gabry, J. (2017). “Practical
Bayesian model evaluation using leave-one-out cross-validation and
WAIC.” *Statistics and Computing*, 27(5), 1413–1432.
[doi:10.1007/s11222-016-9696-4](https://doi.org/10.1007/s11222-016-9696-4)

**\[3\]** Vehtari, A., Simpson, D., Gelman, A., Yao, Y., and Gabry, J.
(2024). “Pareto smoothed importance sampling.” *Journal of Machine
Learning Research*, 25(72), 1–58.

**\[4\]** Zhang, J. and Stephens, M. A. (2009). “A new and efficient
estimation method for the generalized Pareto distribution.”
*Technometrics*, 51(3), 316–325.
[doi:10.1198/tech.2009.08017](https://doi.org/10.1198/tech.2009.08017)

**\[5\]** Watanabe, S. (2010). “Asymptotic equivalence of Bayes cross
validation and widely applicable information criterion in singular
learning theory.” *Journal of Machine Learning Research*, 11, 3571–3594.

**\[6\]** Burkner, P.-C., Gabry, J., and Vehtari, A. (2020).
“Approximate leave-future-out cross-validation for Bayesian time series
models.” *Journal of Statistical Computation and Simulation*, 90(14),
2499–2523.
[doi:10.1080/00949655.2020.1783262](https://doi.org/10.1080/00949655.2020.1783262)

**\[7\]** Roberts, D. R., Bahn, V., Ciuti, S., Boyce, M. S., Elith, J.,
Guillera-Arroita, G., et al. (2017). “Cross-validation strategies for
data with temporal, spatial, hierarchical, or phylogenetic structure.”
*Ecography*, 40(8), 913–929.
[doi:10.1111/ecog.02881](https://doi.org/10.1111/ecog.02881)

**\[8\]** Vehtari, A. and Lampinen, J. (2002). “Bayesian model
assessment and comparison using cross-validation predictive densities.”
*Neural Computation*, 14(10), 2439–2468.
[doi:10.1162/089976602760407935](https://doi.org/10.1162/089976602760407935)

**\[9\]** Yao, Y., Vehtari, A., Simpson, D., and Gelman, A. (2018).
“Using stacking to average Bayesian predictive distributions (with
discussion).” *Bayesian Analysis*, 13(3), 917–1007.
[doi:10.1214/17-BA1091](https://doi.org/10.1214/17-BA1091)
