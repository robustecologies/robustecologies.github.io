# Prior and posterior predictive checks

``` r

library(lucifer)
library(ggplot2)
```

## The predictive checking lifecycle

A fitted Bayesian model is not an end in itself but a generator of
predictions. Whether those predictions are plausible, given what we know
before and after observing data, determines whether the model deserves
our trust. Posterior predictive checking, introduced by Rubin
[\[1\]](#ref1) and formalized by Gelman, Meng, and Stern [\[2\]](#ref2),
asks a deceptively simple question: if the model were true, would it
produce data that look like the data we actually observed? Systematic
discrepancies between observed and replicated data reveal model
misspecification that convergence diagnostics alone cannot detect. A
chain may mix perfectly and still generate predictions that miss the
tails, undercount zeros, or ignore spatial structure.

The checking process has two stages. Prior predictive checks occur
before fitting: they sample parameters from the prior and push them
through the likelihood to generate synthetic datasets, testing whether
the joint model produces outcomes that are at least physically or
scientifically plausible [\[4\]](#ref4). If the prior predictive
distribution places substantial mass on absurd outcomes (negative
counts, temperatures above 1000 degrees, probabilities outside the unit
interval), the prior is too diffuse and should be tightened before any
data are conditioned upon. Posterior predictive checks occur after
fitting: they draw parameter values from the posterior distribution and
generate replicated datasets, comparing their properties against the
observed data. When the replicated data systematically deviate from the
observations, some aspect of the data-generating process has been
misspecified.

Different model families exhibit different failure modes, and no single
PPC diagnostic is universally informative. A density overlay reveals
distributional misfit (e.g., Gaussian errors applied to heavy-tailed
data) but says nothing about temporal autocorrelation. A Durbin-Watson
test detects residual autocorrelation but ignores spatial structure. A
rootogram diagnoses count data calibration but is meaningless for
continuous outcomes. This vignette works through five progressively
complex examples, each with synthetic data and known ground truth, to
demonstrate how the full suite of 39 PPC plot styles, 18 discrepancy
measures, and 6 convenience wrappers in lucifer addresses each model
class. The examples move from simple linear regression through
hierarchical models, count data, time series, and finally
spatio-temporal fields, and at each stage the PPC tools that are most
informative for that model family are introduced and interpreted.

## Theoretical foundations

### The posterior predictive distribution

The posterior predictive distribution integrates the likelihood over the
posterior uncertainty in the parameters. Given observed data \\y\\ and a
model parameterized by \\\theta\\, the distribution of a new or
replicated dataset \\\tilde{y}\\ is

\\p(\tilde{y} \mid y) = \int p(\tilde{y} \mid \theta)\\ p(\theta \mid
y)\\ d\theta.\\

Two interpretations are important. When \\\tilde{y}\\ is generated at
the same covariate values as the observed data, it is called a
replicated dataset \\y^{\text{rep}}\\ and the comparison against \\y\\
asks whether the model can reproduce the data that were used to fit it.
When \\\tilde{y}\\ is generated at new covariate values, it is a genuine
prediction \\y^{\text{new}}\\ and the comparison addresses
generalization. Both are instances of the same integral, but they answer
different questions: \\y^{\text{rep}}\\ diagnoses internal consistency;
\\y^{\text{new}}\\ diagnoses predictive adequacy.

Given \\S\\ posterior samples \\\theta^{(1)}, \ldots, \theta^{(S)}\\,
the posterior predictive distribution is approximated by the empirical
distribution of the \\S\\ replicated datasets: for each sample \\s\\,
generate \\\tilde{y}^{(s)} \sim p(\cdot \mid \theta^{(s)})\\ and collect
the results. This produces an \\N \times S\\ matrix of predictions,
where \\N\\ is the number of observations, that is the basis for all
graphical and numerical PPC diagnostics. The
[`predict()`](https://rdrr.io/r/stats/predict.html) method in lucifer
constructs this matrix automatically from any fitted object.

### Prior predictive checks

The prior predictive distribution marginalizes over the prior rather
than the posterior:

\\p(\tilde{y}) = \int p(\tilde{y} \mid \theta)\\ p(\theta)\\ d\theta.\\

Because this distribution depends only on the model specification and
the prior, not on the observed data, it can be computed before fitting.
The logic is simple: if the prior assigns high probability to parameter
values that produce absurd predictions, the model as a whole is poorly
specified regardless of what the data look like. The
[`prior_predictive_check()`](https://robustecologies.github.io/lucifer/reference/prior_predictive_check.md)
function in lucifer automates this process by sampling parameter vectors
from a user-supplied prior function `rprior(n)`, evaluating the Model at
each draw, and comparing the resulting predictive distribution against
the observed data.

Gabry et al. [\[4\]](#ref4) argue that prior predictive checking should
be a routine step in every Bayesian workflow, not an afterthought. A
common failure mode is priors that are technically proper but
scientifically absurd: a \\\text{Normal}(0, 10^6)\\ prior on a
regression coefficient in a model where \\y\\ ranges from 0 to 10 will
generate prior predictive datasets with values in the millions,
indicating that the prior is uninformative in a pathological rather than
a useful sense.

### Test statistics and discrepancy measures

A graphical comparison of replicated and observed data can be summarized
numerically through test statistics. Let \\T(\cdot)\\ be a scalar
function of the data. The posterior predictive p-value is

\\p_B = \Pr\bigl(T(y^{\text{rep}}) \geq T(y) \mid y\bigr) \approx
\frac{1}{S} \sum\_{s=1}^{S} \mathbf{1}\bigl\[T(\tilde{y}^{(s)}) \geq
T(y)\bigr\],\\

where values near 0 or 1 indicate that the observed test statistic is
extreme relative to the model’s predictions [\[2\]](#ref2),
[\[3\]](#ref3). The choice of \\T\\ depends on the aspect of the data
one wishes to probe: the mean detects location bias, the standard
deviation detects variance misspecification, the maximum detects tail
inadequacy, and the proportion of zeros detects zero-inflation in count
data.

The [`summary()`](https://rdrr.io/r/base/summary.html) method for PPC
objects in lucifer computes 18 discrepancy measures. For continuous
outcomes these include Chi-Square (Pearson residuals), Chi-Square2
(replicated vs expected), Kurtosis and Skewness (distributional shape),
L-criterion (Laud and Ibrahim predictive loss [\[8\]](#ref8)), MASE
(mean absolute scaled error, standard in time series), MSE and RMSE
(prediction error), PPL (posterior predictive loss), Quadratic Loss and
Quadratic Utility (decision-theoretic), and several logical tests that
probe whether specific summary statistics of the replicated data exceed
those of the observed data. For categorical outcomes, the
misclassification rate \\p(\hat{y}\_i \neq y_i)\\ is available.

### Graphical calibration

Graphical PPC tools fall into several families, each suited to different
diagnostic questions. Density and histogram overlays show whether the
marginal distribution of replicated data matches the observed
distribution; systematic departures in the tails, modes, or skewness are
immediately visible. ECDF overlays provide a cumulative view of the same
information, with the advantage that departures are easier to quantify
visually. Interval plots display the posterior predictive credible
interval for each observation, with the observed value superimposed;
observations that fall outside their intervals are potential outliers or
regions of model misfit. Fan chart ribbons extend this idea to ordered
data (time series, spatial transects) by showing nested credible bands
(50%, 80%, 95%) as a function of the ordering variable.

For count data, the hanging rootogram [\[10\]](#ref10) is the diagnostic
of choice: it compares observed and expected frequencies on a
square-root scale, with bars hanging from the expected curve so that
departures from zero indicate miscalibration at specific counts. For
temporal data, the Durbin-Watson statistic [\[11\]](#ref11) detects
residual autocorrelation that indicates a missing lag or structural
component. The Jarque-Bera test [\[12\]](#ref12) diagnoses non-normality
in residuals, and Mardia’s test [\[13\]](#ref13) extends this to the
multivariate case.

The LOO-PIT (leave-one-out probability integral transform) provides a
unified calibration diagnostic applicable to any model [\[9\]](#ref9).
For each observation \\i\\, the PIT value is \\u_i = P(y_i^{\text{rep}}
\leq y_i \mid y\_{-i})\\, approximated via PSIS-LOO weights. If the
model is well-calibrated, the PIT values are uniformly distributed, and
their empirical CDF follows the diagonal. Departures from uniformity
indicate specific miscalibration patterns: a U-shaped histogram
indicates underdispersion, an inverted-U indicates overdispersion, and
asymmetry indicates location bias.

### Quantitative summaries

Three scalar summaries complement the graphical diagnostics. The
Bayesian Predictive Information Criterion (BPIC) [\[7\]](#ref7) is
defined as \\\text{BPIC} = \bar{D} + 2 p_D\\, where \\\bar{D}\\ is the
posterior mean deviance and \\p_D = \text{Var}(D)/2\\ is the effective
number of parameters; lower BPIC indicates better predictive fit. The
L-criterion [\[8\]](#ref8) is \\L = \sum_i
\sqrt{\text{Var}(\tilde{y}\_i) + (y_i - E\[\tilde{y}\_i\])^2}\\,
combining posterior predictive variance (model uncertainty) and bias
(systematic misfit) into a single measure; models with lower \\L\\ are
preferred. Concordance measures the fraction of observations whose
predictive quantile \\q_i = P(\tilde{y}\_i \geq y_i)\\ falls within the
interval \\(0.025, 0.975)\\; a well-calibrated model should achieve
concordance near 0.95.

## The lucifer PPC toolkit

### Generating predictions

Posterior predictive samples are generated by the
[`predict()`](https://rdrr.io/r/stats/predict.html) method, which is
available for all eight inference families. The method iterates through
posterior parameter samples, evaluates the Model function at each sample
to obtain `yhat`, and collects the results into an \\N \times S\\
matrix.

``` r

# Posterior predictive samples from any fitted object
ppc <- predict(fit, Model, Data)
```

Prior predictive samples are generated by
[`prior_predictive_check()`](https://robustecologies.github.io/lucifer/reference/prior_predictive_check.md),
which takes a user-supplied sampling function `rprior(n)` that returns
an \\n \times K\\ matrix of prior draws.

``` r

# Prior predictive check with density overlay
prior_predictive_check(Model, Data, rprior = rprior, n = 500, type = "density")

# Prior predictive with interval bands
prior_predictive_check(Model, Data, rprior = rprior, n = 500, type = "intervals")

# Prior predictive fan chart
prior_predictive_check(Model, Data, rprior = rprior, n = 500, type = "ribbon")

# Prior predictive test statistic
prior_predictive_check(Model, Data, rprior = rprior, n = 500,
                       type = "stat", stat_fun = sd)
```

Two additional functions compare prior and posterior distributions. The
[`prior_vs_posterior()`](https://robustecologies.github.io/lucifer/reference/prior_vs_posterior.md)
function overlays the prior and posterior densities for selected
parameters, with optional ground truth reference lines. The
[`prior_sensitivity()`](https://robustecologies.github.io/lucifer/reference/prior_sensitivity.md)
function overlays posteriors obtained under different prior
specifications to visualize how sensitive the inference is to the choice
of prior.

``` r

# Prior vs posterior with ground truth
prior_vs_posterior(fit, prior = rprior,
                   Parms = c("beta", "sigma"),
                   ground_truth = c(beta = 2, sigma = 1.5))

# Prior sensitivity across three specifications
prior_sensitivity(
    list("N(0,1)" = fit_narrow, "N(0,10)" = fit_default, "N(0,100)" = fit_vague),
    ground_truth = c(beta = 2)
)
```

### Plotting predictions: 39 styles

The [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method for
PPC objects accepts a `Style` argument that selects from 39
visualization types. The following table groups them by diagnostic
family.

``` r

styles <- data.frame(
    Family = c(
        rep("Distribution comparison", 5),
        rep("Fitted vs observed", 3),
        rep("Intervals and bands", 3),
        rep("Residual diagnostics", 4),
        rep("Statistical tests", 5),
        rep("Count data", 2),
        rep("Categorical and grouped", 3),
        rep("Spatial and temporal", 5),
        rep("Calibration", 1),
        rep("Multivariate variants", 8)
    ),
    Style = c(
        "Density", "Density Overlay", "Histogram Overlay",
        "ECDF Overlay", "ECDF",
        "Fitted", "Scatter Average", "Bars",
        "Intervals", "Ribbon", "Predictive Quantiles",
        "Residuals", "Residual Density", "Error Scatter", "Error Histogram",
        "DW", "Jarque-Bera", "Mardia", "Stat", "Stat 2D",
        "Rootogram", "Bars",
        "Covariates", "Covariates, Categorical DV", "Violin Grouped",
        "Spatial", "Spatial Uncertainty", "Space-Time by Space",
        "Space-Time by Time", "Time-Series",
        "LOO-PIT",
        "Fitted, Multivariate, C", "Fitted, Multivariate, R",
        "Residuals, Multivariate, C", "Residuals, Multivariate, R",
        "Residual Density, Multivariate, C", "Residual Density, Multivariate, R",
        "Time-Series, Multivariate, C", "Time-Series, Multivariate, R"
    ),
    Requirements = c(
        "None", "None", "None", "None", "None",
        "None", "None", "None",
        "None", "None (forecast_start optional)", "None",
        "None", "None", "None", "None",
        "None", "None", "Data$Y matrix", "stat_fun", "stat_fun, stat_fun2",
        "None (count yhat)", "None",
        "Data$X", "Data$X", "Group vector",
        "Data$longitude, Data$latitude",
        "Data$longitude, Data$latitude",
        "Data$longitude, Data$latitude, Data$S, Data$T",
        "Data$longitude, Data$latitude, Data$S, Data$T",
        "None",
        "loo object",
        "Data$Y", "Data$Y", "Data$Y", "Data$Y",
        "Data$Y", "Data$Y", "Data$Y", "Data$Y"
    )
)

knitr::kable(styles, caption = "PPC plot styles available in lucifer")
```

### Summarizing predictions

The [`summary()`](https://rdrr.io/r/base/summary.html) method computes
per-observation predictive summaries and aggregate diagnostics. The
output includes the posterior predictive mean, standard deviation, 95%
credible interval, predictive quantile, and an optional discrepancy
statistic for each observation.

``` r

# Basic summary with concordance and BPIC
s <- summary(ppc)

# Summary with a specific discrepancy measure
s <- summary(ppc, Discrep = "Chi-Square")
s$BPIC          # Bayesian Predictive Information Criterion
s$Concordance   # Fraction of y within 95% posterior predictive interval
s$L.criterion   # L-criterion (Laud and Ibrahim, 1995)
s$Summary       # Per-observation table: y, Mean, SD, LB, Median, UB, PQ, Discrep
```

### Convenience wrappers

Six convenience functions provide a bayesplot-compatible calling
convention for the most commonly used PPC plot styles. They accept a PPC
object and delegate to the internal plotting engine with the appropriate
`Style` argument.

``` r

ppc_dens_overlay(ppc)                                     # Density Overlay
ppc_intervals(ppc)                                        # Intervals
ppc_ribbon(ppc, forecast_start = 201)                     # Ribbon (fan chart)
ppc_rootogram(ppc)                                        # Rootogram
ppc_stat(ppc, stat_fun = function(x) mean(x > 0))        # Test statistic
ppc_loo_pit(ppc, loo = loo_obj)                           # LOO-PIT
```

### Converting validation objects

When out-of-sample validation has been performed, the
[`as.ppc()`](https://robustecologies.github.io/lucifer/reference/as.ppc.md)
function converts a `demonoid.val` object to a PPC object for plotting.
The `set` argument controls which data are included: modeled data only
(`set = 1`), validation data only (`set = 2`), or both (`set = 3`, the
default).

``` r

ppc_val <- as.ppc(validation_result, set = 2)
plot(ppc_val, Style = "Fitted")
```

## Worked examples

### Example 1: normal linear regression

This first example introduces the full PPC lifecycle on the simplest
possible model: a Gaussian linear regression with known ground truth. It
demonstrates all four prior predictive plot types, 16 posterior
predictive styles, and the detection of distributional misspecification
when the error distribution is wrong.

#### Model specification

The model is a standard normal linear regression with three predictors,
weakly informative priors on the coefficients, and a half-Cauchy prior
on the residual standard deviation:

\\y_i \sim \mathcal{N}(\mu_i, \sigma^2), \quad \mu_i = \mathbf{X}\_i
\boldsymbol{\beta}, \quad i = 1, \ldots, N\\ \\\beta_j \sim
\mathcal{N}(0, 100), \quad j = 1, \ldots, 3\\ \\\sigma \sim
\text{Half-Cauchy}(25)\\

``` r

# Ground truth
set.seed(42)
N <- 200
beta_true <- c(2.0, -1.0, 0.5)
sigma_true <- 1.5

# Design matrix: intercept + 2 predictors
X <- cbind(1, matrix(rnorm(N * 2), N, 2))
mu_true <- X %*% beta_true
y <- rnorm(N, mu_true, sigma_true)
```

The
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL produces a ready-to-use Model function, data template builder, and
initial value generator from a declarative specification.

``` r

spec <- model_spec("
    y ~ Normal(mu, sigma)
    mu = X %*% beta
    beta[j] ~ Normal(0, 100), j = 1,...,J
    sigma ~ HalfCauchy(25)
")

Data <- spec$data_template(y = y, X = X, J = 3)
Initial.Values <- spec$initial_values(Data)
```

#### Prior predictive checks

Before fitting, we verify that the prior produces predictions that are
scientifically plausible. The `rprior` function draws from the joint
prior, returning an \\n \times 4\\ matrix (three coefficients plus one
scale parameter).

``` r

rprior <- function(n) {
    cbind(
        rnorm(n, 0, 10),           # beta[1]
        rnorm(n, 0, 10),           # beta[2]
        rnorm(n, 0, 10),           # beta[3]
        abs(rcauchy(n, 0, 25))     # sigma
    )
}
```

The four prior predictive plot types each probe a different aspect of
the prior specification. The density overlay shows the marginal
distribution of prior predictive \\\tilde{y}\\ against the observed
\\y\\: if the observed data sit comfortably within the prior predictive
envelope, the prior is compatible with the data without being
excessively concentrated. The interval plot shows observation-level
prior predictive credible bands, revealing whether the prior allows
sufficient variation at each data point. The ribbon plot extends this to
an ordered sequence, useful for detecting whether the prior allows the
regression to tilt in implausible directions. The test statistic
compares a scalar summary of the prior predictive data against the
observed value.

``` r

# Density overlay: prior predictive vs observed
prior_predictive_check(spec$Model, Data, rprior = rprior, n = 300,
                       type = "density")
```

``` r

# Interval plot: observation-level prior predictive bands
prior_predictive_check(spec$Model, Data, rprior = rprior, n = 300,
                       type = "intervals")
```

``` r

# Fan chart: prior predictive bands as ribbon
prior_predictive_check(spec$Model, Data, rprior = rprior, n = 300,
                       type = "ribbon")
```

``` r

# Test statistic: prior predictive standard deviation vs observed
prior_predictive_check(spec$Model, Data, rprior = rprior, n = 300,
                       type = "stat", stat_fun = sd)
```

With \\\text{Normal}(0, 100)\\ priors on the coefficients and a
\\\text{Half-Cauchy}(25)\\ prior on \\\sigma\\, the prior predictive
distribution is wide enough to cover the data but does not extend to
absurd values, confirming that the priors are weakly informative rather
than pathologically vague.

#### Fitting and convergence

``` r

rx <- Prescribe(spec$Model, Data, Initial.Values)
print(rx)
```

``` r

fit <- lucifer(spec$Model, Data, Initial.Values,
               Iterations = 10000, Status = 2000,
               Thinning = 5, Algorithm = "NUTS")
```

``` r

cx <- Consort(fit)
print(cx)
```

#### Prior vs posterior comparison

Before proceeding to predictive checks, a direct comparison of the prior
and posterior distributions shows how much the data updated each
parameter. The `ground_truth` argument adds vertical reference lines at
the known true values.

``` r

truth <- c("beta[1]" = 2.0, "beta[2]" = -1.0,
           "beta[3]" = 0.5, sigma = 1.5)

prior_vs_posterior(
    fit,
    prior = rprior,
    ground_truth = truth
)
```

For all four parameters the posterior concentrates tightly relative to
the prior, and the true values fall within the posterior mass. With \\N
= 200\\ observations, the data dominate the inference, as expected for a
well-identified linear model.

#### Posterior predictive checks: core diagnostics

The [`predict()`](https://rdrr.io/r/stats/predict.html) method generates
the \\N \times S\\ matrix of posterior predictive samples.

``` r

ppc <- predict(fit, Model = spec$Model, Data = Data)
```

The default plot style is “Density Overlay”, which superimposes density
estimates of replicated datasets (one per posterior sample) against the
observed data density. For a well-specified model, the observed density
should be a typical member of the replicated ensemble.

``` r

plot(ppc, Style = "Density")
```

``` r

# Convenience wrapper (equivalent to Style = "Density Overlay")
ppc_dens_overlay(ppc)
```

``` r

plot(ppc, Style = "Histogram Overlay")
```

``` r

plot(ppc, Style = "ECDF Overlay")
```

The ECDF overlay provides a cumulative view of the same comparison.
Departures between the observed and replicated ECDFs indicate
miscalibration at specific quantiles.

``` r

plot(ppc, Style = "ECDF")
```

The fitted vs observed plot shows the posterior predictive mean for each
observation against its observed value, with 95% credible intervals.
Points falling on the diagonal indicate accurate prediction.

``` r

plot(ppc, Style = "Fitted", Data = Data)
```

``` r

plot(ppc, Style = "Scatter Average")
```

``` r

# Convenience wrapper (equivalent to Style = "Intervals")
ppc_intervals(ppc)
```

The predictive quantile plot shows \\q_i = P(\tilde{y}\_i \geq y_i)\\
for each observation. Under a well-calibrated model, these quantiles
should be approximately uniformly distributed; extreme values (near 0
or 1) flag potential outliers.

``` r

plot(ppc, Style = "Predictive Quantiles")
```

Residual diagnostics examine the distribution of \\y_i -
E\[\tilde{y}\_i\]\\ and probe for patterns that indicate model misfit.

``` r

plot(ppc, Style = "Residuals")
```

``` r

plot(ppc, Style = "Residual Density")
```

``` r

plot(ppc, Style = "Error Scatter")
```

``` r

plot(ppc, Style = "Error Histogram")
```

The covariates plot shows the posterior predictive mean as a function of
each covariate column in the design matrix, overlaid with a loess
smooth. Systematic curvature in the smooth indicates that the linear
specification is inadequate for that predictor.

``` r

plot(ppc, Style = "Covariates", Data = Data)
```

Test statistics provide a scalar summary of model-data agreement. The
“Stat” style computes \\T(y^{\text{rep}})\\ for each replicated dataset
and displays the resulting histogram with the observed \\T(y)\\ as a
reference line.

``` r

plot(ppc, Style = "Stat", stat_fun = mean)
```

The “Stat 2D” style plots two test statistics jointly, showing whether
the observed pair \\(T_1(y), T_2(y))\\ falls within the cloud of
replicated pairs.

``` r

plot(ppc, Style = "Stat 2D", stat_fun = mean, stat_fun2 = sd)
```

#### Quantitative summaries

The [`summary()`](https://rdrr.io/r/base/summary.html) method computes
BPIC, concordance, L-criterion, and per-observation summaries. Several
discrepancy measures probe different aspects of the model-data
relationship.

``` r

s1 <- summary(ppc, Discrep = "Chi-Square")
```

``` r

s2 <- summary(ppc, Discrep = "MSE")
```

``` r

s3 <- summary(ppc, Discrep = "RMSE")
```

``` r

s4 <- summary(ppc, Discrep = "L.criterion")
```

``` r

s5 <- summary(ppc, Discrep = "mean(yhat[i,]) > mean(y)")
```

``` r

s6 <- summary(ppc, Discrep = "sd(yhat[i,]) > sd(y)")
```

``` r

s7 <- summary(ppc, Discrep = "min(yhat[i,]) < min(y)")
```

The concordance should be close to 0.95 for a well-calibrated model.
Values substantially below 0.95 indicate that too many observations fall
outside their predictive intervals (the model is overconfident or
misspecified), while values above 0.95 indicate unnecessarily wide
intervals (the model is conservative).

``` r

diag_table <- data.frame(
    Metric = c("BPIC", "Concordance", "L-criterion"),
    Value = c(s1$BPIC[1, "BPIC"], s1$Concordance, s4$L.criterion)
)
knitr::kable(diag_table, digits = 3,
             caption = "Quantitative PPC summaries for the well-specified model")
```

#### Diagnosing misspecification: Gaussian vs Student-t errors

To demonstrate how PPC detects distributional misspecification, we
generate a second dataset from Student-t errors with 4 degrees of
freedom (heavy tails) but fit the same Gaussian model.

``` r

set.seed(123)
y_heavy <- mu_true + sigma_true * rt(N, df = 4) / sqrt(4 / (4 - 2))
Data_heavy <- spec$data_template(y = y_heavy, X = X, J = 3)
```

``` r

fit_heavy <- lucifer(spec$Model, Data_heavy, Initial.Values,
                     Iterations = 10000, Status = 2000,
                     Thinning = 5, Algorithm = "NUTS")
```

``` r

ppc_heavy <- predict(fit_heavy, Model = spec$Model, Data = Data_heavy)
```

The density overlay immediately reveals the misfit: the Gaussian model
produces replicated datasets that are too concentrated near the center,
failing to reproduce the heavy tails of the observed data.

``` r

# Well-specified model: Gaussian data, Gaussian model
ppc_dens_overlay(ppc)

# Misspecified model: heavy-tailed data, Gaussian model
ppc_dens_overlay(ppc_heavy)
```

The quantitative comparison confirms what the eye detects.

``` r

s_good <- summary(ppc, Discrep = "Chi-Square")
s_bad <- summary(ppc_heavy, Discrep = "Chi-Square")

comparison <- data.frame(
    Model = c("Well-specified", "Misspecified"),
    BPIC = c(s_good$BPIC[1, "BPIC"], s_bad$BPIC[1, "BPIC"]),
    Concordance = c(s_good$Concordance, s_bad$Concordance),
    Chi_Square = c(s_good$Discrepancy.Statistic, s_bad$Discrepancy.Statistic)
)
knitr::kable(comparison, digits = 3,
             caption = "PPC comparison: well-specified vs misspecified error distribution")
```

The misspecified model shows lower concordance (observations falling
outside their predictive intervals more often than expected) and a
higher Chi-Square discrepancy statistic. The BPIC difference quantifies
the predictive cost of using the wrong error distribution.

### Example 2: hierarchical eight schools

Hierarchical models produce posterior predictive distributions that
reflect partial pooling across groups. In this example, PPC reveals how
shrinkage toward the group mean affects predictions, and the
multivariate plot variants decompose discrepancies by school and by
variable. Prior sensitivity analysis shows how the choice of hyperprior
on the between-school standard deviation \\\tau\\ affects the posterior
and, consequently, the predictive distribution.

#### Model specification

The eight schools model [\[2\]](#ref2) is a canonical example of
hierarchical shrinkage:

\\y_j \sim \mathcal{N}(\theta_j, \sigma_j^2), \quad j = 1, \ldots, 8\\
\\\theta_j \sim \mathcal{N}(\mu, \tau^2)\\ \\\mu \sim \mathcal{N}(0,
100), \quad \tau \sim \text{Half-Cauchy}(25)\\

where \\y_j\\ and \\\sigma_j\\ are the observed treatment effect and
known standard error in school \\j\\.

``` r

# Ground truth: hierarchical structure
set.seed(42)
mu_true <- 8.0
tau_true <- 6.0
theta_true <- rnorm(8, mu_true, tau_true)

# Known standard errors (from the original eight schools data)
sigma_j <- c(15, 10, 16, 11, 9, 11, 10, 18)

# Observed treatment effects
y_schools <- rnorm(8, theta_true, sigma_j)
school_names <- paste("School", LETTERS[1:8])
```

``` r

Model_schools <- function(parm, Data) {
    mu <- parm[1]
    tau <- interval(parm[2], 1e-100, Inf)
    parm[2] <- tau
    theta <- parm[3:10]

    # Log-prior
    theta.prior <- sum(dnorm(theta, mu, tau, log = TRUE))
    mu.prior <- dnormv(mu, 0, 100, log = TRUE)
    tau.prior <- dhalfcauchy(tau, 25, log = TRUE)

    # Log-likelihood
    LL <- sum(dnorm(Data$y, theta, Data$sigma_j, log = TRUE))

    LP <- LL + theta.prior + mu.prior + tau.prior

    # Fitted values
    yhat <- rnorm(Data$J, theta, Data$sigma_j)

    Monitor <- LP
    return(list(LP = LP, Dev = -2 * LL, Monitor = Monitor,
                yhat = yhat, parm = parm))
}

Data_schools <- list(
    y = y_schools,
    Y = matrix(y_schools, nrow = 8, ncol = 1),
    sigma_j = sigma_j,
    J = 8,
    mon.names = "LP",
    parm.names = c("mu", "tau", paste0("theta[", 1:8, "]"))
)

IV_schools <- c(0, 5, y_schools)
```

#### Prior predictive checks

``` r

rprior_schools <- function(n) {
    mu <- rnorm(n, 0, 10)
    tau <- abs(rcauchy(n, 0, 25))
    theta <- matrix(NA, n, 8)
    for (j in 1:8) theta[, j] <- rnorm(n, mu, tau)
    cbind(mu, tau, theta)
}
```

``` r

prior_predictive_check(Model_schools, Data_schools,
                       rprior = rprior_schools, n = 300, type = "density")
```

The test statistic \\T(y) = \max(y)\\ checks whether the prior allows
treatment effects as large as the observed maximum. If the prior
predictive maximum rarely exceeds the observed maximum, the prior is too
tight.

``` r

prior_predictive_check(Model_schools, Data_schools,
                       rprior = rprior_schools, n = 300,
                       type = "stat", stat_fun = max)
```

#### Fitting and convergence

``` r

fit_schools <- lucifer(Model_schools, Data_schools, IV_schools,
                       Iterations = 20000, Status = 5000,
                       Thinning = 10, Algorithm = "NUTS")

cx_schools <- Consort(fit_schools)
print(cx_schools)
```

#### Prior vs posterior with ground truth

``` r

truth_schools <- c(mu = mu_true, tau = tau_true,
                   setNames(theta_true, paste0("theta[", 1:8, "]")))

prior_vs_posterior(
    fit_schools,
    prior = rprior_schools,
    Parms = c("mu", "tau"),
    ground_truth = truth_schools
)
```

#### Posterior predictive checks: shrinkage and group structure

``` r

ppc_schools <- predict(fit_schools, Model = Model_schools, Data = Data_schools)
```

The interval plot makes the hierarchical shrinkage pattern immediately
visible: the posterior predictive intervals for all eight schools are
more similar in width and location than the raw data would suggest,
because the hierarchical model borrows strength across schools.

``` r

ppc_intervals(ppc_schools)
```

The Violin Grouped style shows the full posterior predictive
distribution for each school as a violin plot, grouped by school label.

``` r

plot(ppc_schools, Style = "Violin Grouped",
     Group = school_names)
```

``` r

plot(ppc_schools, Style = "Bars")
```

``` r

plot(ppc_schools, Style = "ECDF")
```

Test statistics summarize the overall agreement. The “Stat” style with
\\T = \text{mean}\\ shows whether the model reproduces the grand mean of
the treatment effects.

``` r

ppc_stat(ppc_schools, stat_fun = mean)
```

``` r

plot(ppc_schools, Style = "Stat 2D", stat_fun = mean, stat_fun2 = sd)
```

The Jarque-Bera test probes whether the residuals are normally
distributed, as the model assumes.

``` r

plot(ppc_schools, Style = "Jarque-Bera")
```

#### Multivariate diagnostics

The multivariate PPC styles decompose the model-data comparison by
column (variable) and by row (observation/school). For the eight schools
model, the “Multivariate, R” variants produce one panel per school,
showing how well the model predicts each school individually.

``` r

plot(ppc_schools, Style = "Fitted, Multivariate, C", Data = Data_schools)
```

``` r

plot(ppc_schools, Style = "Fitted, Multivariate, R", Data = Data_schools)
```

``` r

plot(ppc_schools, Style = "Residuals, Multivariate, C", Data = Data_schools)
```

``` r

plot(ppc_schools, Style = "Residuals, Multivariate, R", Data = Data_schools)
```

``` r

plot(ppc_schools, Style = "Residual Density, Multivariate, C",
     Data = Data_schools)
```

``` r

plot(ppc_schools, Style = "Residual Density, Multivariate, R",
     Data = Data_schools)
```

#### Discrepancy measures

``` r

s_schools1 <- summary(ppc_schools, Discrep = "Chi-Square2")
s_schools2 <- summary(ppc_schools, Discrep = "Kurtosis")
s_schools3 <- summary(ppc_schools, Discrep = "Skewness")
```

#### Prior sensitivity

Prior sensitivity analysis shows how the choice of hyperprior on
\\\tau\\ affects the posterior. Three specifications are compared: a
tight half-Cauchy(5), the default half-Cauchy(25), and a diffuse
half-Cauchy(100).

``` r

# Tight prior on tau
Model_tight <- Model_schools  # same model, different prior
Data_tight <- Data_schools
# (In practice, modify the tau.prior line; here we show the workflow)

fit_tight <- lucifer(Model_schools, Data_schools, IV_schools,
                     Iterations = 20000, Status = 5000,
                     Thinning = 10, Algorithm = "NUTS")

# Repeat with modified tau prior (HC(5) and HC(100))
# fit_narrow_tau, fit_wide_tau obtained similarly
```

``` r

prior_sensitivity(
    list("HC(5)" = fit_tight,
         "HC(25)" = fit_schools,
         "HC(100)" = fit_schools),  # placeholder; use actual fits
    Parms = c("mu", "tau"),
    ground_truth = c(mu = mu_true, tau = tau_true)
)
```

The prior sensitivity plot reveals whether the posterior is dominated by
the data or by the prior. For the eight schools problem, the posterior
on \\\tau\\ is notoriously sensitive to the hyperprior because only 8
groups contribute information about the between-group variance; the
posterior on \\\mu\\ is more robust because all observations contribute
to the grand mean estimate.

### Example 3: count data (Poisson vs zero-inflated Poisson)

Count data require special PPC diagnostics because the discreteness of
the outcome creates failure modes that continuous diagnostics miss. The
most informative tool for count data is the hanging rootogram
[\[10\]](#ref10), which compares observed and expected frequencies on a
square-root scale. This example simulates zero-inflated count data, fits
both a standard Poisson model (misspecified) and a zero-inflated Poisson
model (correctly specified), and shows how PPC reveals the excess zeros
that the Poisson model cannot accommodate.

#### Model specification

The data are generated from a zero-inflated Poisson (ZIP) process where
30% of observations are structural zeros.

``` r

set.seed(42)
N <- 300
J <- 3
beta_true <- c(1.0, 0.3, -0.5)
gamma_true <- c(-0.8, 0.5)  # zero-inflation logistic coefficients
pi_true <- 0.3               # average zero-inflation probability

X <- cbind(1, matrix(rnorm(N * (J - 1)), N, J - 1))
Z <- cbind(1, rnorm(N))

lambda_true <- exp(X %*% beta_true)
pi_i <- plogis(Z %*% gamma_true)

# Generate zero-inflated Poisson data
zero_indicator <- rbinom(N, 1, pi_i)
y_count <- ifelse(zero_indicator == 1, 0L,
                  rpois(N, lambda_true))
```

The Poisson model ignores zero-inflation:

\\y_i \sim \text{Poisson}(\lambda_i), \quad \log(\lambda_i) =
\mathbf{X}\_i \boldsymbol{\beta}\\ \\\beta_j \sim \mathcal{N}(0, 100),
\quad j = 1, \ldots, J\\

``` r

spec_pois <- model_spec("
    y ~ Poisson(lambda)
    lambda = exp(X %*% beta)
    beta[j] ~ Normal(0, 100), j = 1,...,J
")

Data_pois <- spec_pois$data_template(y = y_count, X = X, J = J)
IV_pois <- spec_pois$initial_values(Data_pois)
```

The ZIP model accounts for the structural zeros:

\\P(y_i = 0) = \pi_i + (1 - \pi_i) e^{-\lambda_i}\\ \\P(y_i = k) = (1 -
\pi_i) \frac{\lambda_i^k e^{-\lambda_i}}{k!}, \quad k \> 0\\
\\\text{logit}(\pi_i) = \mathbf{Z}\_i \boldsymbol{\gamma}\\
\\\log(\lambda_i) = \mathbf{X}\_i \boldsymbol{\beta}\\

``` r

Model_zip <- function(parm, Data) {
    beta <- parm[Data$pos.beta]
    gamma <- parm[Data$pos.gamma]

    lambda <- as.vector(exp(Data$X %*% beta))
    pi_i <- as.vector(plogis(Data$Z %*% gamma))

    # ZIP log-likelihood
    idx0 <- which(Data$y == 0)
    idx1 <- which(Data$y > 0)
    LL <- 0
    if (length(idx0) > 0)
        LL <- LL + sum(log(pi_i[idx0] +
                           (1 - pi_i[idx0]) * dpois(0, lambda[idx0])))
    if (length(idx1) > 0)
        LL <- LL + sum(log((1 - pi_i[idx1]) *
                           dpois(Data$y[idx1], lambda[idx1])))

    # Log-prior
    beta.prior <- sum(dnormv(beta, 0, 100, log = TRUE))
    gamma.prior <- sum(dnormv(gamma, 0, 100, log = TRUE))
    LP <- LL + beta.prior + gamma.prior

    # Fitted values: simulate from ZIP
    zi <- rbinom(Data$N, 1, pi_i)
    yhat <- ifelse(zi == 1, 0L, rpois(Data$N, lambda))

    Monitor <- LP
    return(list(LP = LP, Dev = -2 * LL, Monitor = Monitor,
                yhat = yhat, parm = parm))
}

Data_zip <- list(
    y = y_count, X = X, Z = Z, N = N,
    pos.beta = 1:J, pos.gamma = (J + 1):(J + 2),
    mon.names = "LP",
    parm.names = c(paste0("beta[", 1:J, "]"),
                   paste0("gamma[", 1:2, "]"))
)

IV_zip <- c(rep(0, J), rep(0, 2))
```

#### Prior predictive checks

The most informative prior predictive test statistic for count data is
the proportion of zeros. This checks whether the prior allows
zero-inflation rates compatible with the observed data.

``` r

rprior_pois <- function(n) {
    matrix(rnorm(n * J, 0, 10), n, J)
}
```

``` r

prior_predictive_check(spec_pois$Model, Data_pois,
                       rprior = rprior_pois, n = 300, type = "density")
```

``` r

prior_predictive_check(spec_pois$Model, Data_pois,
                       rprior = rprior_pois, n = 300,
                       type = "stat",
                       stat_fun = function(x) mean(x == 0))
```

#### Fitting and convergence

``` r

fit_pois <- lucifer(spec_pois$Model, Data_pois, IV_pois,
                    Iterations = 10000, Status = 2000,
                    Thinning = 5, Algorithm = "NUTS")
```

``` r

fit_zip <- lucifer(Model_zip, Data_zip, IV_zip,
                   Iterations = 10000, Status = 2000,
                   Thinning = 5, Algorithm = "NUTS")
```

``` r

Consort(fit_pois)
Consort(fit_zip)
```

#### Prior vs posterior with ground truth

``` r

truth_zip <- c(setNames(beta_true, paste0("beta[", 1:J, "]")),
               setNames(gamma_true, paste0("gamma[", 1:2, "]")))

rprior_zip <- function(n) {
    cbind(matrix(rnorm(n * J, 0, 10), n, J),
          matrix(rnorm(n * 2, 0, 10), n, 2))
}

prior_vs_posterior(fit_zip, prior = rprior_zip, ground_truth = truth_zip)
```

#### Posterior predictive checks: discrete distributions

``` r

ppc_pois <- predict(fit_pois, Model = spec_pois$Model, Data = Data_pois)
ppc_zip <- predict(fit_zip, Model = Model_zip, Data = Data_zip)
```

The rootogram is the definitive count data PPC diagnostic. For the
Poisson model, the hanging bars at count zero will be noticeably above
zero, indicating that the model underpredicts the number of zeros. For
the ZIP model, the rootogram should be well-calibrated across all
counts.

``` r

# Poisson: visible excess at zero
ppc_rootogram(ppc_pois)
```

``` r

# ZIP: well-calibrated
ppc_rootogram(ppc_zip)
```

``` r

plot(ppc_pois, Style = "Bars")
plot(ppc_zip, Style = "Bars")
```

``` r

ppc_dens_overlay(ppc_pois)
ppc_dens_overlay(ppc_zip)
```

The test statistic \\T(y) = \text{proportion of zeros}\\ directly
targets the zero-inflation mechanism.

``` r

ppc_stat(ppc_pois, stat_fun = function(x) mean(x == 0))
ppc_stat(ppc_zip, stat_fun = function(x) mean(x == 0))
```

The maximum count test statistic checks whether the model can reproduce
the tail of the count distribution.

``` r

plot(ppc_pois, Style = "Stat", stat_fun = max)
plot(ppc_zip, Style = "Stat", stat_fun = max)
```

``` r

plot(ppc_pois, Style = "Covariates, Categorical DV", Data = Data_pois)
```

``` r

plot(ppc_pois, Style = "Scatter Average")
plot(ppc_zip, Style = "Scatter Average")
```

``` r

plot(ppc_pois, Style = "Predictive Quantiles")
```

#### Discrepancy measures and model comparison

``` r

s_pois <- summary(ppc_pois, Discrep = "Chi-Square")
s_zip <- summary(ppc_zip, Discrep = "Chi-Square")

s_pois_d <- summary(ppc_pois, Discrep = "mean(yhat[i,] > d)", d = 0)
s_zip_d <- summary(ppc_zip, Discrep = "mean(yhat[i,] > d)", d = 0)

s_pois_r <- summary(ppc_pois, Discrep = "round(yhat[i,]) = d", d = 0)
s_zip_r <- summary(ppc_zip, Discrep = "round(yhat[i,]) = d", d = 0)

s_pois_c <- summary(ppc_pois, Categorical = TRUE,
                    Discrep = "p(yhat[i,] != y[i])")
```

``` r

model_comp <- data.frame(
    Model = c("Poisson", "ZIP"),
    BPIC = c(s_pois$BPIC[1, "BPIC"], s_zip$BPIC[1, "BPIC"]),
    Concordance = c(s_pois$Concordance, s_zip$Concordance),
    Chi_Square = c(s_pois$Discrepancy.Statistic, s_zip$Discrepancy.Statistic)
)
knitr::kable(model_comp, digits = 3,
             caption = "PPC comparison: Poisson vs ZIP for zero-inflated data")
```

The formal model comparison via Arena confirms what the rootogram
reveals visually: the ZIP model has lower BPIC, higher concordance, and
lower Chi-Square discrepancy.

``` r

arena_counts <- Arena(x = list(Poisson = fit_pois, ZIP = fit_zip))
print(arena_counts)
plot(arena_counts, type = "efficiency")
```

### Example 4: autoregressive time series

Time series data require PPC diagnostics that probe temporal structure:
residual autocorrelation, forecast calibration, and the ability to
reproduce the observed dynamics. This example generates an AR(2) process
with known parameters, fits both AR(1) (misspecified) and AR(2)
(correctly specified) models, and shows how the Durbin-Watson test and
ribbon plots detect the missing lag in the AR(1) model.

#### Model specification

The AR(2) data-generating process is

\\y_t \sim \mathcal{N}(\mu_t, \sigma^2), \quad \mu_t = \alpha + \phi_1
y\_{t-1} + \phi_2 y\_{t-2}, \quad t = 3, \ldots, T\\

``` r

set.seed(42)
T_total <- 300
alpha_true <- 0.1
phi1_true <- 0.6
phi2_true <- -0.2
sigma_true <- 0.5

# Simulate AR(2) process
y_ts <- numeric(T_total)
y_ts[1:2] <- rnorm(2, alpha_true / (1 - phi1_true - phi2_true), sigma_true)
for (t in 3:T_total) {
    y_ts[t] <- alpha_true + phi1_true * y_ts[t - 1] +
               phi2_true * y_ts[t - 2] + rnorm(1, 0, sigma_true)
}

# Hold out last 50 for forecast evaluation
T_train <- 250
y_train <- y_ts[1:T_train]
y_test <- y_ts[(T_train + 1):T_total]
```

The AR(2) model:

``` r

Model_ar2 <- function(parm, Data) {
    alpha <- parm[1]
    phi1 <- parm[2]
    phi2 <- parm[3]
    sigma <- interval(parm[4], 1e-100, Inf)
    parm[4] <- sigma

    # Conditional mean
    mu <- alpha + phi1 * Data$y[2:(Data$T - 1)] +
                  phi2 * Data$y[1:(Data$T - 2)]
    y_obs <- Data$y[3:Data$T]

    # Log-likelihood
    LL <- sum(dnorm(y_obs, mu, sigma, log = TRUE))

    # Log-prior
    LP <- LL +
        dnormv(alpha, 0, 100, log = TRUE) +
        dnormv(phi1, 0, 1, log = TRUE) +
        dnormv(phi2, 0, 1, log = TRUE) +
        dhalfcauchy(sigma, 25, log = TRUE)

    # Fitted values (full length T for plotting)
    yhat <- rep(NA, Data$T)
    yhat[1:2] <- Data$y[1:2]
    for (t in 3:Data$T) {
        yhat[t] <- rnorm(1, alpha + phi1 * Data$y[t - 1] +
                         phi2 * Data$y[t - 2], sigma)
    }

    Monitor <- LP
    return(list(LP = LP, Dev = -2 * LL, Monitor = Monitor,
                yhat = yhat, parm = parm))
}

Data_ar2 <- list(
    y = y_train, T = T_train,
    mon.names = "LP",
    parm.names = c("alpha", "phi1", "phi2", "sigma")
)

IV_ar2 <- c(0, 0.5, 0, 1)
```

The AR(1) model (misspecified, omitting the second lag):

``` r

Model_ar1 <- function(parm, Data) {
    alpha <- parm[1]
    phi1 <- parm[2]
    sigma <- interval(parm[3], 1e-100, Inf)
    parm[3] <- sigma

    mu <- alpha + phi1 * Data$y[1:(Data$T - 1)]
    y_obs <- Data$y[2:Data$T]

    LL <- sum(dnorm(y_obs, mu, sigma, log = TRUE))

    LP <- LL +
        dnormv(alpha, 0, 100, log = TRUE) +
        dnormv(phi1, 0, 1, log = TRUE) +
        dhalfcauchy(sigma, 25, log = TRUE)

    yhat <- rep(NA, Data$T)
    yhat[1] <- Data$y[1]
    for (t in 2:Data$T) {
        yhat[t] <- rnorm(1, alpha + phi1 * Data$y[t - 1], sigma)
    }

    Monitor <- LP
    return(list(LP = LP, Dev = -2 * LL, Monitor = Monitor,
                yhat = yhat, parm = parm))
}

Data_ar1 <- list(
    y = y_train, T = T_train,
    mon.names = "LP",
    parm.names = c("alpha", "phi1", "sigma")
)

IV_ar1 <- c(0, 0.5, 1)
```

#### Prior predictive checks

The ribbon-style prior predictive check shows whether the prior allows
plausible time series dynamics. With \\\phi \sim \mathcal{N}(0, 1)\\,
the prior encompasses both stationary and non-stationary processes,
which is appropriate for an exploratory analysis.

``` r

rprior_ar2 <- function(n) {
    cbind(
        rnorm(n, 0, 10),        # alpha
        rnorm(n, 0, 1),         # phi1
        rnorm(n, 0, 1),         # phi2
        abs(rcauchy(n, 0, 25))  # sigma
    )
}
```

``` r

prior_predictive_check(Model_ar2, Data_ar2,
                       rprior = rprior_ar2, n = 200, type = "ribbon")
```

#### Fitting and convergence

``` r

rx_ar2 <- Prescribe(Model_ar2, Data_ar2, IV_ar2)
print(rx_ar2)
```

``` r

fit_ar2 <- lucifer(Model_ar2, Data_ar2, IV_ar2,
                   Iterations = 10000, Status = 2000,
                   Thinning = 5, Algorithm = "NUTS")

fit_ar1 <- lucifer(Model_ar1, Data_ar1, IV_ar1,
                   Iterations = 10000, Status = 2000,
                   Thinning = 5, Algorithm = "NUTS")
```

``` r

cx_ar2 <- Consort(fit_ar2)
print(cx_ar2)
```

#### Prior vs posterior with ground truth

``` r

truth_ar <- c(alpha = alpha_true, phi1 = phi1_true,
              phi2 = phi2_true, sigma = sigma_true)

prior_vs_posterior(
    fit_ar2,
    prior = rprior_ar2,
    ground_truth = truth_ar
)
```

#### Posterior predictive checks: temporal diagnostics

``` r

ppc_ar2 <- predict(fit_ar2, Model = Model_ar2, Data = Data_ar2)
ppc_ar1 <- predict(fit_ar1, Model = Model_ar1, Data = Data_ar1)
```

The time-series plot shows the observed series with the posterior
predictive mean and 95% credible ribbon. For a well-specified model, the
observed series should track within the ribbon.

``` r

plot(ppc_ar2, Style = "Time-Series")
```

The ribbon (fan chart) plot uses nested credible bands (50%, 80%, 95%)
to communicate the predictive uncertainty structure. The
`forecast_start` parameter marks the boundary between in-sample fit and
out-of-sample forecast.

``` r

ppc_ribbon(ppc_ar2)
```

To visualize out-of-sample forecast performance, we refit on the
training set and plot with `forecast_start` marking the beginning of the
holdout period. Note that this requires predictions that extend to the
full series length.

``` r

# For illustration: generate predictions over the full series
# by fitting on training data but generating yhat for the full T
ppc_ribbon(ppc_ar2, forecast_start = 201)
```

The Durbin-Watson test is the most informative temporal diagnostic. For
the AR(2) model fitted to AR(2) data, the DW statistic should be near 2
(no residual autocorrelation). For the AR(1) model fitted to AR(2) data,
the DW statistic will deviate from 2 because the missing second lag
leaves autocorrelation in the residuals.

``` r

# AR(2) model: DW near 2 (correct specification)
plot(ppc_ar2, Style = "DW")
```

``` r

# AR(1) model: DW deviates from 2 (missing lag)
plot(ppc_ar1, Style = "DW")
```

``` r

plot(ppc_ar2, Style = "Jarque-Bera")
```

``` r

plot(ppc_ar2, Style = "Residuals")
```

``` r

plot(ppc_ar2, Style = "Residual Density")
```

``` r

plot(ppc_ar2, Style = "Fitted")
```

``` r

plot(ppc_ar2, Style = "Error Scatter")
```

#### Discrepancy measures

``` r

s_ar2_mase <- summary(ppc_ar2, Discrep = "MASE")
s_ar2_mse <- summary(ppc_ar2, Discrep = "MSE")
s_ar2_ppl <- summary(ppc_ar2, Discrep = "PPL")

s_ar1_mase <- summary(ppc_ar1, Discrep = "MASE")
```

The MASE (mean absolute scaled error) is the standard accuracy metric
for time series because it is scale-free and interpretable: a MASE below
1 indicates that the model outperforms a naive random-walk forecast.

``` r

ts_comp <- data.frame(
    Model = c("AR(1)", "AR(2)"),
    BPIC = c(summary(ppc_ar1)$BPIC[1, "BPIC"],
             summary(ppc_ar2)$BPIC[1, "BPIC"]),
    Concordance = c(summary(ppc_ar1)$Concordance,
                    summary(ppc_ar2)$Concordance),
    MASE = c(s_ar1_mase$Discrepancy.Statistic,
             s_ar2_mase$Discrepancy.Statistic)
)
knitr::kable(ts_comp, digits = 3,
             caption = "PPC comparison: AR(1) vs AR(2) for AR(2) data")
```

The AR(2) model shows lower BPIC, higher concordance, and lower MASE,
confirming the correct specification. The DW test provides the most
direct evidence: the AR(1) model’s residuals retain the autocorrelation
at lag 2 that the AR(2) model captures.

``` r

arena_ts <- Arena(x = list("AR(1)" = fit_ar1, "AR(2)" = fit_ar2))
print(arena_ts)
```

### Example 5: spatio-temporal field

Spatial and spatio-temporal data require PPC diagnostics that probe
geographic patterns: whether the model reproduces spatial gradients,
whether prediction uncertainty varies sensibly across locations, and
whether temporal evolution is captured at each site. This example
generates a synthetic spatio-temporal field with known spatial
covariates and temporal random effects, then demonstrates all spatial
and multivariate temporal PPC styles, culminating in the LOO-PIT
calibration diagnostic.

#### Model specification

The data-generating process is a hierarchical spatio-temporal
regression:

\\y\_{s,t} \sim \mathcal{N}(\mu\_{s,t}, \sigma\_\epsilon^2)\\
\\\mu\_{s,t} = \beta_0 + \beta_1 \text{Elev}\_s + \beta_2
\text{Region}\_s + \beta_3 \text{Lat}\_s + \gamma_t\\ \\\gamma_t \sim
\mathcal{N}(0, \sigma\_\gamma^2)\\

where \\s = 1, \ldots, S\\ indexes spatial sites and \\t = 1, \ldots,
T\\ indexes time points.

``` r

set.seed(42)
side <- 6           # 6x6 grid; defining the side explicitly avoids
S <- side * side    # the seq(length.out = sqrt(S)) rounding trap
Ti <- 12            # time points

# Ground truth
beta0_true <- 20.0
beta1_true <- -0.5   # elevation effect
beta2_true <- 2.0    # region effect
beta3_true <- -0.3   # latitude gradient
sigma_gamma_true <- 1.5
sigma_eps_true <- 0.8

# Spatial coordinates on a regular grid
lon <- rep(seq(0, 10, length.out = side), times = side)
lat <- rep(seq(0, 10, length.out = side), each = side)

# Synthetic spatial covariates
# Elevation: increases with distance from origin
elevation <- sqrt(lon^2 + lat^2) + rnorm(S, 0, 0.5)
# Region: binary (left vs right half)
region <- as.numeric(lon > 5)
# Latitude is used directly

# Temporal random effects
gamma_true <- rnorm(Ti, 0, sigma_gamma_true)

# Generate spatio-temporal data
# y is vectorized as: y[1:S] for t=1, y[(S+1):(2*S)] for t=2, etc.
mu_st <- matrix(NA, S, Ti)
for (t in 1:Ti) {
    mu_st[, t] <- beta0_true + beta1_true * elevation +
                  beta2_true * region + beta3_true * lat + gamma_true[t]
}

y_st <- rnorm(S * Ti, as.vector(mu_st), sigma_eps_true)
Y_mat <- matrix(y_st, nrow = S, ncol = Ti)
```

``` r

Model_st <- function(parm, Data) {
    beta0 <- parm[1]
    beta1 <- parm[2]
    beta2 <- parm[3]
    beta3 <- parm[4]
    sigma_gamma <- interval(parm[5], 1e-100, Inf)
    parm[5] <- sigma_gamma
    sigma_eps <- interval(parm[6], 1e-100, Inf)
    parm[6] <- sigma_eps
    gamma <- parm[7:(6 + Data$T)]

    # Prior
    beta.prior <- sum(dnormv(c(beta0, beta1, beta2, beta3), 0, 1000,
                             log = TRUE))
    sigma_gamma.prior <- dhalfcauchy(sigma_gamma, 25, log = TRUE)
    sigma_eps.prior <- dhalfcauchy(sigma_eps, 25, log = TRUE)
    gamma.prior <- sum(dnorm(gamma, 0, sigma_gamma, log = TRUE))

    # Spatio-temporal mean: vectorized as y (S*T length)
    mu <- numeric(Data$S * Data$T)
    for (t in 1:Data$T) {
        idx <- ((t - 1) * Data$S + 1):(t * Data$S)
        mu[idx] <- beta0 + beta1 * Data$elevation +
                   beta2 * Data$region + beta3 * Data$latitude + gamma[t]
    }

    # Log-likelihood
    LL <- sum(dnorm(Data$y, mu, sigma_eps, log = TRUE))

    LP <- LL + beta.prior + sigma_gamma.prior + sigma_eps.prior + gamma.prior

    # Fitted values
    yhat <- rnorm(Data$S * Data$T, mu, sigma_eps)

    # Monitor: LP and per-observation log-likelihoods for LOO-PIT
    ll_obs <- dnorm(Data$y, mu, sigma_eps, log = TRUE)
    Monitor <- c(LP, ll_obs)

    return(list(LP = LP, Dev = -2 * LL, Monitor = Monitor,
                yhat = yhat, parm = parm))
}

Data_st <- list(
    y = y_st,
    Y = Y_mat,
    S = S,
    T = Ti,
    longitude = lon,
    latitude = lat,
    elevation = elevation,
    region = region,
    X = cbind(1, elevation, region, lat),
    mon.names = c("LP", paste0("ll[", 1:(S * Ti), "]")),
    parm.names = c("beta0", "beta1", "beta2", "beta3",
                   "sigma_gamma", "sigma_eps",
                   paste0("gamma[", 1:Ti, "]"))
)

IV_st <- c(0, 0, 0, 0, 1, 1, rep(0, Ti))
```

#### Prior predictive checks

The interval-style prior predictive check verifies that the prior
produces physically plausible predictions. For a temperature-like field
with values around 15-25, the prior predictive intervals should
encompass this range without extending to absurd values.

``` r

rprior_st <- function(n) {
    np <- 6 + Ti
    mat <- matrix(NA, n, np)
    mat[, 1:4] <- matrix(rnorm(n * 4, 0, 10), n, 4)   # beta
    mat[, 5] <- abs(rcauchy(n, 0, 25))                  # sigma_gamma
    mat[, 6] <- abs(rcauchy(n, 0, 25))                  # sigma_eps
    for (t in 1:Ti) {
        mat[, 6 + t] <- rnorm(n, 0, mat[, 5])           # gamma[t]
    }
    mat
}
```

``` r

prior_predictive_check(Model_st, Data_st,
                       rprior = rprior_st, n = 200, type = "intervals")
```

#### Fitting and convergence

``` r

rx_st <- Prescribe(Model_st, Data_st, IV_st)
print(rx_st)
```

``` r

fit_st <- lucifer(Model_st, Data_st, IV_st,
                  Iterations = 15000, Status = 3000,
                  Thinning = 5, Algorithm = "NUTS")
```

``` r

cx_st <- Consort(fit_st)
print(cx_st)
```

#### Prior vs posterior with ground truth

``` r

truth_st <- c(beta0 = beta0_true, beta1 = beta1_true,
              beta2 = beta2_true, beta3 = beta3_true,
              sigma_gamma = sigma_gamma_true,
              sigma_eps = sigma_eps_true)

prior_vs_posterior(
    fit_st,
    prior = rprior_st,
    Parms = c("beta0", "beta1", "beta2", "beta3",
              "sigma_gamma", "sigma_eps"),
    ground_truth = truth_st
)
```

#### Posterior predictive checks: spatial diagnostics

``` r

ppc_st <- predict(fit_st, Model = Model_st, Data = Data_st)
```

The spatial plot maps the posterior predictive mean (or observed values)
at each site, colored by the predicted value. Geographic patterns in the
predictions should match the known spatial covariates.

``` r

plot(ppc_st, Style = "Spatial", Data = Data_st)
```

The spatial uncertainty plot maps the width of the 95% posterior
predictive interval at each site. Sites with wider intervals have higher
prediction uncertainty, which may correlate with geographic features
(e.g., sites at extreme elevations or on the boundary between regions).

``` r

plot(ppc_st, Style = "Spatial Uncertainty", Data = Data_st)
```

The space-time decomposition plots show the temporal evolution at each
site (Space-Time by Space) and the spatial field at each time point
(Space-Time by Time). These are the most informative PPC styles for
spatio-temporal data because they reveal whether misfit is concentrated
at specific sites, specific time points, or both.

``` r

plot(ppc_st, Style = "Space-Time by Space", Data = Data_st)
```

``` r

plot(ppc_st, Style = "Space-Time by Time", Data = Data_st)
```

#### Multivariate temporal diagnostics

The multivariate time-series styles decompose the temporal comparison by
column (time point) and by row (site), providing detailed diagnostics of
how the model performs across the spatial and temporal dimensions.

``` r

plot(ppc_st, Style = "Time-Series, Multivariate, C", Data = Data_st)
```

``` r

plot(ppc_st, Style = "Time-Series, Multivariate, R", Data = Data_st)
```

``` r

plot(ppc_st, Style = "DW, Multivariate, C", Data = Data_st)
```

``` r

plot(ppc_st, Style = "Jarque-Bera, Multivariate, C", Data = Data_st)
```

Mardia’s test checks multivariate normality of the residuals across the
spatial dimension. If the residuals exhibit spatial clustering (positive
spatial autocorrelation not captured by the covariates), Mardia’s
skewness statistic will be inflated.

``` r

plot(ppc_st, Style = "Mardia", Data = Data_st)
```

``` r

plot(ppc_st, Style = "Fitted, Multivariate, C", Data = Data_st)
```

``` r

plot(ppc_st, Style = "Residuals, Multivariate, C", Data = Data_st)
```

``` r

plot(ppc_st, Style = "ECDF")
```

``` r

ppc_dens_overlay(ppc_st)
```

#### LOO-PIT calibration

The LOO-PIT diagnostic provides a single summary of model calibration.
For each observation, the PIT value is the probability that a replicated
observation is less than or equal to the observed value, computed using
PSIS-LOO importance weights. If the model is well-calibrated, the PIT
values are uniformly distributed and the ECDF follows the diagonal.

To compute LOO-PIT, we need the pointwise log-likelihoods that were
monitored during fitting. These are extracted from the Monitor matrix
and reshaped into the \\N \times S\\ format required by
[`LOO()`](https://robustecologies.github.io/lucifer/reference/LOO.md).

``` r

# Extract pointwise log-likelihoods from Monitor
# Monitor columns: LP, ll[1], ll[2], ..., ll[S*T]
posterior_monitor <- fit_st$Monitor
ll_cols <- 2:(S * Ti + 1)  # skip LP column
ll_matrix <- posterior_monitor[, ll_cols]

# Compute LOO
loo_st <- LOO(ll_matrix)

# LOO-PIT plot
ppc_loo_pit(ppc_st, loo = loo_st)
```

A well-calibrated model produces a LOO-PIT ECDF that falls within the
95% simultaneous confidence band around the diagonal. Departures
indicate specific miscalibration patterns: a U-shaped PIT histogram
indicates underdispersion (the model is overconfident), an inverted-U
indicates overdispersion (the model is too conservative), and asymmetry
indicates systematic location bias.

#### Discrepancy measures

``` r

s_st1 <- summary(ppc_st, Discrep = "L.criterion")
s_st2 <- summary(ppc_st, Discrep = "Quadratic Loss")
s_st3 <- summary(ppc_st, Discrep = "Quadratic Utility")
s_st4 <- summary(ppc_st, Discrep = "max(yhat[i,]) > max(y)")
```

``` r

st_diag <- data.frame(
    Metric = c("BPIC", "Concordance", "L-criterion"),
    Value = c(s_st1$BPIC[1, "BPIC"], s_st1$Concordance, s_st1$L.criterion)
)
knitr::kable(st_diag, digits = 3,
             caption = "Quantitative PPC summaries for the spatio-temporal model")
```

The L-criterion decomposes into a variance component (posterior
predictive uncertainty) and a bias component (systematic misfit). For a
well-specified model, the variance component dominates; a large bias
component indicates that the model is systematically over- or
under-predicting at certain locations or time points. The Quadratic Loss
and Quadratic Utility provide decision-theoretic perspectives on the
same trade-off.

### Example 6: spike-and-slab variable selection (SSVS)

Spike-and-slab priors [\[14\]](#ref14) are two-component mixture priors
that perform automatic variable selection within the Bayesian framework.
Each coefficient receives a binary inclusion indicator \\\gamma_j\\ that
switches between a narrow spike component (concentrated near zero,
effectively excluding the variable) and a broad slab component (allowing
the coefficient to take non-negligible values). The posterior
distribution of \\\gamma_j\\ directly yields the posterior inclusion
probability (PIP), which quantifies the evidence for including variable
\\j\\ in the model. This example demonstrates the dedicated SSVS
visualization tools in lucifer, which go beyond standard PPC to probe
the spike-slab structure, visualize the prior-to-posterior transition
for both the mixture components and the inclusion probabilities, and
diagnose the variable selection decision through Bayesian FDR curves and
network plots.

#### Model specification

The SSVS linear regression model assigns each coefficient (except the
intercept) a spike-and-slab prior controlled by a Bernoulli indicator:

\\y_i \sim \mathcal{N}(\mu_i, \sigma^2), \quad \mu_i = \mathbf{X}\_i
\boldsymbol{\beta}, \quad i = 1, \ldots, N\\ \\(\beta_j \mid \gamma_j)
\sim (1 - \gamma_j)\\\mathcal{N}(0, \nu_0) + \gamma_j\\\mathcal{N}(0,
\nu_1), \quad j = 2, \ldots, J\\ \\\gamma_j \sim
\text{Bernoulli}(\pi_0), \quad \beta_1 \sim \mathcal{N}(0, 1000)\\
\\\sigma \sim \text{Half-Cauchy}(25)\\

where \\\nu_0 = 0.01\\ is the spike variance (strong shrinkage toward
zero) and \\\nu_1 = 100\\ is the slab variance (diffuse prior allowing
large coefficients). The prior inclusion probability \\\pi_0\\ is set to
match the expected sparsity.

``` r

set.seed(42)
N <- 200
J <- 10  # total predictors (including intercept)

# Ground truth: 4 active predictors out of 9 (intercept always included)
beta_true <- c(3.0, 1.5, -2.0, 0, 0, 0.8, 0, 0, 0, -1.2)
gamma_true <- c(1, 1, 0, 0, 0, 1, 0, 0, 0, 1)  # indicators for beta[2:10]
sigma_true <- 1.0

X <- cbind(1, matrix(rnorm(N * (J - 1)), N, J - 1))
mu_true <- X %*% beta_true
y_ssvs <- rnorm(N, mu_true, sigma_true)
```

``` r

Model_ssvs <- function(parm, Data) {
    # Extract parameters
    beta <- parm[Data$pos.beta]
    gamma <- parm[Data$pos.gamma]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma

    # Spike-and-slab prior on beta[2:J]
    beta_sigma <- rep(sqrt(Data$slab_var), Data$J - 1)
    beta_sigma[gamma < 0.5] <- sqrt(Data$spike_var)

    beta.prior <- dnormv(beta[1], 0, 1000, log = TRUE) +
        sum(dnorm(beta[2:Data$J], 0, beta_sigma, log = TRUE))

    # Bernoulli prior on indicators
    gamma.prior <- sum(dbern(gamma, Data$prior_prob, log = TRUE))

    # Observation model
    sigma.prior <- dhalfcauchy(sigma, 25, log = TRUE)
    mu <- as.vector(Data$X %*% beta)
    LL <- sum(dnorm(Data$y, mu, sigma, log = TRUE))
    LP <- LL + beta.prior + gamma.prior + sigma.prior

    yhat <- rnorm(Data$N, mu, sigma)
    Monitor <- LP

    return(list(LP = LP, Dev = -2 * LL, Monitor = Monitor,
                yhat = yhat, parm = parm))
}

Data_ssvs <- list(
    y = y_ssvs, X = X, N = N, J = J,
    spike_var = 0.01,
    slab_var = 100,
    prior_prob = 4 / (J - 1),  # expected 4 active out of 9
    pos.beta = 1:J,
    pos.gamma = (J + 1):(J + J - 1),
    pos.sigma = 2 * J,
    ## Declare gamma as discrete for Gibbs sampling
    dparm = (J + 1):(J + J - 1),
    dsupport = NULL,  # default: {0, 1}
    mon.names = "LP",
    parm.names = c(paste0("beta[", 1:J, "]"),
                   paste0("gamma[", 1:(J - 1), "]"),
                   "sigma")
)

IV_ssvs <- c(rep(0, J), rep(0, J - 1), 1)
```

#### Prior predictive checks

The prior predictive check for an SSVS model is especially informative:
it reveals whether the spike-slab prior specification, combined with the
Bernoulli inclusion probability, produces predictions that are plausible
and appropriately uncertain.

``` r

rprior_ssvs <- function(n) {
    mat <- matrix(NA, n, 2 * J - 1 + 1)
    # beta[1] (intercept, always included)
    mat[, 1] <- rnorm(n, 0, sqrt(1000))
    # gamma[j] and beta[j] for j=2,...,J
    for (j in 2:J) {
        gj <- rbinom(n, 1, 4 / (J - 1))
        mat[, J + j - 1] <- gj  # gamma[j-1]
        spike_sd <- sqrt(0.01)
        slab_sd <- sqrt(100)
        mat[, j] <- ifelse(gj == 1,
                           rnorm(n, 0, slab_sd),
                           rnorm(n, 0, spike_sd))
    }
    # sigma
    mat[, 2 * J] <- abs(rcauchy(n, 0, 25))
    mat
}
```

``` r

prior_predictive_check(Model_ssvs, Data_ssvs,
                       rprior = rprior_ssvs, n = 300, type = "density")
```

``` r

prior_predictive_check(Model_ssvs, Data_ssvs,
                       rprior = rprior_ssvs, n = 300,
                       type = "stat", stat_fun = sd)
```

#### Fitting with Gibbs sampling

The SSVS model requires a sampler that can handle the discrete inclusion
indicators. The Gibbs sampler with auto-FC (automatic full conditionals)
is the natural choice: it uses exhaustive enumeration for the binary
indicators and slice sampling for the continuous parameters.

``` r

fit_ssvs <- lucifer(Model_ssvs, Data_ssvs, IV_ssvs,
                    Iterations = 20000, Status = 5000,
                    Thinning = 10, Algorithm = "Gibbs")
```

``` r

Consort(fit_ssvs)
```

#### SSVS-specific diagnostics with ssvs_summary()

The
[`ssvs_summary()`](https://robustecologies.github.io/lucifer/reference/ssvs_summary.md)
function is designed specifically for spike-and-slab models. It extracts
posterior inclusion probabilities, identifies the median probability
model, and provides 7 specialized plot types. The function can
auto-detect binary indicators or accept explicit parameter patterns.

``` r

# Build ground truth vector for the indicators
gt_indicators <- setNames(gamma_true[2:J],
                          paste0("gamma[", 1:(J - 1), "]"))

# Also provide ground truth for coefficients (for spike_slab overlay)
gt_coefficients <- setNames(beta_true[2:J],
                            paste0("beta[", 2:J, "]"))
gt_all <- c(gt_indicators, gt_coefficients)

ss <- ssvs_summary(
    fit_ssvs,
    indicators = "gamma",
    coefficients = paste0("beta[", 2:J, "]"),
    spike_var = 0.01,
    slab_var = 100,
    prior_prob = 4 / (J - 1),
    ground_truth = gt_all,
    labels = paste0("X", 2:J)
)

print(ss)
summary(ss)
```

#### Posterior inclusion probability bar chart

The default plot type shows the posterior inclusion probability for each
variable as a horizontal bar, with the decision threshold (0.5 for the
median probability model) as a dashed line and the prior probability as
a dotted line. Variables above the threshold are colored as selected;
ground truth status is marked with symbols.

``` r

plot(ss)
```

Variables whose PIP substantially exceeds the prior probability have
been reinforced by the data; variables whose PIP has dropped below the
prior have been penalized. The four truly active variables (\\X_2\\,
\\X_3\\, \\X_6\\, \\X\_{10}\\) should show PIPs near 1, while the five
null variables should show PIPs near 0. The transition from prior to
posterior is the variable selection in action.

#### Spike-slab prior vs posterior density

The most distinctive visualization for spike-and-slab models overlays
the spike component \\\mathcal{N}(0, 0.01)\\, the slab component
\\\mathcal{N}(0, 100)\\, and their mixture (the marginal prior) against
the marginal posterior density of each coefficient. For active
variables, the posterior should concentrate away from zero, escaping the
spike; for null variables, the posterior should collapse onto the spike.

``` r

plot(ss, type = "spike_slab")
```

This visualization answers the question that
[`prior_vs_posterior()`](https://robustecologies.github.io/lucifer/reference/prior_vs_posterior.md)
cannot address for mixture priors: how does the posterior relate to each
mixture component separately? For an active variable like \\X_2\\ (true
\\\beta = 1.5\\), the posterior mass sits entirely in the slab region,
far from the spike. For a null variable like \\X_4\\ (true \\\beta =
0\\), the posterior is indistinguishable from the spike, concentrated
tightly around zero. The ground truth values appear as vertical dashed
lines.

#### Conditional posterior density

The conditional plot splits the posterior density of each coefficient by
inclusion status: the distribution conditional on \\\gamma = 1\\
(included) and the distribution conditional on \\\gamma = 0\\
(excluded). This decomposition reveals the mechanism of the
spike-and-slab prior: when the indicator is on, the coefficient is free
to take large values; when it is off, the coefficient is constrained to
near-zero values.

``` r

plot(ss, type = "conditional")
```

For active variables, the included distribution (blue) carries most of
the posterior mass and is centered near the true value, while the
excluded distribution (grey) is negligible. For null variables, the
excluded distribution dominates and is concentrated at zero. Variables
near the decision boundary (\\\text{PIP} \approx 0.5\\) show substantial
mass in both distributions, indicating genuine uncertainty about
inclusion.

#### Running PIP trajectory

The trajectory plot shows the cumulative posterior inclusion probability
as a function of MCMC iteration. This monitors convergence of the
variable selection decision and detects indicators that are slow to mix.

``` r

plot(ss, type = "trajectory")
```

Well-mixing indicators converge quickly to their final PIP, while
poorly-mixing ones show persistent oscillation. Indicators near the
decision boundary naturally show more variability, but prolonged drift
indicates inadequate chain length.

#### Bayesian FDR and FOR curves

When ground truth is available (as in simulation studies), the Bayesian
false discovery rate (FDR = proportion of selected variables that are
truly null) and false omission rate (FOR = proportion of excluded
variables that are truly active) can be computed as functions of the
inclusion threshold. This helps calibrate the threshold beyond the
default 0.5.

``` r

plot(ss, type = "fdr")
```

The FDR curve decreases as the threshold increases (stricter selection
admits fewer false positives), while the FOR curve increases (stricter
selection misses more true positives). The optimal threshold balances
these two errors. For well-separated signals, the FDR drops to zero well
before threshold = 0.5, and the median probability model achieves zero
FDR and zero FOR.

#### Standard PPC for the SSVS model

The SSVS model also produces `yhat` for standard posterior predictive
checks. These diagnostics evaluate the predictive adequacy of the entire
model, including the variable selection.

``` r

ppc_ssvs <- predict(fit_ssvs, Model = Model_ssvs, Data = Data_ssvs)
```

``` r

ppc_dens_overlay(ppc_ssvs)
```

``` r

plot(ppc_ssvs, Style = "Fitted", Data = Data_ssvs)
```

``` r

plot(ppc_ssvs, Style = "Covariates", Data = Data_ssvs)
```

The covariate plots are particularly informative for SSVS: active
predictors show a clear relationship between the covariate and the
posterior predictive mean, while null predictors show flat (zero-slope)
relationships, confirming that the spike-and-slab prior has correctly
zeroed out the irrelevant coefficients.

``` r

s_ssvs <- summary(ppc_ssvs, Discrep = "Chi-Square")
```

#### VAR-SSVS: heatmap and network visualizations

For multivariate time series with SSVS priors on the autoregressive
coefficients, the inclusion indicators form a matrix (or array, for
multiple lags). The heatmap and network plot types are designed for this
structure, visualizing which directed connections between variables are
supported by the data. We simulate a sparse VAR(1) process with known
\\\Gamma\\ structure, fit it with Gibbs sampling using auto-FC for the
binary indicators, and reshape the inclusion probabilities back into the
original \\J \times J\\ array via the `dims` argument of
[`ssvs_summary()`](https://robustecologies.github.io/lucifer/reference/ssvs_summary.md).

``` r

set.seed(42)
J_var <- 4
TT <- 150

# Ground truth: sparse VAR(1) connectivity (6 active edges out of 16)
Gamma_true <- matrix(c(
    1, 1, 0, 0,
    0, 1, 1, 0,
    0, 0, 1, 0,
    0, 0, 0, 1
), J_var, J_var, byrow = TRUE)

Phi_true <- Gamma_true * matrix(c(
     0.50,  0.40,  0.00, 0.00,
     0.00,  0.50,  0.30, 0.00,
     0.00,  0.00,  0.50, 0.00,
     0.00,  0.00,  0.00, 0.50), J_var, J_var, byrow = TRUE)

# Verify stationarity
stopifnot(all(Mod(eigen(Phi_true)$values) < 1))

# Simulate VAR(1)
Y_var <- matrix(0, TT, J_var)
Y_var[1, ] <- rnorm(J_var)
for (t in 2:TT) {
    Y_var[t, ] <- Phi_true %*% Y_var[t - 1, ] + rnorm(J_var, 0, 0.3)
}
```

The model places a spike-and-slab prior on each autoregressive
coefficient \\\Phi\_{ij}\\ controlled by a binary indicator
\\\Gamma\_{ij}\\, exactly as in the univariate case but now over the
full \\J \times J\\ matrix. The product \\\Gamma\_{ij}\Phi\_{ij}\\ in
the mean function ensures that excluded coefficients contribute nothing
to the likelihood regardless of their drawn value.

``` r

Model_var_ssvs <- function(parm, Data) {
    Gamma_vec <- parm[Data$pos.Gamma]
    Phi_vec   <- parm[Data$pos.Phi]
    sigma     <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma

    Gamma <- matrix(Gamma_vec, Data$J, Data$J)
    Phi   <- matrix(Phi_vec,   Data$J, Data$J)

    # Spike-and-slab prior on Phi
    Phi_sigma <- ifelse(Gamma == 1,
                        sqrt(Data$slab_var),
                        sqrt(Data$spike_var))

    Gamma.prior <- sum(dbern(Gamma_vec, Data$prior_prob, log = TRUE))
    Phi.prior   <- sum(dnorm(Phi_vec, 0,
                             as.vector(Phi_sigma), log = TRUE))
    sigma.prior <- sum(dhalfcauchy(sigma, 25, log = TRUE))

    # Likelihood: Y[t,] = (Gamma * Phi) %*% Y[t-1,] + N(0, sigma^2 I)
    A  <- Gamma * Phi
    mu <- Data$Y[1:(Data$TT - 1), ] %*% t(A)
    Sigma_mat <- matrix(sigma, Data$TT - 1, Data$J, byrow = TRUE)
    LL <- sum(dnorm(Data$Y[2:Data$TT, ], mu, Sigma_mat, log = TRUE))

    LP <- LL + Gamma.prior + Phi.prior + sigma.prior
    yhat <- rnorm(length(mu), mu, Sigma_mat)
    list(LP = LP, Dev = -2 * LL, Monitor = LP,
         yhat = yhat, parm = parm)
}

parm.names_var <- c(
    paste0("Gamma[", rep(1:J_var, J_var), ",",
           rep(1:J_var, each = J_var), "]"),
    paste0("Phi[",   rep(1:J_var, J_var), ",",
           rep(1:J_var, each = J_var), "]"),
    paste0("sigma[", 1:J_var, "]")
)

K <- J_var * J_var
Data_var <- list(
    Y = Y_var, J = J_var, TT = TT,
    spike_var = 0.01, slab_var = 1.0,
    prior_prob = 0.5,
    pos.Gamma = 1:K,
    pos.Phi   = (K + 1):(2 * K),
    pos.sigma = (2 * K + 1):(2 * K + J_var),
    parm.names = parm.names_var,
    mon.names = "LP",
    ## Declare Gamma indicators as discrete Bernoulli for auto-FC Gibbs
    dparm = 1:K,
    dsupport = NULL
)

IV_var <- c(rep(1, K),       # start with dense Gamma (all on)
            rep(0.05, K),    # Phi near zero
            rep(0.5, J_var)) # sigma
```

``` r

fit_var_ssvs <- lucifer(Model_var_ssvs, Data_var, IV_var,
                        Iterations = 2000, Status = 1000,
                        Thinning = 5, Algorithm = "Gibbs")
```

After fitting the VAR-SSVS model,
[`ssvs_summary()`](https://robustecologies.github.io/lucifer/reference/ssvs_summary.md)
extracts the posterior inclusion probabilities and, when supplied with
`dims = c(J_var, J_var)`, reshapes them into the original
adjacency-matrix layout. The ground truth vector is constructed in
column-major order to match the parameter naming convention used
throughout lucifer.

``` r

gt_var <- setNames(as.vector(Gamma_true),
                   paste0("Gamma[",
                          rep(1:J_var, J_var), ",",
                          rep(1:J_var, each = J_var), "]"))

ss_var <- ssvs_summary(
    fit_var_ssvs,
    indicators   = "Gamma",
    coefficients = paste0("Phi[",
                          rep(1:J_var, J_var), ",",
                          rep(1:J_var, each = J_var), "]"),
    spike_var = 0.01, slab_var = 1.0,
    prior_prob = 0.5,
    ground_truth = gt_var,
    dims = c(J_var, J_var)
)
print(ss_var)
```

``` r

plot(ss_var, type = "heatmap")
```

The heatmap uses a diverging color scale centered at 0.5 (the decision
boundary). Tiles above 0.5 are colored blue (selected), tiles below 0.5
are white (excluded). Text overlays show the numeric PIP and, when
ground truth is available, the true indicator value in parentheses. For
the simulated process, the four diagonal (own-lag) connections and the
two true off-diagonal connections (\\\Gamma\_{1,2}\\ and
\\\Gamma\_{2,3}\\) should appear in blue, while the remaining tiles
should remain near zero.

``` r

plot(ss_var, type = "network")
```

The network plot represents each time series as a node on a circle and
draws directed edges between nodes whose connection has \\\text{PIP} \>
0.1\\. Edge opacity and width are proportional to the PIP, making strong
connections visually prominent and weak connections faint. Edges above
the threshold are colored blue; edges below are grey. The resulting
graph immediately communicates the Granger-causal structure identified
by the VAR-SSVS model: arrows from \\X_2 \to X_1\\ and \\X_3 \to X_2\\
should appear, while the diagonal self-loops are absorbed by the node
markers themselves.

## Comparative analysis

### When PPC reveals misspecification

Across the six examples, PPC revealed misspecification through distinct
signatures. In Example 1, the density overlay showed that a Gaussian
model applied to heavy-tailed data produces replicated datasets that are
too concentrated near the center, missing the extremes; the concordance
dropped below 0.95 because observations in the tails fell outside their
predictive intervals. In Example 3, the rootogram showed a
characteristic “spike at zero” pattern where the Poisson model
systematically underproduced zeros; the test statistic \\T(y) =
\text{proportion of zeros}\\ gave a posterior predictive p-value near 0,
providing a scalar summary of the same failure. In Example 4, the
Durbin-Watson statistic deviated from 2 for the AR(1) model, indicating
that the residuals retained the lag-2 autocorrelation that the missing
term could not absorb.

In Example 6, the spike-slab density overlay showed the mechanism of
Bayesian variable selection in action: the posterior of active
coefficients escaped the spike and concentrated in the slab region,
while null coefficients collapsed onto the spike. The PIP bar chart,
conditional density, and Bayesian FDR curves provided complementary
perspectives on the same selection decision.

These patterns share a common structure: the PPC diagnostic that is most
informative depends on the aspect of the data that the model is failing
to capture. Distributional misfit shows up in density overlays and ECDF
comparisons. Zero-inflation shows up in rootograms and zero-proportion
test statistics. Temporal structure shows up in DW tests and time-series
plots. Spatial structure shows up in spatial maps and Mardia’s test.
Variable selection structure shows up in spike-slab overlays and PIP
trajectories. The lesson is that no single PPC diagnostic is sufficient;
a comprehensive evaluation requires diagnostics targeted at each
structural aspect of the model.

### Integrating PPC with model comparison

PPC provides model-specific diagnostics (how well does this model
reproduce the data?), while Arena provides cross-model comparisons
(which model reproduces the data best?). The two are complementary.
BPIC, computed from the posterior predictive deviance, serves as a
bridge: it is a PPC-based measure (derived from the predictive
distribution) that also functions as a model selection criterion (lower
is better, with a penalty for effective complexity).

The Arena framework in lucifer unifies PPC-based metrics (BPIC,
concordance, L-criterion) with efficiency metrics (ESS per second),
convergence diagnostics (Rhat, ESS), and information-theoretic measures
(DIC, WAIC) into a single comparison table. This integration ensures
that model selection accounts for both statistical adequacy and
computational practicality. A model that achieves marginally better BPIC
but requires 100 times more computation may not be worth the cost, and
Arena’s Pareto frontier visualization makes this trade-off explicit.

### The role of LOO-PIT in calibration assessment

The LOO-PIT diagnostic stands apart from the other PPC tools because it
provides a genuinely out-of-sample assessment: each observation’s PIT
value is computed using the leave-one-out posterior predictive
distribution, which conditions on all observations except the one being
evaluated. This makes LOO-PIT immune to the criticism that posterior
predictive checks use the data twice (once to fit, once to check),
because the LOO posterior predictive distribution has not seen the
observation it is being asked to predict.

In practice, the PSIS approximation [\[9\]](#ref9) avoids the
computational cost of refitting the model \\N\\ times by
importance-weighting the existing posterior samples. The quality of this
approximation is monitored by the Pareto \\k\\ diagnostic: values below
0.7 indicate reliable LOO-PIT values; values above 0.7 suggest that some
observations are highly influential and the PSIS approximation may be
unreliable for those points.

The LOO-PIT ECDF is the single most informative summary of overall model
calibration. If it falls within the 95% confidence band, the model is
adequately calibrated in the aggregate; if it deviates systematically,
the pattern of deviation diagnoses the nature of the miscalibration.

## Summary of all PPC plot styles

The following table lists all 39 posterior predictive plot styles, their
data requirements, the example in this vignette where each is
demonstrated, and guidance on when each style is most informative.

``` r

summary_table <- data.frame(
    Style = c(
        "Density", "Density Overlay", "Histogram Overlay", "ECDF Overlay",
        "ECDF", "Fitted", "Scatter Average", "Bars",
        "Intervals", "Ribbon", "Predictive Quantiles",
        "Residuals", "Residual Density", "Error Scatter", "Error Histogram",
        "DW", "Jarque-Bera", "Mardia", "Stat", "Stat 2D",
        "Rootogram",
        "Covariates", "Covariates, Categorical DV", "Violin Grouped",
        "Time-Series",
        "Spatial", "Spatial Uncertainty",
        "Space-Time by Space", "Space-Time by Time",
        "LOO-PIT",
        "Fitted, Multivariate, C", "Fitted, Multivariate, R",
        "Residuals, Multivariate, C", "Residuals, Multivariate, R",
        "Residual Density, Multivariate, C",
        "Residual Density, Multivariate, R",
        "Time-Series, Multivariate, C",
        "Time-Series, Multivariate, R",
        "DW, Multivariate, C",
        "Jarque-Bera, Multivariate, C"
    ),
    Example = c(
        "1", "1,3,5", "1", "1",
        "1,2,5", "1,4", "1,3", "2,3",
        "1,2", "4", "1,3",
        "1,4", "1,4", "1,4", "1",
        "4", "2,4", "5", "1,2,3", "1,2",
        "3",
        "1", "3", "2",
        "4",
        "5", "5",
        "5", "5",
        "5",
        "2,5", "2",
        "2,5", "2",
        "2",
        "2",
        "5",
        "5",
        "5",
        "5"
    ),
    When_to_use = c(
        "Per-observation predictive density",
        "Overall distributional comparison",
        "Overall distributional comparison (discrete-friendly)",
        "Cumulative distributional comparison",
        "Aggregate calibration check",
        "Prediction accuracy per observation",
        "Overall prediction accuracy",
        "Category-level comparison (count/discrete data)",
        "Credible intervals per observation (shrinkage patterns)",
        "Time series with nested credible bands",
        "Outlier detection via predictive quantiles",
        "Residual patterns (heteroscedasticity, trends)",
        "Residual distribution shape",
        "Residuals vs fitted (heteroscedasticity)",
        "Error distribution shape",
        "Residual autocorrelation (time series)",
        "Residual normality",
        "Multivariate residual normality",
        "Scalar test statistic (customizable)",
        "Bivariate test statistic (customizable)",
        "Count data calibration (hanging rootogram)",
        "Covariate effects on predictions",
        "Covariate effects for categorical outcomes",
        "Group-level predictive distributions",
        "Temporal fit with uncertainty",
        "Spatial prediction map",
        "Spatial uncertainty map",
        "Per-site temporal dynamics",
        "Per-time spatial snapshot",
        "Leave-one-out calibration (requires LOO object)",
        "Per-variable fit (multivariate response)",
        "Per-observation fit (multivariate response)",
        "Per-variable residuals (multivariate)",
        "Per-observation residuals (multivariate)",
        "Per-variable residual density",
        "Per-observation residual density",
        "Per-variable temporal dynamics",
        "Per-observation temporal dynamics",
        "Per-variable autocorrelation test",
        "Per-variable normality test"
    )
)

knitr::kable(summary_table,
             col.names = c("Style", "Example", "When to use"),
             caption = "Complete reference: all 39 PPC plot styles in lucifer")
```

In addition, four prior predictive plot types are available via
[`prior_predictive_check()`](https://robustecologies.github.io/lucifer/reference/prior_predictive_check.md):
density overlay, intervals, ribbon, and test statistic. Six convenience
wrappers provide a bayesplot-compatible API:
[`ppc_dens_overlay()`](https://robustecologies.github.io/lucifer/reference/ppc_dens_overlay.md),
[`ppc_intervals()`](https://robustecologies.github.io/lucifer/reference/ppc_intervals.md),
[`ppc_ribbon()`](https://robustecologies.github.io/lucifer/reference/ppc_ribbon.md),
[`ppc_rootogram()`](https://robustecologies.github.io/lucifer/reference/ppc_rootogram.md),
[`ppc_stat()`](https://robustecologies.github.io/lucifer/reference/ppc_stat.md),
and
[`ppc_loo_pit()`](https://robustecologies.github.io/lucifer/reference/ppc_loo_pit.md).
The
[`prior_vs_posterior()`](https://robustecologies.github.io/lucifer/reference/prior_vs_posterior.md)
function compares prior and posterior densities with optional ground
truth, and
[`prior_sensitivity()`](https://robustecologies.github.io/lucifer/reference/prior_sensitivity.md)
overlays posteriors under different prior specifications.

## References

**\[1\]** Rubin, D.B. (1984). Bayesianly justifiable and relevant
frequency calculations for the applied statistician. *Annals of
Statistics*, 12(4), 1151-1172.
[doi:10.1214/aos/1176346785](https://doi.org/10.1214/aos/1176346785)

**\[2\]** Gelman, A., Meng, X.L., and Stern, H. (1996). Posterior
predictive assessment of model fitness via realized discrepancies.
*Statistica Sinica*, 6, 733-807.

**\[3\]** Meng, X.L. (1994). Posterior predictive p-values. *Annals of
Statistics*, 22, 1142-1160.

**\[4\]** Gabry, J., Simpson, D., Vehtari, A., Betancourt, M., and
Gelman, A. (2019). Visualization in Bayesian workflow. *Journal of the
Royal Statistical Society: Series A*, 182(2), 389-402.
[doi:10.1111/rssa.12378](https://doi.org/10.1111/rssa.12378)

**\[5\]** Gelfand, A. (1996). Model determination using sampling-based
methods. In W.R. Gilks, S. Richardson, and D.J. Spiegelhalter (Eds.),
*Markov Chain Monte Carlo in Practice*, pp. 145-162. Chapman and Hall.

**\[6\]** Gelfand, A. and Ghosh, S. (1998). Model choice: a minimum
posterior predictive loss approach. *Biometrika*, 85(1), 1-11.
[doi:10.1093/biomet/85.1.1](https://doi.org/10.1093/biomet/85.1.1)

**\[7\]** Ando, T. (2007). Bayesian predictive information criterion for
the evaluation of hierarchical Bayesian and empirical Bayes models.
*Biometrika*, 94(2), 443-458.
[doi:10.1093/biomet/asm017](https://doi.org/10.1093/biomet/asm017)

**\[8\]** Laud, P.W. and Ibrahim, J.G. (1995). Predictive model
selection. *Journal of the Royal Statistical Society: Series B*, 57,
247-262.

**\[9\]** Vehtari, A., Gelman, A., and Gabry, J. (2017). Practical
Bayesian model evaluation using leave-one-out cross-validation and WAIC.
*Statistics and Computing*, 27(5), 1413-1432.
[doi:10.1007/s11222-016-9696-4](https://doi.org/10.1007/s11222-016-9696-4)

**\[10\]** Kleiber, C. and Zeileis, A. (2016). Visualizing count data
regressions using rootograms. *The American Statistician*, 70(3),
296-303.
[doi:10.1080/00031305.2016.1173590](https://doi.org/10.1080/00031305.2016.1173590)

**\[11\]** Durbin, J. and Watson, G.S. (1950). Testing for serial
correlation in least squares regression. *Biometrika*, 37, 409-428.
[doi:10.1093/biomet/37.3-4.409](https://doi.org/10.1093/biomet/37.3-4.409)

**\[12\]** Jarque, C.M. and Bera, A.K. (1980). Efficient tests for
normality, homoscedasticity and serial independence of regression
residuals. *Economics Letters*, 6(3), 255-259.
[doi:10.1016/0165-1765(80)90024-5](https://doi.org/10.1016/0165-1765(80)90024-5)

**\[13\]** Mardia, K.V. (1970). Measures of multivariate skewness and
kurtosis with applications. *Biometrika*, 57(3), 519-530.
[doi:10.1093/biomet/57.3.519](https://doi.org/10.1093/biomet/57.3.519)

**\[14\]** George, E.I. and McCulloch, R.E. (1993). Variable selection
via Gibbs sampling. *Journal of the American Statistical Association*,
88(423), 881-889.
[doi:10.1080/01621459.1993.10476353](https://doi.org/10.1080/01621459.1993.10476353)

**\[15\]** Barbieri, M.M. and Berger, J.O. (2004). Optimal predictive
model selection. *Annals of Statistics*, 32(3), 870-897.
[doi:10.1214/009053604000000238](https://doi.org/10.1214/009053604000000238)
