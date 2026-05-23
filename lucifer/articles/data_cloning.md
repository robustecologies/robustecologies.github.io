# Data cloning: maximum likelihood estimation via Bayesian MCMC

## Introduction

Data cloning is a computational method introduced by Lele, Dennis, and
Lutscher [\[1\]](#ref1) that redirects the Bayesian MCMC machinery to
produce frequentist maximum likelihood estimates and their standard
errors. The method is particularly valuable for complex hierarchical
models where the likelihood integral has no closed-form solution, where
numerical optimization is unreliable due to multimodality or high
dimensionality, and where standard asymptotic approximations (such as
those underlying `optim` or the Laplace method) break down. Data cloning
requires no gradient computation, no analytical evaluation of the
normalizing constant, and no reformulation of the model; any model that
can be expressed in the Bayesian framework can be used directly.

The
[`DataCloning()`](https://robustecologies.github.io/lucifer/reference/DataCloning.md)
function in lucifer implements data cloning with automatic convergence
diagnostics, identifiability testing, and profile likelihood
computation. It leverages the full suite of 55+ MCMC algorithms,
multi-chain parallel infrastructure, and modern diagnostics already
available in the package.

## Theoretical foundations

### The data cloning principle

Consider observed data \\y = (y_1, \ldots, y_n)\\ and a parametric model
with likelihood \\L(\theta \mid y) = \prod\_{i=1}^n f(y_i \mid \theta)\\
and prior \\\pi(\theta)\\. The standard Bayesian posterior is

\\\pi(\theta \mid y) \propto L(\theta \mid y) \\ \pi(\theta).\\

The data cloning idea is to construct an artificial dataset \\y^{(K)}\\
by replicating \\y\\ exactly \\K\\ times. The likelihood of the cloned
data is \\L(\theta \mid y^{(K)}) = \[L(\theta \mid y)\]^K\\, and the
corresponding posterior is

\\\pi_K(\theta \mid y) \propto \[L(\theta \mid y)\]^K \\ \pi(\theta).\\

As \\K\\ increases, the factor \\\[L(\theta \mid y)\]^K\\ concentrates
exponentially around the maximum likelihood estimate
\\\hat{\theta}\_{\mathrm{MLE}}\\, overwhelming the prior contribution
\\\pi(\theta)\\. In the limit, the posterior becomes a point mass at the
MLE regardless of the prior, provided the MLE exists and is unique.

### Asymptotic results

The theoretical justification rests on three asymptotic properties
established by Walker [\[2\]](#ref2) in the context of increasing sample
sizes and applied to data cloning by Lele et al. [\[1\]](#ref1):

**Convergence of the posterior mean.** The posterior mean under cloning
converges to the MLE:

\\\bar{\theta}\_K = E\[\theta \mid y^{(K)}\] \xrightarrow{K \to \infty}
\hat{\theta}\_{\mathrm{MLE}}.\\

**Convergence of the scaled variance.** The posterior covariance matrix,
scaled by \\K\\, converges to the inverse of the observed Fisher
information matrix:

\\K \cdot \mathrm{Var}\_K(\theta) \xrightarrow{K \to \infty}
I(\hat{\theta})^{-1}\\

where \\I(\theta) = -\nabla^2 \log L(\theta \mid y)\big\|\_{\theta =
\hat{\theta}}\\ is the observed Fisher information. This provides an
automatic estimate of the standard errors of the MLE without requiring
computation of second derivatives.

**Asymptotic normality.** The cloned posterior converges in distribution
to a multivariate normal:

\\\theta \mid y^{(K)} \xrightarrow{d}
N\\\left(\hat{\theta}\_{\mathrm{MLE}},\\ K^{-1}
I(\hat{\theta})^{-1}\right).\\

This result mirrors the classical Bernstein–von Mises theorem for
increasing sample sizes, except that data cloning achieves the
concentration by repeating the same observations rather than collecting
new ones.

### Connection to prior feedback and deterministic annealing

Robert [\[3\]](#ref3) introduced the prior feedback approach, in which
the posterior from one Bayesian analysis is used as the prior for the
next, iteratively concentrating the posterior on the MLE. Data cloning
achieves the same effect in a single step by raising the likelihood to
the \\K\\th power. The connection to deterministic annealing is also
direct: define the tempered posterior \\\pi\_\beta(\theta) \propto
L(\theta)^\beta \\ \pi(\theta)\\ with inverse temperature \\\beta\\.
Data cloning corresponds to setting \\\beta = K\\, where \\K\\ is an
integer. Unlike simulated annealing, data cloning uses a discrete
sequence of temperatures \\K = 1, 2, 4, \ldots\\ with full MCMC
convergence at each level, and the target is the MLE rather than a
global minimum of an arbitrary objective.

### The LP/Dev decomposition

The lucifer Model interface guarantees that every model evaluation
returns the log-posterior \\\mathrm{LP} = \ell(\theta) + \log
\pi(\theta)\\ and the deviance \\\mathrm{Dev} = -2\\\ell(\theta)\\. This
algebraic structure allows exact decomposition of the log-likelihood and
the log-prior without any modification to the user’s Model function:

\\\ell(\theta) = -\frac{\mathrm{Dev}}{2}, \qquad \log \pi(\theta) =
\mathrm{LP} + \frac{\mathrm{Dev}}{2}.\\

The cloned model then computes

\\\mathrm{LP}\_K = K \cdot \ell(\theta) + \log \pi(\theta), \qquad
\mathrm{Dev}\_K = K \cdot \mathrm{Dev}.\\

This is implemented as a thin wrapper around the Model function. No
actual data duplication is performed, and all MCMC algorithms,
multi-chain support, and diagnostics operate transparently on the
wrapped model.

## Basic usage: normal model

Consider the simplest possible case: a normal model with unknown mean
\\\mu\\ and known variance \\\sigma^2 = 1\\, observed on \\n = 100\\
independent observations. The MLE is the sample mean \\\hat{\mu} =
\bar{y}\\ and the asymptotic standard error is \\\mathrm{SE} =
1/\sqrt{n} = 0.1\\.

``` r

library(lucifer)

set.seed(42)
n <- 100
true_mu <- 3.0
y <- rnorm(n, mean = true_mu, sd = 1)

Data <- list(
    mon.names = "LP",
    parm.names = "mu",
    y = y, n = n
)

Model <- function(parm, Data) {
    mu <- parm[1]
    LL <- sum(dnorm(Data$y, mu, 1, log = TRUE))
    LP <- LL + dnorm(mu, 0, 100, log = TRUE)
    list(LP = LP, Dev = -2 * LL, Monitor = LP,
         yhat = rep(mu, Data$n), parm = parm)
}
```

The prior is intentionally diffuse (\\\mu \sim N(0, 100^2)\\) to
demonstrate that data cloning overcomes prior influence. Running
[`DataCloning()`](https://robustecologies.github.io/lucifer/reference/DataCloning.md)
with five clone sizes:

``` r

dc <- DataCloning(
    Model = Model,
    Data = Data,
    Initial.Values = 0,
    n_clones = c(1, 2, 4, 8, 16),
    Iterations = 10000,
    Thinning = 10,
    Algorithm = "NUTS",
    Chains = 3L,
    CPUs = 3L,
    update_prior = TRUE,
    verbose = TRUE
)
```

### Interpreting the output

The `print` method provides a concise summary showing the MLE, standard
error, and Wald confidence interval:

``` r

print(dc)
```

We can verify against the analytical solution:

``` r

verify_tbl <- data.frame(
    Quantity = c("MLE", "SE", "CI lower", "CI upper"),
    Analytical = c(mean(y), 1 / sqrt(n),
                   mean(y) - 1.96 / sqrt(n),
                   mean(y) + 1.96 / sqrt(n)),
    `Data cloning` = c(dc$MLE, dc$SE, dc$Wald.CI[1], dc$Wald.CI[2]),
    check.names = FALSE
)
knitr::kable(verify_tbl, digits = 4, align = "lrr",
             caption = "Comparison of analytical and data cloning estimates for the normal model.")
```

The `summary` method provides extended diagnostics including per-K
tables, eigenvalue ratios, scaled variance convergence, and the
MSE/R-squared statistic:

``` r

summary(dc)
```

## Convergence diagnostics

Data cloning provides three complementary diagnostics for assessing
whether the sequence of clone sizes is sufficient. All three are
computed automatically and displayed in the `summary` output and the
`plot` methods.

### Convergence of posterior means

The primary diagnostic is the convergence plot, which shows the
posterior mean \\\pm\\ standard error as a function of \\K\\ on a
log-scale axis. As \\K\\ increases, the posterior should concentrate at
the MLE (horizontal dashed line) with monotonically shrinking
uncertainty:

``` r

plot(dc, type = "convergence")
```

If the posterior mean has not stabilized at the largest \\K\\, more
clone sizes should be added to the sequence. The `update_prior = TRUE`
option (warm starting) generally accelerates convergence by using the
posterior from \\K\\ as the starting point for \\K+1\\.

### Eigenvalue diagnostic

Lele, Nadeem, and Schmuland [\[4\]](#ref4) proposed monitoring the
largest eigenvalue of the posterior covariance matrix, standardized by
its value at \\K = 1\\:

\\\lambda_S(K) =
\frac{\lambda\_{\max}(\Sigma_K)}{\lambda\_{\max}(\Sigma_1)}.\\

Under identifiability, \\\lambda_S(K)\\ should decrease at rate \\1/K\\.
The eigenvalue plot shows \\\lambda_S(K)\\ against \\K\\ on log-log
axes, with a dashed \\1/K\\ reference line. Departure from the reference
line indicates potential non-identifiability: if \\\lambda_S(K)\\ levels
off rather than continuing to decrease, at least one direction in
parameter space is not being constrained by the data.

``` r

plot(dc, type = "eigenvalue")
```

### Scaled variance

The quantity \\K \cdot \mathrm{Var}\_K(\theta_j)\\ should stabilize
(become approximately constant) as \\K\\ increases, converging to the
\\j\\th diagonal element of \\I(\hat{\theta})^{-1}\\. The scaled
variance plot shows this quantity against \\K\\ for each parameter. A
horizontal trajectory indicates that the asymptotic regime has been
reached; a decreasing or increasing trajectory suggests that more clones
are needed.

``` r

plot(dc, type = "scaled_variance")
```

### Density overlay

The density overlay shows the posterior density for each \\K\\ value,
superimposed on the same axes. As \\K\\ increases, the density should
narrow and center on the MLE (vertical dashed line). This provides an
intuitive visual confirmation that the prior is being overwhelmed by the
likelihood amplification:

``` r

plot(dc, type = "density")
```

### Diagnostic dashboard

The diagnostics dashboard combines all four plots in a single 2x2 panel:

``` r

plot(dc, type = "diagnostics")
```

## Example: logistic regression

Data cloning is not limited to simple models. Consider a logistic
regression with two parameters (intercept and slope):

``` r

set.seed(42)
n <- 200
x <- rnorm(n)
true_beta <- c(0.5, 1.2)
prob <- plogis(true_beta[1] + true_beta[2] * x)
y <- rbinom(n, 1, prob)

Data_lr <- list(
    mon.names = "LP",
    parm.names = c("beta0", "beta1"),
    X = cbind(1, x),
    y = y, n = n
)

Model_lr <- function(parm, Data) {
    beta <- parm
    eta <- Data$X %*% beta
    p <- plogis(eta)
    p <- pmin(pmax(p, 1e-10), 1 - 1e-10)
    LL <- sum(dbinom(Data$y, 1, p, log = TRUE))
    LP <- LL + sum(dnorm(beta, 0, 10, log = TRUE))
    list(LP = LP, Dev = -2 * LL, Monitor = LP,
         yhat = as.numeric(p), parm = parm)
}

dc_lr <- DataCloning(
    Model = Model_lr,
    Data = Data_lr,
    Initial.Values = c(0, 0),
    n_clones = c(1, 2, 4, 8, 16),
    Iterations = 10000,
    Algorithm = "NUTS",
    Chains = 3L,
    CPUs = 3L,
    update_prior = TRUE,
    verbose = TRUE
)
```

The data cloning MLE should closely match the
[`glm()`](https://rdrr.io/r/stats/glm.html) output:

``` r

glm_fit <- glm(y ~ x, family = binomial)
glm_coefs <- summary(glm_fit)$coefficients
compare_tbl <- data.frame(
    Parameter = dc_lr$parm.names,
    `glm MLE` = glm_coefs[, 1],
    `DC MLE` = dc_lr$MLE,
    `glm SE` = glm_coefs[, 2],
    `DC SE` = dc_lr$SE,
    check.names = FALSE
)
knitr::kable(compare_tbl, digits = 4, align = "lrrrr",
             caption = "Comparison of glm() and data cloning estimates for the logistic model.")
```

The convergence and pairs plots for the logistic regression:

``` r

plot(dc_lr, type = "convergence")
plot(dc_lr, type = "pairs")
```

## Identifiability and estimability

A critical advantage of data cloning over standard numerical MLE methods
is its ability to detect structural non-identifiability. When a
parameter is not identifiable, the posterior does not concentrate as
\\K\\ increases; instead, it remains diffuse and prior-dependent.
Campbell and Lele [\[5\]](#ref5) formalized this observation through a
two-stage ANOVA test.

### The ANOVA estimability test

The test requires running data cloning under \\P \geq 2\\ distinct
priors and examining how the posterior mean of each parameter depends on
the clone size \\K\\ and the prior index \\p\\. For each parameter
\\\theta_j\\, the posterior means \\\bar{\theta}\_{j,k,p}\\ across clone
sizes \\k = 1, \ldots, K\_{\max}\\ and priors \\p = 1, \ldots, P\\ are
analyzed in a two-stage ANOVA:

**Stage 1 (clone effect).** Test \\H_0\\: the posterior mean does not
depend on \\K\\. Under this null, cloning has no effect, which implies
the parameter cannot be estimated from the data. If the F-test is
significant (\\p \< 0.05\\), the posterior is responding to the
increasing likelihood amplification, a necessary condition for
estimability.

**Stage 2 (prior effect).** Conditional on the clone effect, test
\\H_0\\: the posterior mean does not depend on the prior. If the F-test
is significant, the prior still influences the estimate even at large
\\K\\, meaning the data alone cannot determine the parameter value. This
indicates non-identifiability.

The classification of each parameter follows from these two tests:

- **Estimable**: significant clone effect, non-significant prior effect.
  The posterior concentrates on a unique MLE regardless of the prior.
- **Weakly estimable**: significant clone effect, significant prior
  effect. The data provides some information but the prior still
  matters.
- **Not estimable**: non-significant clone effect. The posterior does
  not concentrate regardless of \\K\\. Structural non-identifiability.

### Non-identifiable model example

Consider a model where only the sum \\\mu_1 + \mu_2\\ is identifiable,
but neither \\\mu_1\\ nor \\\mu_2\\ individually:

\\y_i \sim N(\mu_1 + \mu_2, 1), \quad i = 1, \ldots, n.\\

The likelihood depends only on \\\mu_1 + \mu_2\\, so there are
infinitely many MLEs along the ridge \\\mu_1 + \mu_2 = \bar{y}\\.

``` r

set.seed(999)
n <- 100
true_sum <- 5.0
y <- rnorm(n, mean = true_sum, sd = 1)

Data_ni <- list(
    mon.names = "LP",
    parm.names = c("mu1", "mu2"),
    y = y, n = n
)

Model_ni <- function(parm, Data) {
    mu1 <- parm[1]
    mu2 <- parm[2]
    mu_sum <- mu1 + mu2
    LL <- sum(dnorm(Data$y, mu_sum, 1, log = TRUE))
    LP <- LL + dnorm(mu1, 0, 10, log = TRUE) + dnorm(mu2, 0, 10, log = TRUE)
    list(LP = LP, Dev = -2 * LL, Monitor = LP,
         yhat = rep(mu_sum, Data$n), parm = parm)
}

dc_ni <- DataCloning(
    Model = Model_ni,
    Data = Data_ni,
    Initial.Values = c(2.5, 2.5),
    n_clones = c(1, 2, 4, 8, 16),
    Iterations = 10000,
    Status = 1000,
    Algorithm = "NUTS",
    Chains = 1L,
    CPUs = 1L,
    n_priors = 3,
    verbose = FALSE
)

summary(dc_ni)
```

The ANOVA test should flag both \\\mu_1\\ and \\\mu_2\\ as weakly
estimable or not estimable, while the sum \\\mu_1 + \mu_2\\ (were it
monitored) would be fully estimable.

### Visualizing estimability

The estimability plot shows posterior means from different priors at
each \\K\\. For identifiable parameters, the means converge to the same
value regardless of prior; for non-identifiable parameters, the means
remain separated:

``` r

plot(dc_ni, type = "estimability")
```

The prior sensitivity plot provides a forest plot showing the MLE
\\\pm\\ SE for each prior at the highest \\K\\. Large variation across
priors confirms non-identifiability:

``` r

plot(dc_ni, type = "prior_sensitivity")
```

The eigenvalue plot will show that \\\lambda_S(K)\\ levels off rather
than continuing to decrease at rate \\1/K\\, reflecting the fact that
one direction in parameter space (the \\\mu_1 - \mu_2\\ direction) is
not being constrained:

``` r

plot(dc_ni, type = "eigenvalue")
```

The “diagnostics” plot provides a comprehensive view of the structural
unidentifiability of the model:

``` r

plot(dc_ni, type = "diagnostics")
```

## Profile likelihood

Profile likelihood provides confidence intervals that are more accurate
than Wald intervals for nonlinear models and small samples
[\[6\]](#ref6). For a single parameter \\\theta_j\\, the profile
log-likelihood is defined as

\\\ell_P(\theta_j) = \max\_{\theta\_{-j}} \ell(\theta_j, \theta\_{-j})\\

where the maximization is over all parameters except \\\theta_j\\. The
profile likelihood confidence interval at level \\1 - \alpha\\ is the
set of values \\\theta_j\\ satisfying

\\2\left\[\ell_P(\hat{\theta}\_j) - \ell_P(\theta_j)\right\] \leq
\chi^2\_{1, 1-\alpha}\\

where \\\chi^2\_{1, 1-\alpha}\\ is the \\(1-\alpha)\\ quantile of the
chi-squared distribution with 1 degree of freedom.

The `plot(dc, type = "profile")` method computes profile likelihood
curves on demand. This requires the original Model and Data to be passed
as additional arguments, since the profile computation performs
numerical optimization at each grid point:

``` r

plot(dc_lr, type = "profile",
     Model = Model_lr, Data = Data_lr,
     n_grid = 60)
```

The horizontal dashed line marks the \\-\chi^2_1(0.95)/2\\ threshold.
Where the profile curve crosses this line defines the 95% profile
likelihood confidence interval, which is invariant to reparameterization
and does not assume normality of the estimator.

The profile curve for parameters \\\mu_1\\ and \\\mu_2\\ of the
structurally unidentifiable model clearly shows the problem of the
infinitely many MLEs along the ridge \\\mu_1 + \mu_2 = \bar{y}\\:

``` r

plot(dc_ni, type = "profile",
     Model = Model_ni, Data = Data_ni,
     n_grid = 60)
```

## Model comparison via data cloning

Data cloning provides a natural route to likelihood-based model
comparison. The maximized log-likelihood at the MLE can be extracted
from the final-K run, since \\\mathrm{Dev}\_K = K \cdot \mathrm{Dev}\\
implies

\\\hat{\ell} = -\frac{\mathrm{Dev}\_{K\_{\max}}}{2 K\_{\max}}\\

evaluated at the last iteration (which should be near the MLE). This
enables computation of the AIC, BIC, or the data cloning likelihood
ratio (DCLR) test described by Ponciano et al. [\[6\]](#ref6):

``` r

K_max <- max(dc$n_clones)
Dev_final <- dc$final_fit$Deviance[length(dc$final_fit$Deviance)]
LL_at_MLE <- -Dev_final / (2 * K_max)

AIC_dc <- -2 * LL_at_MLE + 2 * dc$LIV
BIC_dc <- -2 * LL_at_MLE + log(Data$n) * dc$LIV

ic_tbl <- data.frame(
    Criterion = c("Log-likelihood at MLE", "AIC", "BIC"),
    Value = c(LL_at_MLE, AIC_dc, BIC_dc)
)
knitr::kable(ic_tbl, digits = 4, align = "lr",
             caption = "Information criteria extracted from the data cloning fit.")
```

For comparing two nested models \\M_1 \subset M_2\\, the DCLR test
statistic is

\\\Lambda = 2\left\[\hat{\ell}\_2 - \hat{\ell}\_1\right\]\\

which has a \\\chi^2\\ distribution with \\p_2 - p_1\\ degrees of
freedom under the null hypothesis that \\M_1\\ is adequate.

## Practical guidance

### Choosing clone sizes

The sequence of clone sizes should be geometric (e.g., \\K = 1, 2, 4, 8,
16, 32\\) to span a wide range efficiently. Starting with \\K = 1\\
provides the unmodified Bayesian posterior as a reference. The largest
\\K\\ should be chosen so that the scaled variance \\K \cdot
\mathrm{Var}\_K\\ has stabilized; monitoring the convergence diagnostics
across the sequence reveals whether more clones are needed. In practice,
\\K\_{\max}\\ between 16 and 64 suffices for most well-identified
models. Non-identifiable models will show failure to converge at any
\\K\\.

### Choosing the MCMC algorithm

Adaptive algorithms (NUTS, AMWG, RAM) handle the narrowing posterior at
high \\K\\ automatically through their adaptation mechanisms. The
No-U-Turn Sampler (NUTS) is the default and generally the best choice,
since it adapts its step size and trajectory length to the local
curvature. For high-dimensional models, Hamiltonian Monte Carlo variants
are strongly preferred over random-walk Metropolis. Setting
`update_prior = TRUE` (warm starting) is recommended for all algorithms,
as it substantially reduces the adaptation burden at each new \\K\\
level.

### Warm starting

When `update_prior = TRUE`, the posterior mean and covariance from clone
\\K\\ serve as the initial values and proposal covariance for clone
\\K+1\\. This is particularly important at large \\K\\ where the
posterior becomes very narrow: a sampler initialized at the prior mean
would waste most of its iterations exploring the tails before finding
the high-density region. Warm starting ensures that sampling begins near
the mode with an appropriately scaled proposal, reducing the effective
burn-in to a fraction of what would otherwise be required.

### Automatic prior generation for estimability testing

When `n_priors > 1` and no `prior_generator` is supplied,
[`DataCloning()`](https://robustecologies.github.io/lucifer/reference/DataCloning.md)
generates alternative priors automatically. A natural question is
whether these should span different distributional families (normal,
uniform, Cauchy, etc.) to maximize the diversity of prior
specifications. In practice, distributional family is irrelevant to the
estimability test. The Campbell and Lele [\[5\]](#ref5) ANOVA checks
whether the posterior mean converges to the same value regardless of the
prior; what matters is that the alternative priors place mass at
different locations in parameter space, not that the mass has a
particular shape. A normal prior centered at \\\bar{\theta}\_1 +
3\\\mathrm{SD}\_1\\ and a uniform prior with similar location and spread
would produce essentially the same posterior shift at low \\K\\ and the
same convergence at high \\K\\. The identifiability signal comes from
location and informativeness, not functional form.

The current implementation uses location-shifted normal priors derived
from the K=1 posterior of the original model: after running the first
clone level, the posterior mean \\\bar{\theta}\_1\\ and standard
deviation \\\mathrm{SD}\_1\\ are computed, and alternative priors are
centered at \\\bar{\theta}\_1 \pm 3\\\mathrm{SD}\_1\\ with standard
deviation \\\mathrm{SD}\_1\\. This guarantees a detectable prior effect
at low \\K\\ (the informative priors pull the posterior away from the
MLE) that vanishes for estimable parameters as \\K\\ grows and the
cloned likelihood overwhelms the prior. The ANOVA classification is
further cross-validated against the eigenvalue and scaled-variance
diagnostics, so even when the ANOVA lacks statistical power (e.g.,
because a very diffuse original prior already places the posterior near
the MLE at \\K = 1\\), false “not estimable” classifications are caught.

The one situation where custom prior families matter is when the
parameter space has constraints (e.g., variance parameters that must be
positive, correlation parameters in \\(-1, 1)\\, or simplex-valued
parameters). In such cases, the user should supply a `prior_generator`
function that respects the parameter constraints and places mass in
meaningfully different regions of the feasible space.

### Parallelization across clone sizes

The current implementation parallelizes chains within each clone level
(via `Chains` and `CPUs`) but runs the clone sizes sequentially. One
might ask whether the \\K\\ values could also be parallelized, running
\\K = 1, 4, 16\\ simultaneously on separate groups of cores for a
potential speedup proportional to the number of clone levels.

This is architecturally problematic for two reasons. First, warm
starting (`update_prior = TRUE`, the recommended default) creates a
sequential dependency: clone \\K+1\\ uses the posterior mean and
covariance from clone \\K\\ as its initial values and proposal, so
\\K+1\\ cannot start until \\K\\ finishes. Parallelizing across \\K\\
only works when `update_prior = FALSE`, which sacrifices significant
sampling efficiency at high \\K\\ where the posterior is very narrow and
adaptation from scratch is expensive. Second, resource accounting
becomes fragile. If `Chains = 3` and three \\K\\ values run in parallel,
nine CPU-bound processes compete for cores. On a typical 8-core machine
this causes thrashing, not speedup; the optimal allocation (how many
\\K\\ batches \\\times\\ how many chains per batch) depends on the
hardware, the algorithm’s per-iteration cost, and the operating system’s
scheduling behaviour. A flat worker pool with dynamic task assignment
would solve this but requires restructuring the chain infrastructure.

For most workflows (5-6 clone levels with warm starting), sequential
execution across \\K\\ adds negligible overhead because each MCMC run is
short (the posterior is narrow and converges quickly thanks to warm
starting). If speed is a concern, the most effective levers are: (a)
increasing `CPUs` for more parallel chains within each \\K\\, (b) using
NUTS rather than random-walk Metropolis (fewer iterations needed for the
same effective sample size), and (c) enabling `update_prior = TRUE` to
reduce burn-in at each successive clone level.

### Diagnostics to check

After running
[`DataCloning()`](https://robustecologies.github.io/lucifer/reference/DataCloning.md),
the following should be verified:

1.  The convergence plot (`type = "convergence"`) should show the
    posterior mean stabilizing at a horizontal asymptote (the MLE) with
    shrinking error bars.
2.  The eigenvalue plot (`type = "eigenvalue"`) should show
    \\\lambda_S(K)\\ decreasing at rate \\1/K\\ (parallel to the dashed
    reference line on log-log axes). Leveling off indicates
    non-identifiability.
3.  The scaled variance plot (`type = "scaled_variance"`) should show
    \\K \cdot \mathrm{Var}\_K\\ becoming flat. A non-constant trajectory
    indicates that the asymptotic regime has not been reached.
4.  The R-hat values reported per clone level should all be below 1.05,
    confirming adequate MCMC convergence within each run.

If any of these diagnostics fail, the remedies are: increase
`Iterations`, increase the largest clone size, enable
`update_prior = TRUE`, or switch to a more efficient algorithm.

## Limitations

The DIC computed from cloned runs is not interpretable as a model
comparison criterion, since \\\mathrm{Dev}\_K = K \cdot \mathrm{Dev}\\
inflates the deviance by a factor of \\K\\. Use the DCLR approach
described above for model comparison instead.

At very large \\K\\, the cloned log-likelihood \\K \cdot \ell(\theta)\\
can underflow to \\-\infty\\ at parameter values far from the MLE.
[`DataCloning()`](https://robustecologies.github.io/lucifer/reference/DataCloning.md)
validates that the cloned LP is finite at the initial values and skips
any \\K\\ that produces non-finite values, with a warning. Warm starting
mitigates this issue by ensuring that the initial values remain close to
the mode.

Data cloning assumes that the MLE exists and is unique. For models with
flat ridges in the likelihood surface (e.g., mixture models with label
switching, or models with redundant parameters), the posterior will not
concentrate and the diagnostics will correctly flag non-identifiability.
In such cases, the ANOVA estimability test (`n_priors > 1`) provides a
formal assessment.

## References

**\[1\]** Lele, S.R., Dennis, B. & Lutscher, F. (2007). Data cloning:
easy maximum likelihood estimation for complex ecological models using
Bayesian Markov chain Monte Carlo methods. *Ecology Letters*, 10,
551-563.
[doi:10.1111/j.1461-0248.2007.01047.x](https://doi.org/10.1111/j.1461-0248.2007.01047.x)

**\[2\]** Walker, A.M. (1969). On the asymptotic behaviour of posterior
distributions. *J. Roy. Statist. Soc. B*, 31, 80-88.
[doi:10.1111/j.2517-6161.1969.tb00767.x](https://doi.org/10.1111/j.2517-6161.1969.tb00767.x)

**\[3\]** Robert, C.P. (1993). Prior feedback: a Bayesian approach to
maximum likelihood estimation. *Comput. Statist.*, 8, 279-294.

**\[4\]** Lele, S.R., Nadeem, K. & Schmuland, B. (2010). Estimability
and likelihood inference for generalized linear mixed models using data
cloning. *J. Amer. Statist. Assoc.*, 105(492), 1617-1625.
[doi:10.1198/jasa.2010.tm09757](https://doi.org/10.1198/jasa.2010.tm09757)

**\[5\]** Campbell, D. & Lele, S. (2014). An ANOVA test for parameter
estimability using data cloning with application to statistical
inference for dynamic systems. *Comput. Statist. Data Anal.*, 70,
257-267.
[doi:10.1016/j.csda.2013.09.013](https://doi.org/10.1016/j.csda.2013.09.013)

**\[6\]** Ponciano, J.M., Taper, M.L., Dennis, B. & Lele, S.R. (2009).
Hierarchical models in ecology: confidence intervals, hypothesis
testing, and model selection using data cloning. *Ecology*, 90, 356-362.
[doi:10.1890/07-1960.1](https://doi.org/10.1890/07-1960.1)
