# The frequentist bridge: maximum likelihood inference with Bayesian machinery in lucifer

## Introduction

Most applied statisticians learn inference through the frequentist lens:
point estimates from [`lm()`](https://rdrr.io/r/stats/lm.html), standard
errors from the Hessian, \\p\\-values from \\t\\-tests, confidence
intervals from the normal approximation. This vignette shows that every
one of those quantities can be computed by lucifer, a package whose
primary purpose is Bayesian inference, without altering any underlying
theory. The result is a computational bridge: the same MCMC, Laplace,
and quadrature algorithms that sample posteriors can, under specific
conditions, deliver maximum likelihood estimates, Fisher information,
Wald tests, and likelihood ratio statistics that match the output of
[`lm()`](https://rdrr.io/r/stats/lm.html),
[`glm()`](https://rdrr.io/r/stats/glm.html),
[`nls()`](https://rdrr.io/r/stats/nls.html), and
[`nlme::lme()`](https://rdrr.io/pkg/nlme/man/lme.html) to several
decimal places.

Why would a frequentist use Bayesian machinery? Three reasons. First,
lucifer’s 55 MCMC algorithms, variational methods, and Laplace
optimizers handle models that admit no closed-form likelihood gradient
or Hessian, including hierarchical, censored, and nonlinear structures
where classical optimizers frequently fail. Second, data cloning
[\[7\]](#ref7) provides automatic identifiability diagnostics, Fisher
information estimation from posterior variance, and profile likelihood
curves, none of which standard optimizers offer out of the box. Third,
the bridge is pedagogically valuable: a user who verifies that
flat-prior MCMC reproduces the [`lm()`](https://rdrr.io/r/stats/lm.html)
output has a concrete foundation from which to introduce informative
priors and discover when and why the two paradigms diverge.

The vignette proceeds as follows. Section 2 develops the theoretical
foundations: the Bernstein-von Mises theorem, Fisher information, and
the Wald-LR-Score test trinity. Section 3 demonstrates three routes to
the MLE on a single linear regression problem, comparing against
[`lm()`](https://rdrr.io/r/stats/lm.html). Sections 4 through 7 extend
the bridge to generalized linear models, nonlinear regression, mixed
effects, and survival analysis. Section 8 examines three scenarios where
frequentist and Bayesian answers genuinely diverge. Section 9 validates
coverage via Monte Carlo simulation. Section 10 offers practical
guidelines for users transitioning between paradigms.

## Theoretical foundations

### The Bernstein-von Mises theorem

The central result connecting frequentist and Bayesian inference is the
Bernstein-von Mises theorem. Under standard regularity conditions (a
correctly specified, identifiable model with a smooth, positive prior at
the true parameter value), the posterior distribution converges in total
variation to a normal distribution centered at the MLE with covariance
equal to the inverse observed Fisher information, scaled by the sample
size:

\\\left\\ \pi_n(\theta \mid y\_{1:n}) -
\mathcal{N}\\\left(\hat\theta\_{\text{MLE}},\\ n^{-1}
I(\hat\theta)^{-1}\right) \right\\\_{\text{TV}}
\xrightarrow{P\_{\theta_0}} 0 \quad \text{as} \quad n \to \infty\\

where \\\pi_n\\ denotes the posterior density,
\\\hat\theta\_{\text{MLE}}\\ is the maximum likelihood estimator, and
\\I(\theta) = -E\[\nabla^2 \log f(Y \mid \theta)\]\\ is the expected
Fisher information matrix [\[8\]](#ref8).

The practical consequence is immediate: for sufficiently large \\n\\,
Bayesian credible intervals and frequentist confidence intervals have
the same asymptotic coverage. The posterior mean converges to the MLE,
the posterior standard deviation converges to the standard error, and
the posterior quantiles converge to the normal quantiles used for Wald
intervals. This equivalence holds regardless of the prior, provided the
prior is positive and continuous in a neighborhood of the true
parameter; the data overwhelm the prior as \\n\\ grows.

This theorem justifies the bridge. When a lucifer model specifies flat
priors (so that the log-posterior equals the log-likelihood), the
posterior IS the likelihood surface, and any mode-finding or sampling
algorithm that explores it is performing maximum likelihood estimation.
The posterior mean is the MLE, the posterior covariance is the inverse
observed information, and credible intervals are confidence intervals.

The theorem fails in several instructive cases. When the true parameter
lies on the boundary of the parameter space (as when a variance
component equals zero), the posterior may be concentrated on the
boundary half-space rather than converging to a normal. When the model
is misspecified, the posterior concentrates on the pseudo-true value
(the Kullback-Leibler minimizer), but its width does not match the
sandwich standard error that frequentists use for robust inference
[\[10\]](#ref10). When the parameter space is infinite-dimensional
(nonparametric models), the theorem may fail entirely [\[5\]](#ref5).
These failures are examined in Section 8.

### Fisher information and the Hessian

The frequentist covariance of the MLE is estimated by the inverse of the
observed Fisher information:

\\\hat{V} = J(\hat\theta)^{-1}, \qquad J(\hat\theta) = -\nabla^2
\ell(\hat\theta)\\

where \\\ell(\theta) = \sum\_{i=1}^n \log f(y_i \mid \theta)\\ is the
log-likelihood. lucifer implements three methods for estimating this
matrix through the
[`CovEstim()`](https://robustecologies.github.io/lucifer/reference/Matrices.md)
function, accessible via the `CovEst` argument of
[`LaplaceApproximation()`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md):

The **Hessian** method (default) computes \\J(\hat\theta)\\ via
Richardson extrapolation using the
[`Hessian()`](https://robustecologies.github.io/lucifer/reference/Matrices.md)
function. This is the most accurate approach for moderate dimensionality
and is equivalent to what [`nlm()`](https://rdrr.io/r/stats/nlm.html)
and `optim(method="BFGS", hessian=TRUE)` provide.

The **OPG** (outer product of gradients) method estimates the
information as \\\hat{I}(\hat\theta) = \sum\_{i=1}^n s_i(\hat\theta)
s_i(\hat\theta)^T\\, where \\s_i(\theta) = \nabla \log f(y_i \mid
\theta)\\ is the score contribution of observation \\i\\ [\[2\]](#ref2).
This requires \\N + NJ\\ model evaluations and is asymptotically
equivalent to the Hessian estimator under correct specification.

The **Sandwich** estimator \\\hat{V}\_S = J^{-1} \hat{I} J^{-1}\\ is
robust to model misspecification [\[10\]](#ref10). Under correct
specification, \\I = J\\ and the sandwich reduces to the Hessian
inverse. Under misspecification, the sandwich gives asymptotically valid
standard errors while the Hessian inverse does not. This distinction
becomes important in Section 8.3.

### The Wald-LR-Score trinity

Three classical tests arise naturally from the likelihood function. All
are asymptotically equivalent under the null hypothesis \\H_0\\: R\theta
= r\\, but differ in finite samples.

The **Wald test** evaluates the constraint at the unrestricted MLE: \\W
= (R\hat\theta - r)^T \left\[R \hat{V} R^T\right\]^{-1} (R\hat\theta -
r) \\\sim\\ \chi^2(q)\\

where \\q = \text{rank}(R)\\. This is the test reported by
`summary(lm(...))` and `summary(glm(...))` for individual coefficients
(with \\R\\ being a row selector and \\r = 0\\). lucifer provides this
through
[`wald_test()`](https://robustecologies.github.io/lucifer/reference/wald_test.md).

The **likelihood ratio test** compares the maximized log-likelihoods of
the full and reduced models: \\\Lambda =
2\left(\hat\ell\_{\text{full}} - \hat\ell\_{\text{reduced}}\right)
\\\sim\\ \chi^2(q)\\

This is the most reliable of the three in finite samples
[\[11\]](#ref11). lucifer provides it through
[`lr_test()`](https://robustecologies.github.io/lucifer/reference/lr_test.md),
which extracts the maximized log-likelihood from the deviance component
of any fit object (recall that every lucifer Model returns
`Dev = -2\ell`).

The **score test** evaluates the gradient of the log-likelihood at the
restricted estimate: \\S = U(\hat\theta_0)^T I(\hat\theta_0)^{-1}
U(\hat\theta_0) \\\sim\\ \chi^2(q)\\

where \\U(\theta) = \nabla \ell(\theta)\\ is the score vector. Its
advantage is that only the restricted model needs to be fitted. lucifer
provides this through
[`score_test()`](https://robustecologies.github.io/lucifer/reference/score_test.md),
which computes the gradient via
[`partial()`](https://robustecologies.github.io/lucifer/reference/Math.md)
and the information via
[`CovEstim()`](https://robustecologies.github.io/lucifer/reference/Matrices.md).
Because the gradient is estimated by finite differences, the score test
is less numerically stable than the other two and should be treated as
experimental.

Under \\H_0\\, all three statistics converge to the same \\\chi^2(q)\\
distribution. In finite samples, the Wald test tends to overreject (its
rejection region is approximately elliptical in parameter space, which
can be a poor approximation of the true confidence region for nonlinear
models), the LR test has the best overall behavior, and the score test
lies between them [\[4\]](#ref4). For routine use,
[`wald_test()`](https://robustecologies.github.io/lucifer/reference/wald_test.md)
and
[`lr_test()`](https://robustecologies.github.io/lucifer/reference/lr_test.md)
are recommended.

### Data cloning as a bridge

Data cloning provides the most principled frequentist extraction from
Bayesian machinery. Rather than eliminating priors (which forces the
user to write a non-standard Model function), data cloning retains the
prior but amplifies the likelihood by replicating the data \\K\\ times.
As \\K \to \infty\\, the posterior concentrates on the MLE regardless of
the prior, and the scaled posterior covariance \\K \cdot
\text{Var}\_K(\theta)\\ converges to the inverse Fisher information
[\[9\]](#ref9), [\[7\]](#ref7).

Data cloning is covered in detail in the [data cloning
vignette](https://robustecologies.github.io/lucifer/articles/data_cloning.md).
The key advantage over flat-prior MCMC is that it provides automatic
identifiability testing via the eigenvalue diagnostic and the ANOVA
estimability test [\[3\]](#ref3), and profile likelihood curves for
asymmetric confidence intervals. Any Model that works with
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
works with
[`DataCloning()`](https://robustecologies.github.io/lucifer/reference/DataCloning.md)
without modification.

## Three roads to the MLE

We demonstrate the three mechanisms on a simple normal linear regression
with known ground truth, comparing against
[`lm()`](https://rdrr.io/r/stats/lm.html) throughout.

### Reference: `lm()`

``` r

set.seed(42)
n <- 200
x1 <- rnorm(n)
x2 <- rnorm(n)
true_beta <- c(1.5, 0.8, -0.6)
true_sigma <- 0.5
y <- true_beta[1] + true_beta[2] * x1 + true_beta[3] * x2 + rnorm(n, 0, true_sigma)
lm_fit <- lm(y ~ x1 + x2)
lm_tab <- coef(summary(lm_fit))
knitr::kable(lm_tab, digits = 4, caption = "Reference: lm() output")
```

### Data and model for lucifer

The Model function below sets `LP = LL` (no priors), making the
posterior identical to the likelihood surface. We parameterize
\\\sigma\\ on the log scale to ensure an unconstrained parameter space,
which yields a well-conditioned Hessian.

``` r

library(lucifer)
X <- cbind(1, x1, x2)
Data <- list(
    parm.names = c("beta0", "beta1", "beta2", "log.sigma"),
    mon.names  = "LL",
    X = X, y = y, n = n, N = n
)
Model <- function(parm, Data) {
    beta  <- parm[1:3]
    sigma <- exp(parm[4])
    mu    <- Data$X %*% beta
    LL    <- sum(dnorm(Data$y, mu, sigma, log = TRUE))
    list(LP = LL, Dev = -2 * LL, Monitor = LL,
         yhat = rnorm(length(mu), mu, sigma), parm = parm)
}
Initial.Values <- c(0, 0, 0, 0)
```

### Laplace approximation

The fastest route to the MLE is
[`LaplaceApproximation()`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
which performs deterministic optimization. With `CovEst = "Hessian"`, it
estimates the covariance via Richardson-extrapolated second derivatives,
exactly the same computation as
[`nlm()`](https://rdrr.io/r/stats/nlm.html) with `hessian = TRUE`.

``` r

fit_la <- LaplaceApproximation(Model, Initial.Values, Data,
                               Iterations = 200, Method = "SPG",
                               CovEst = "Hessian", sir = FALSE)
fs_la <- freq_summary(fit_la, df.residual = n - 4)
print(fs_la)
```

### Flat-prior MCMC

The Model above has `LP = LL`, so any MCMC sampler applied to it will
sample from the likelihood surface. The posterior mean converges to the
MLE and the posterior covariance converges to the inverse observed
information. This route is slower than Laplace but provides the full
sampling distribution and MCMC diagnostics.

``` r

fit_mcmc <- lucifer(Model, Data, Initial.Values,
                    Algorithm = "NUTS", Iterations = 10000, Status = 2000,
                    Specs = list(A = 500, delta = 0.8, epsilon = NULL,
                                 Lmax = 20))
Consort(fit_mcmc)
fs_mcmc <- freq_summary(fit_mcmc, df.residual = n - 4)
print(fs_mcmc)
plot(fit_mcmc, Data=Data, type = 'all')
```

### Data cloning

Data cloning retains the Model unmodified (priors can be present or
absent) and progressively amplifies the likelihood. The MLE and Fisher
information are extracted automatically.

``` r

fit_dc <- DataCloning(Model, Data, Initial.Values,
                      n_clones = c(1, 2, 4, 8, 16),
                      Algorithm = "NUTS",
                      Iterations = 5000, Status = 1000,
                      Specs = list(A = 500, delta = 0.8,
                                   epsilon = NULL, Lmax = 20))
print(fit_dc)
summary(fit_dc)
fs_dc <- freq_summary(fit_dc, df.residual = n - 4)
print(fs_dc)
plot(fit_dc, type = "diagnostics")
```

### Comparison

The table below compares all four routes to the MLE. All estimates
should agree to 2-3 significant figures for a well-specified linear
regression with \\n = 200\\.

``` r

## Build comparison from all four fits
comp <- data.frame(
    lm.Estimate      = lm_tab[, 1],
    lm.SE            = lm_tab[, 2],
    Laplace.Estimate = fs_la$coefficients[1:3, "Estimate"],
    Laplace.SE       = fs_la$coefficients[1:3, "Std.Error"],
    MCMC.Estimate    = fs_mcmc$coefficients[1:3, "Estimate"],
    MCMC.SE          = fs_mcmc$coefficients[1:3, "Std.Error"],
    DC.Estimate      = fs_dc$coefficients[1:3, "Estimate"],
    DC.SE            = fs_dc$coefficients[1:3, "Std.Error"]
)
rownames(comp) <- c("beta0", "beta1", "beta2")
knitr::kable(comp, digits = 4,
             caption = "Comparison of lm(), Laplace, MCMC, and data cloning")
```

The estimates match because the Bernstein-von Mises theorem guarantees
asymptotic equivalence, and with \\n = 200\\ the asymptotics have kicked
in. The standard errors agree because all methods invert the same
Hessian at the same mode (or, for MCMC, estimate the posterior
covariance from samples, which converges to the inverse Hessian). The
only systematic difference arises from the parameterization:
[`lm()`](https://rdrr.io/r/stats/lm.html) estimates \\\sigma\\ directly
(via RSS\\/n\\), while the lucifer Model estimates \\\log\sigma\\ and
the SE refers to that transformed parameter.

The
[`confint_compare()`](https://robustecologies.github.io/lucifer/reference/confint_compare.md)
function provides a visual comparison. It produces a forest plot where
each method’s point estimate and confidence interval are overlaid for
every parameter, with the true values marked as diamonds.

``` r

confint_compare(
    list(Laplace = fit_la, MCMC = fit_mcmc, DataCloning = fit_dc),
    true_values = c(beta0 = true_beta[1], beta1 = true_beta[2],
                    beta2 = true_beta[3], log.sigma = log(true_sigma)),
    parm = c("beta0", "beta1", "beta2")
)
```

Individual `freq_summary` objects also have a
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) method that
displays coefficient estimates with confidence intervals. The Laplace
result below includes the true parameter values as red crosses for
reference.

``` r

plot(fs_la, true_values = c(beta0 = true_beta[1], beta1 = true_beta[2],
                            beta2 = true_beta[3],
                            log.sigma = log(true_sigma)))
```

### Residual diagnostics

The
[`freq_residuals()`](https://robustecologies.github.io/lucifer/reference/freq_residuals.md)
function produces the classical four-panel residual diagnostic display
familiar from `plot(lm(...))`. It re-evaluates the Model at the MLE
(averaging over multiple calls to smooth out stochastic `yhat`
components), computes raw and standardized residuals, leverage values
from the hat matrix, and Cook’s distances.

``` r

rd <- freq_residuals(fit_la, Model, Data)
print(rd)
plot(rd)
```

The four panels show residuals vs fitted values (to detect non-linearity
or heteroscedasticity), a normal Q-Q plot (to assess normality of
residuals), a scale-location plot (to check homoscedasticity), and
residuals vs leverage with Cook’s distance contours (to identify
influential observations). For a well-specified linear regression the
residuals should scatter randomly around zero, the Q-Q plot should lie
close to the diagonal, and no observations should exhibit extreme
leverage or Cook’s distance.

### Profile likelihood

Wald confidence intervals assume the log-likelihood is quadratic around
the MLE, which holds well for linear models but can be a poor
approximation for nonlinear or boundary parameters. The
[`profile_likelihood()`](https://robustecologies.github.io/lucifer/reference/profile_likelihood.md)
function computes profile log-likelihood curves by fixing each parameter
at a grid of values and re-optimizing over the nuisance parameters,
providing confidence intervals that respect the actual curvature of the
likelihood surface.

``` r

prof <- profile_likelihood(fit_la, Model, Data)
print(prof)
plot(prof)
```

For the linear regression example the profile and Wald intervals nearly
coincide, confirming that the quadratic approximation is adequate. The
summary method reports an asymmetry ratio for each parameter; values
near 1.0 indicate symmetric confidence regions where the Wald
approximation suffices, while substantial departures (common in
nonlinear models; see the Michaelis-Menten example below) signal that
profile intervals should be preferred.

``` r

summary(prof)
```

## Generalized linear models

### Logistic regression vs `glm(family = binomial)`

Binary outcomes with a logit link provide the first non-trivial test.
The log-likelihood is \\\ell(\beta) = \sum_i \left\[y_i \log p_i + (1 -
y_i) \log(1 - p_i)\right\]\\ where \\p_i = \text{logit}^{-1}(X_i^T
\beta)\\.

``` r

set.seed(123)
n <- 300
x1 <- rnorm(n)
x2 <- rnorm(n)
true_beta_logistic <- c(0.5, 1.2, -0.8)
p <- plogis(true_beta_logistic[1] + true_beta_logistic[2] * x1 +
            true_beta_logistic[3] * x2)
y_bin <- rbinom(n, 1, p)
```

``` r

glm_fit <- glm(y_bin ~ x1 + x2, family = binomial)
knitr::kable(coef(summary(glm_fit)), digits = 4,
             caption = "Reference: glm() logistic regression")
```

``` r

X_logistic <- cbind(1, x1, x2)
Data_logistic <- list(
    parm.names = c("beta0", "beta1", "beta2"),
    mon.names  = "LL",
    X = X_logistic, y = y_bin, n = n, N = n
)
Model_logistic <- function(parm, Data) {
    eta <- Data$X %*% parm
    p   <- plogis(eta)
    p   <- pmin(pmax(p, 1e-10), 1 - 1e-10)
    LL  <- sum(dbinom(Data$y, 1, p, log = TRUE))
    list(LP = LL, Dev = -2 * LL, Monitor = LL,
         yhat = rbinom(Data$n, 1, p), parm = parm)
}
fit_logistic <- LaplaceApproximation(Model_logistic, c(0, 0, 0),
                                     Data_logistic, Iterations = 200,
                                     Method = "SPG", sir = FALSE)
fs_logistic <- freq_summary(fit_logistic)
print(fs_logistic)
plot(fit_logistic, Data = Data_logistic)
```

``` r

comp_logistic <- data.frame(
    glm.Estimate     = coef(summary(glm_fit))[, 1],
    glm.SE           = coef(summary(glm_fit))[, 2],
    Laplace.Estimate = fs_logistic$coefficients[, "Estimate"],
    Laplace.SE       = fs_logistic$coefficients[, "Std.Error"]
)
knitr::kable(comp_logistic, digits = 4,
             caption = "Logistic regression: glm() vs Laplace")
```

The estimates and standard errors match because
[`glm()`](https://rdrr.io/r/stats/glm.html) uses iteratively reweighted
least squares (Fisher scoring), which converges to the same mode that
lucifer’s SPG optimizer finds, and both invert the same Hessian.

#### Hypothesis testing

We now demonstrate the three classical tests. The Wald test on
individual coefficients reproduces
`summary(glm(...))$coefficients[, "Pr(>|z|)"]`.

``` r

wt <- wald_test(fit_logistic)
print(wt)
```

The likelihood ratio test compares the full model against an
intercept-only model:

``` r

Data_null <- Data_logistic
Data_null$parm.names <- "beta0"
Model_null <- function(parm, Data) {
    p  <- plogis(parm[1])
    p  <- pmin(pmax(p, 1e-10), 1 - 1e-10)
    LL <- sum(dbinom(Data$y, 1, p, log = TRUE))
    list(LP = LL, Dev = -2 * LL, Monitor = LL,
         yhat = rbinom(Data$n, 1, p), parm = parm)
}
fit_null <- LaplaceApproximation(Model_null, 0, Data_null,
                                 Iterations = 200, Method = "SPG",
                                 sir = FALSE)
lrt <- lr_test(fit_logistic, fit_null, df = 2)
print(lrt)
```

Compare against `anova(glm)`:

``` r

glm_null <- glm(y_bin ~ 1, family = binomial)
anova_result <- anova(glm_null, glm_fit, test = "LRT")
knitr::kable(anova_result, digits = 4,
             caption = "glm() likelihood ratio test")
```

The LR statistics and \\p\\-values should agree closely.

#### Information criteria

``` r

ic_lucifer <- freq_ic(fit_logistic, n = n)
ic_glm <- c(AIC = AIC(glm_fit), BIC = BIC(glm_fit))
knitr::kable(data.frame(lucifer = ic_lucifer[c("AIC", "BIC")],
                         glm = ic_glm),
             digits = 2, caption = "Information criteria: lucifer vs glm()")
```

### Poisson regression vs `glm(family = poisson)`

Count data with a log link. The structure mirrors the logistic example.

``` r

set.seed(456)
n <- 250
x1 <- rnorm(n)
true_beta_pois <- c(0.3, 0.5, -0.4)
lambda <- exp(true_beta_pois[1] + true_beta_pois[2] * x1)
y_count <- rpois(n, lambda)
```

``` r

glm_pois <- glm(y_count ~ x1, family = poisson)
knitr::kable(coef(summary(glm_pois)), digits = 4,
             caption = "Reference: glm() Poisson regression")
```

``` r

X_pois <- cbind(1, x1)
Data_pois <- list(
    parm.names = c("beta0", "beta1"),
    mon.names  = "LL",
    X = X_pois, y = y_count, n = n, N = n
)
Model_pois <- function(parm, Data) {
    eta    <- Data$X %*% parm
    lambda <- exp(eta)
    LL     <- sum(dpois(Data$y, lambda, log = TRUE))
    list(LP = LL, Dev = -2 * LL, Monitor = LL,
         yhat = rpois(Data$n, lambda), parm = parm)
}
fit_pois <- LaplaceApproximation(Model_pois, c(0, 0), Data_pois,
                                 Iterations = 200, Method = "SPG",
                                 sir = FALSE)
fs_pois <- freq_summary(fit_pois)
print(fs_pois)
```

``` r

comp_pois <- data.frame(
    glm.Estimate     = coef(summary(glm_pois))[, 1],
    glm.SE           = coef(summary(glm_pois))[, 2],
    Laplace.Estimate = fs_pois$coefficients[, "Estimate"],
    Laplace.SE       = fs_pois$coefficients[, "Std.Error"]
)
knitr::kable(comp_pois, digits = 4,
             caption = "Poisson regression: glm() vs Laplace")
```

``` r

ic_pois_lucifer <- freq_ic(fit_pois, n = n)
ic_pois_glm <- c(AIC = AIC(glm_pois), BIC = BIC(glm_pois))
knitr::kable(data.frame(lucifer = ic_pois_lucifer[c("AIC", "BIC")],
                         glm = ic_pois_glm),
             digits = 2, caption = "Information criteria: lucifer vs glm()")
```

## Nonlinear regression vs `nls()`

Nonlinear regression is where the Bayesian machinery begins to add
genuine value over classical optimizers. We use Michaelis-Menten
kinetics as the running example:

\\y_i = \frac{V\_{\max} \cdot x_i}{K_m + x_i} + \varepsilon_i, \qquad
\varepsilon_i \sim \mathcal{N}(0, \sigma^2)\\

This model has three parameters (\\V\_{\max}\\, \\K_m\\, \\\sigma\\), a
nonlinear mean function, and potentially asymmetric confidence regions
because the parameterization is intrinsically curved [\[1\]](#ref1).
MCMC handles this naturally; [`nls()`](https://rdrr.io/r/stats/nls.html)
can struggle with starting values and convergence.

``` r

set.seed(789)
x_mm <- seq(0.1, 10, length.out = 50)
true_Vmax <- 5.0
true_Km   <- 2.0
true_sigma_mm <- 0.3
y_mm <- true_Vmax * x_mm / (true_Km + x_mm) + rnorm(50, 0, true_sigma_mm)
```

``` r

nls_fit <- nls(y_mm ~ Vmax * x_mm / (Km + x_mm),
               start = list(Vmax = 4, Km = 1))
knitr::kable(coef(summary(nls_fit)), digits = 4,
             caption = "Reference: nls() Michaelis-Menten fit")
```

``` r

Data_mm <- list(
    parm.names = c("Vmax", "log.Km", "log.sigma"),
    mon.names  = "LL",
    x = x_mm, y = y_mm, n = 50, N = 50
)
Model_mm <- function(parm, Data) {
    Vmax  <- parm[1]
    Km    <- exp(parm[2])
    sigma <- exp(parm[3])
    mu    <- Vmax * Data$x / (Km + Data$x)
    LL    <- sum(dnorm(Data$y, mu, sigma, log = TRUE))
    list(LP = LL, Dev = -2 * LL, Monitor = LL,
         yhat = rnorm(Data$n, mu, sigma), parm = parm)
}
fit_mm <- LaplaceApproximation(Model_mm, c(4, log(1), log(0.5)),
                               Data_mm, Iterations = 300,
                               Method = "SPG", sir = FALSE)
fs_mm <- freq_summary(fit_mm)
print(fs_mm)
```

The Vmax estimate from lucifer should match the
[`nls()`](https://rdrr.io/r/stats/nls.html) output closely. The Km and
sigma parameters are on the log scale in lucifer, so a direct comparison
requires back-transformation; the point estimates can be compared via
[`exp()`](https://rdrr.io/r/base/Log.html) of the lucifer output.

``` r

nls_coef <- coef(nls_fit)
comp_nls <- data.frame(
    nls      = c(nls_coef["Vmax"], nls_coef["Km"]),
    Laplace  = c(fs_mm$coefficients["Vmax", "Estimate"],
                 exp(fs_mm$coefficients["log.Km", "Estimate"]))
)
rownames(comp_nls) <- c("Vmax", "Km")
knitr::kable(comp_nls, digits = 4,
             caption = "Michaelis-Menten: nls() vs Laplace (Km back-transformed)")
```

For nonlinear models, profile likelihood curves provide more accurate
confidence intervals than the Wald approximation. The
[`profile_likelihood()`](https://robustecologies.github.io/lucifer/reference/profile_likelihood.md)
function works directly on any lucifer fit object, computing profiles by
fixing each parameter and re-optimizing over nuisance parameters:

``` r

prof_mm <- profile_likelihood(fit_mm, Model_mm, Data_mm)
summary(prof_mm)
plot(prof_mm)
```

The asymmetry ratios for `log.Km` and `log.sigma` are typically further
from 1.0 than for the linear regression case, indicating that the Wald
approximation is less adequate. The profile intervals (solid vertical
lines) may differ visibly from the Wald intervals (dashed teal lines),
particularly for the rate constant \\K_m\\.

Residual diagnostics for the nonlinear fit:

``` r

rd_mm <- freq_residuals(fit_mm, Model_mm, Data_mm)
plot(rd_mm)
```

Data cloning provides an alternative route to profile likelihood curves
via `plot(fit_dc, type = "profile")`:

``` r

fit_dc_mm <- DataCloning(Model_mm, Data_mm,
                         Initial.Values = c(4, log(1), log(0.5)),
                         n_clones = c(1, 2, 4, 8, 16),
                         Algorithm = "NUTS", Iterations = 5000)
plot(fit_dc_mm, type = "profile", Model = Model_mm, Data = Data_mm)
```

## Mixed effects models vs `nlme`

Hierarchical models present the strongest case for the frequentist
bridge. A random-intercept model

\\y\_{ij} = \beta_0 + \beta_1 x\_{ij} + u_j + \varepsilon\_{ij}, \qquad
u_j \sim \mathcal{N}(0, \sigma_u^2), \quad \varepsilon\_{ij} \sim
\mathcal{N}(0, \sigma^2)\\

can be fitted either by integrating out the random effects (marginal
likelihood) or by treating them as parameters. lucifer can do both. We
compare against [`nlme::lme()`](https://rdrr.io/pkg/nlme/man/lme.html)
with `method = "ML"` (not REML, since flat priors target the marginal
ML, not the restricted ML).

``` r

set.seed(321)
n_groups  <- 20
n_per     <- 10
n_total   <- n_groups * n_per
group     <- rep(1:n_groups, each = n_per)
x_mixed   <- rnorm(n_total)
true_b0   <- 2.0
true_b1   <- 1.5
true_su   <- 0.8
true_se   <- 0.5
u         <- rnorm(n_groups, 0, true_su)
y_mixed   <- true_b0 + true_b1 * x_mixed + u[group] + rnorm(n_total, 0, true_se)
```

``` r

df_mixed <- data.frame(y = y_mixed, x = x_mixed, group = factor(group))
lme_fit <- nlme::lme(y ~ x, random = ~ 1 | group, data = df_mixed,
                     method = "ML")
knitr::kable(summary(lme_fit)$tTable, digits = 4,
             caption = "Fixed effects (nlme::lme, method = ML)")
vc <- data.frame(
    Component = c("sigma_u (random intercept)", "sigma_e (residual)"),
    Estimate  = c(as.numeric(nlme::VarCorr(lme_fit)[1, "StdDev"]),
                  lme_fit$sigma),
    True      = c(true_su, true_se)
)
knitr::kable(vc, digits = 4, row.names = FALSE,
             caption = "Variance components: lme vs true values")
```

In lucifer, we write the marginal likelihood by integrating out \\u_j\\
analytically. For the random-intercept model, the marginal distribution
of \\\mathbf{y}\_j = (y\_{j1}, \ldots, y\_{jn_j})^T\\ is multivariate
normal:

\\\mathbf{y}\_j \sim \mathcal{N}\left(\beta_0 \mathbf{1} + \beta_1
\mathbf{x}\_j,\\ \sigma_u^2 \mathbf{1}\mathbf{1}^T + \sigma^2 I\right)\\

``` r

Data_mixed <- list(
    parm.names = c("beta0", "beta1", "log.sigma_u", "log.sigma_e"),
    mon.names  = "LL",
    y = y_mixed, x = x_mixed, group = group,
    n = n_total, N = n_total,
    n_groups = n_groups, n_per = n_per
)

Model_mixed <- function(parm, Data) {
    beta0   <- parm[1]
    beta1   <- parm[2]
    sigma_u <- exp(parm[3])
    sigma_e <- exp(parm[4])

    ## Marginal log-likelihood: sum over groups
    LL <- 0
    for (j in seq_len(Data$n_groups)) {
        idx <- which(Data$group == j)
        nj  <- length(idx)
        mu_j <- beta0 + beta1 * Data$x[idx]
        ## Marginal covariance: sigma_u^2 * 11' + sigma_e^2 * I
        V_j <- sigma_u^2 * matrix(1, nj, nj) + sigma_e^2 * diag(nj)
        diff <- Data$y[idx] - mu_j
        ## Use Woodbury for efficiency (not critical for nj=10)
        cholV <- chol(V_j)
        log_det <- 2 * sum(log(diag(cholV)))
        z <- backsolve(cholV, diff, transpose = TRUE)
        LL <- LL - 0.5 * (nj * log(2 * pi) + log_det + sum(z^2))
    }

    yhat <- beta0 + beta1 * Data$x
    list(LP = LL, Dev = -2 * LL, Monitor = LL,
         yhat = rnorm(Data$n, yhat, sigma_e), parm = parm)
}

fit_mixed <- LaplaceApproximation(Model_mixed, c(0, 0, 0, 0),
                                  Data_mixed, Iterations = 300,
                                  Method = "SPG", sir = FALSE)
fs_mixed <- freq_summary(fit_mixed)
print(fs_mixed)
```

Two points deserve attention. First, `lme()` with `method = "REML"`
estimates variance components differently from `method = "ML"`. The
flat-prior lucifer model targets ML because the log-posterior equals the
full marginal log-likelihood, not the restricted log-likelihood. For
variance component comparison, always use `method = "ML"` on the `lme()`
side. Second, when \\\sigma_u^2\\ is near zero, the MLE can lie on the
boundary of the parameter space (exactly zero), while the log-scale
parameterization in lucifer maps this to \\-\infty\\. Data cloning’s
identifiability diagnostic can detect this boundary issue.

## Survival analysis vs `survival`

Right-censored survival data with a Weibull model demonstrate the bridge
for a model class where the likelihood takes a non-standard form. For
observation \\i\\ with event time \\t_i\\ and censoring indicator
\\\delta_i\\:

\\\ell(\theta) = \sum\_{i=1}^n \left\[\delta_i \log f(t_i \mid \theta) +
(1 - \delta_i) \log S(t_i \mid \theta)\right\]\\

where \\f\\ is the Weibull density and \\S\\ is the survival function.
The Weibull parameterization used by
[`survival::survreg()`](https://rdrr.io/pkg/survival/man/survreg.html)
is the accelerated failure time (AFT) model: \\\log(T) = \mu + \sigma
W\\ where \\W \sim \text{EV}(0, 1)\\ is standard extreme value.

``` r

set.seed(654)
n <- 200
x_surv <- rnorm(n)
## Weibull AFT: log(T) = mu + sigma * W
true_mu_surv    <- 2.0
true_beta_surv  <- 0.5
true_scale_surv <- 0.8
log_t <- true_mu_surv + true_beta_surv * x_surv +
         true_scale_surv * log(rexp(n))
t_event <- exp(log_t)
## Random right censoring
cens_time <- rexp(n, rate = 0.05)
time_obs  <- pmin(t_event, cens_time)
status    <- as.integer(t_event <= cens_time)
knitr::kable(data.frame(Metric = "Censoring rate",
                        Value = round(1 - mean(status), 3)),
             row.names = FALSE)
```

``` r

library(survival)
surv_fit <- survreg(Surv(time_obs, status) ~ x_surv, dist = "weibull")
summary(surv_fit)
```

``` r

Data_surv <- list(
    parm.names = c("mu", "beta", "log.scale"),
    mon.names  = "LL",
    x = x_surv, time = time_obs, status = status,
    n = n, N = n
)

Model_surv <- function(parm, Data) {
    mu    <- parm[1]
    beta  <- parm[2]
    scale <- exp(parm[3])  # Weibull scale = 1/shape in survreg param
    eta   <- mu + beta * Data$x

    ## Standardized residuals
    z <- (log(Data$time) - eta) / scale

    ## Log-likelihood (extreme value / Gumbel minimum distribution)
    ## f(z) = exp(z - exp(z)),  S(z) = exp(-exp(z))
    LL <- sum(Data$status * (z - exp(z) - log(scale) - log(Data$time)) +
              (1 - Data$status) * (-exp(z)))

    list(LP = LL, Dev = -2 * LL, Monitor = LL,
         yhat = rep(0, Data$n), parm = parm)
}

fit_surv <- LaplaceApproximation(Model_surv, c(2, 0.5, log(0.8)),
                                 Data_surv, Iterations = 300,
                                 Method = "SPG", sir = FALSE)
fs_surv <- freq_summary(fit_surv)
print(fs_surv)

## Compare lucifer vs survreg vs true values
mle_surv <- fs_surv$coefficients[, "Estimate"]
survreg_coef <- coef(surv_fit)
surv_compare <- data.frame(
    Parameter = c("mu (Intercept)", "beta (x)", "log(scale)"),
    True      = c(true_mu_surv, true_beta_surv, log(true_scale_surv)),
    lucifer   = as.numeric(mle_surv),
    survreg   = c(as.numeric(survreg_coef), log(surv_fit$scale))
)
knitr::kable(surv_compare, digits = 4, row.names = FALSE,
             caption = "Survival model: lucifer vs survreg vs true values")
```

The advantage of the lucifer approach is extensibility: the user can
replace the Weibull density with any custom hazard function, add
time-varying covariates, or introduce frailties, all within the same
Model interface.
[`survreg()`](https://rdrr.io/pkg/survival/man/survreg.html) is limited
to a fixed set of parametric families.

## When they diverge

The preceding sections demonstrate convergence; this section examines
three instructive scenarios where frequentist and Bayesian answers
genuinely differ.

### Small samples

With \\n = 15\\ and \\p = 4\\ parameters, the asymptotic approximations
underlying both the MLE and the Bernstein-von Mises theorem are
unreliable. The posterior under flat priors is essentially the
likelihood surface, which may be poorly approximated by a normal
distribution. A proper prior, even a weak one like \\\beta_j \sim
\mathcal{N}(0, 100)\\, regularizes the posterior and can improve
finite-sample coverage.

``` r

set.seed(999)
n_small <- 15
x_small <- rnorm(n_small)
y_small <- 1.5 + 0.8 * x_small + rnorm(n_small, 0, 0.5)

## Flat-prior model (frequentist)
## ... (same structure as before)

## Weakly informative prior model (Bayesian)
Model_bayes <- function(parm, Data) {
    beta  <- parm[1:2]
    sigma <- exp(parm[3])
    mu    <- Data$X %*% beta
    LL    <- sum(dnorm(Data$y, mu, sigma, log = TRUE))
    ## Weak prior: beta ~ N(0, 100), log(sigma) ~ N(0, 10)
    LP.prior <- sum(dnorm(beta, 0, 10, log = TRUE)) +
                dnorm(parm[3], 0, sqrt(10), log = TRUE)
    LP <- LL + LP.prior
    list(LP = LP, Dev = -2 * LL, Monitor = LL,
         yhat = rnorm(length(mu), mu, sigma), parm = parm)
}
```

The key insight is that with \\n = 15\\, the flat-prior confidence
intervals undercover (empirical coverage below 95%) because the normal
approximation to the likelihood surface is poor. The Bayesian credible
intervals with a proper prior maintain closer-to-nominal coverage by
shrinking extreme estimates toward the prior mean, which acts as
regularization. This is not a failure of either paradigm but a
difference in what they optimize: the frequentist interval is centered
on the MLE (the best point estimate), while the Bayesian interval
reflects the posterior, which incorporates prior information about
plausible parameter values.

### Boundary parameters

Variance components provide the canonical example. When the true
\\\sigma_u^2\\ is exactly zero (no group-level variation), the MLE lies
on the boundary of the parameter space \\\[0, \infty)\\. The
Bernstein-von Mises theorem assumes the true parameter is in the
interior of the parameter space and does not apply at boundaries.

The MLE \\\hat\sigma_u^2 = 0\\ is a valid point estimate, but the Wald
confidence interval \\0 \pm 1.96 \cdot \text{SE}\\ extends into negative
territory, which is physically meaningless. A Bayesian analysis with a
proper prior on \\\sigma_u\\ naturally constrains the posterior to
positive values. The 95% credible interval \\\[0.001, 0.15\]\\ (for
example) is coherent even though the mass is concentrated near zero.

Data cloning detects boundary issues through the eigenvalue diagnostic:
if the eigenvalue ratio \\\lambda_S(K)\\ does not decrease at rate
\\1/K\\ for a variance component, that component may not be identifiable
or may be near a boundary.

### Model misspecification

When the true data-generating process is not a member of the assumed
model family, both the MLE and the posterior concentrate on the
pseudo-true value, defined as the parameter value that minimizes the
Kullback-Leibler divergence from the true distribution to the model
family:

\\\theta^\* = \arg\min\_\theta \\ \text{KL}\\\left(p_0 \\\\\\
p\_\theta\right)\\

Both paradigms agree on the point estimate. However, the standard errors
and interval widths differ. The Hessian-based standard error assumes the
model is correctly specified and produces intervals that are too narrow
when it is not. The sandwich estimator \\\hat{V}\_S = J^{-1} \hat{I}
J^{-1}\\ (available in lucifer via `CovEst = "Sandwich"`) accounts for
misspecification and gives asymptotically valid frequentist coverage
[\[10\]](#ref10).

The Bayesian posterior covariance, by contrast, has no built-in
misspecification correction. Under misspecification, the posterior
concentrates at the same rate as the correctly specified case, producing
credible intervals that are typically too narrow and have below-nominal
coverage. Calibrating Bayesian inference under misspecification requires
explicitly modeling the discrepancy, which is beyond the scope of this
vignette.

## Coverage analysis

The
[`coverage_sim()`](https://robustecologies.github.io/lucifer/reference/coverage_sim.md)
function provides a principled way to verify that confidence and
credible intervals achieve their nominal coverage. For each simulation
replicate, data are generated from a known DGP, the model is fitted with
each method, intervals are extracted, and coverage of the true parameter
values is checked.

``` r

dgp <- function() {
    n <- 100
    x <- rnorm(n)
    y <- 1.5 + 0.8 * x + rnorm(n, 0, 0.5)
    X <- cbind(1, x)
    Data <- list(
        parm.names = c("beta0", "beta1", "log.sigma"),
        mon.names  = "LL",
        X = X, y = y, n = n, N = n
    )
    Model <- function(parm, Data) {
        beta  <- parm[1:2]
        sigma <- exp(parm[3])
        mu    <- Data$X %*% beta
        LL    <- sum(dnorm(Data$y, mu, sigma, log = TRUE))
        list(LP = LL, Dev = -2 * LL, Monitor = LL,
             yhat = rnorm(length(mu), mu, sigma), parm = parm)
    }
    list(Model = Model, Data = Data,
         Initial.Values = c(0, 0, 0))
}

methods <- list(
    laplace = function(M, D, IV) {
        LaplaceApproximation(M, IV, D, Iterations = 200,
                             Method = "SPG", sir = FALSE)
    }
)

result <- coverage_sim(dgp, methods,
                       true_parms = c(beta0 = 1.5, beta1 = 0.8,
                                      log.sigma = log(0.5)),
                       n_sims = 500, seed = 666, verbose = FALSE)
summary(result)
plot(result)
```

With \\n = 100\\ and a correctly specified model, all methods should
achieve approximately 95% coverage. The empirical coverage will
fluctuate around the nominal level with a standard error of \\\sqrt{0.95
\times 0.05 / 500} \approx 0.01\\, so values between 93% and 97% are
expected.

## Practical guidelines

The frequentist bridge in lucifer offers three routes to the MLE, each
with distinct strengths and trade-offs.

**Laplace approximation** is the fastest path. It requires a single
optimization run (typically 10-50 iterations) and produces the mode and
Hessian-based covariance in seconds. Use it for standard GLMs and
nonlinear models where the likelihood is unimodal and the sample size is
moderate to large. The main limitation is that it provides only a point
estimate and a normal approximation to the sampling distribution; it
cannot detect multimodality or heavy tails.

**Flat-prior MCMC** is the most pedagogically transparent route. The
user writes `LP = LL` in the Model function and runs lucifer’s MCMC
algorithms as usual. The posterior mean is the MLE, the posterior
covariance is the inverse Fisher information, and all MCMC diagnostics
(R-hat, ESS, trace plots) apply. Use it when you want to verify the
Laplace results, when the likelihood surface is multimodal, or as a
stepping stone toward proper Bayesian analysis: once you confirm that
flat-prior MCMC reproduces [`lm()`](https://rdrr.io/r/stats/lm.html) or
[`glm()`](https://rdrr.io/r/stats/glm.html), you can introduce
informative priors and observe how the posterior changes.

**Data cloning** is the most principled approach. It retains the prior
(avoiding the non-standard `LP = LL` construct), provides automatic
identifiability testing, profile likelihood curves, and convergence
diagnostics. Use it for complex hierarchical models where
identifiability is a concern, or when you need profile-based confidence
intervals for nonlinear parameters. The cost is computational: it runs
the MCMC sampler multiple times at increasing clone sizes.

For users transitioning from frequentist to Bayesian inference, the
recommended workflow is:

1.  Write the flat-prior Model (`LP = LL`) and verify against the known
    frequentist solution (`lm`, `glm`, `nls`)
2.  Introduce weakly informative priors (`LP = LL + LP.prior`) and
    compare the posterior to the flat-prior result
3.  For complex models, use data cloning to obtain frequentist
    benchmarks and identifiability diagnostics
4.  Once comfortable with the Bayesian framework, use proper priors and
    full posterior inference

## References

**\[1\]** Bates, D. M. and Watts, D. G. (1988). *Nonlinear Regression
Analysis and Its Applications*. Wiley. ISBN
[978-0471816430](https://www.wiley.com/en-us/Nonlinear+Regression+Analysis+and+Its+Applications-p-9780471816430).

**\[2\]** Berndt, E. K., Hall, B. H., Hall, R. E. and Hausman, J. A.
(1974). Estimation and inference in nonlinear structural models. *Annals
of Economic and Social Measurement*, 3(4), 653-665.

**\[3\]** Campbell, D. and Lele, S. R. (2014). An ANOVA test for
parameter estimability using data cloning with application to
statistical inference for dynamic systems. *Computational Statistics and
Data Analysis*, 70, 257-267. DOI:
[10.1016/j.csda.2013.09.013](https://doi.org/10.1016/j.csda.2013.09.013).

**\[4\]** Engle, R. F. (1984). Wald, likelihood ratio, and Lagrange
multiplier tests in econometrics. In *Handbook of Econometrics*, Vol. 2,
775-826. North-Holland.

**\[5\]** Freedman, D. A. (1999). On the Bernstein-von Mises theorem
with infinite-dimensional parameters. *Annals of Statistics*, 27(4),
1119-1140. DOI:
[10.1214/aos/1017938917](https://doi.org/10.1214/aos/1017938917).

**\[6\]** Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B.,
Vehtari, A. and Rubin, D. B. (2013). *Bayesian Data Analysis*, 3rd
ed. Chapman & Hall/CRC. ISBN
[978-1439840955](https://www.routledge.com/Bayesian-Data-Analysis/Gelman-Carlin-Stern-Dunson-Vehtari-Rubin/p/book/9781439840955).

**\[7\]** Lele, S. R., Dennis, B. and Lutscher, F. (2007). Data cloning:
easy maximum likelihood estimation for complex ecological models using
Bayesian Markov chain Monte Carlo methods. *Ecology Letters*, 10(7),
551-563. DOI:
[10.1111/j.1461-0248.2007.01047.x](https://doi.org/10.1111/j.1461-0248.2007.01047.x).

**\[8\]** van der Vaart, A. W. (1998). *Asymptotic Statistics*.
Cambridge University Press. ISBN
[978-0521784504](https://doi.org/10.1017/CBO9780511802256).

**\[9\]** Walker, A. M. (1969). On the asymptotic behaviour of posterior
distributions. *Journal of the Royal Statistical Society: Series B*,
31(1), 80-88. DOI:
[10.1111/j.2517-6161.1969.tb00796.x](https://doi.org/10.1111/j.2517-6161.1969.tb00796.x).

**\[10\]** White, H. (1982). Maximum likelihood estimation of
misspecified models. *Econometrica*, 50(1), 1-25. DOI:
[10.2307/1912526](https://doi.org/10.2307/1912526).

**\[11\]** Wilks, S. S. (1938). The large-sample distribution of the
likelihood ratio for testing composite hypotheses. *Annals of
Mathematical Statistics*, 9(1), 60-62. DOI:
[10.1214/aoms/1177732360](https://doi.org/10.1214/aoms/1177732360).
