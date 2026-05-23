# Gibbs sampling in lucifer

  

## Introduction

The Gibbs sampler is one of the foundational algorithms in computational
Bayesian statistics. Unlike generic Metropolis-Hastings schemes that
propose moves in the full parameter space and accept or reject them, the
Gibbs sampler decomposes the joint posterior into a sequence of
univariate (or block) full-conditional distributions and samples from
each in turn. When these full conditionals are available in closed form,
the algorithm achieves an acceptance rate of exactly one, which
eliminates the need for tuning proposal distributions and typically
produces chains with lower autocorrelation than random-walk
alternatives.

The classical limitation of Gibbs sampling is that it requires the user
to derive the full-conditional distributions analytically, a task that
demands specialized knowledge of conjugate families and becomes
increasingly difficult as models grow in complexity. Software like JAGS
and BUGS circumvents this by parsing a restricted model specification
language and automatically identifying conjugacy; lucifer, by contrast,
accepts arbitrary R functions as model specifications, which makes
closed-form derivation intractable at the language level.

This vignette describes lucifer’s solution: **automatic full-conditional
generation (auto-FC)**, a numerical approach that samples from full
conditionals without closed-form expressions. The strategy uses
univariate slice sampling for continuous parameters and exhaustive
enumeration for discrete parameters, producing exact draws from the
target conditionals in both cases. We also show how missing data can be
treated as latent variables within this framework, recovering the
data-augmentation approach of Tanner and Wong [\[5\]](#ref5) without any
special machinery. The core idea of sampling from a multivariate
distribution by cycling through its univariate conditionals was first
proposed by Turchin [\[6\]](#ref6) for Monte Carlo integration, though
the technique did not receive widespread attention until Geman and Geman
[\[2\]](#ref2) independently introduced it under the name “Gibbs
sampling” for image restoration. Its adoption as a general-purpose tool
for Bayesian computation followed the work of Gelfand and Smith
[\[1\]](#ref1).

  

## Theory

### Full-conditional distributions

Consider a posterior \\p(\theta_1, \ldots, \theta_K \mid y)\\ over \\K\\
parameters. The full conditional of parameter \\\theta_j\\ is the
distribution obtained by fixing all other parameters at their current
values:

\\p(\theta_j \mid \theta\_{-j}, y) \propto p(\theta_1, \ldots, \theta_K
\mid y),\\

where \\\theta\_{-j} = (\theta_1, \ldots, \theta\_{j-1}, \theta\_{j+1},
\ldots, \theta_K)\\ is treated as a constant. Since the normalizing
constant of the joint posterior cancels in the full conditional,
\\p(\theta_j \mid \theta\_{-j}, y)\\ is determined entirely by the
unnormalized log-posterior as a function of \\\theta_j\\ alone. This is
precisely what the lucifer Model function evaluates.

The Hammersley-Clifford theorem guarantees that a strictly positive
joint distribution is uniquely determined by its full-conditional
distributions. The Gibbs sampler exploits this by constructing a Markov
chain that cycles through the full conditionals: each transition
replaces one component \\\theta_j\\ with a draw from \\p(\theta_j \mid
\theta\_{-j}, y)\\, which preserves the joint posterior as the
stationary distribution.

### The Gibbs sampler

**Systematic scan.** Initialize \\\theta^{(0)}\\. At iteration \\t\\,
for \\j = 1, \ldots, K\\:

\\\theta_j^{(t)} \sim p\\\left(\theta_j \mid \theta_1^{(t)}, \ldots,
\theta\_{j-1}^{(t)}, \theta\_{j+1}^{(t-1)}, \ldots, \theta_K^{(t-1)},
y\right).\\

**Random scan.** At each iteration, permute the indices \\\\1, \ldots,
K\\\\ uniformly at random and update in that order. This variant
satisfies detailed balance and is theoretically cleaner for convergence
analysis; lucifer uses random scan in auto-FC mode.

### Convergence

Under mild regularity conditions (the joint density is positive on its
support), the Gibbs sampler is ergodic and the chain converges to the
target distribution regardless of initialization. The rate of
convergence depends on the correlation structure of the posterior: high
pairwise correlations between parameters slow mixing because each
conditional update makes only a small adjustment. Reparameterization and
blocking (updating correlated groups jointly) can substantially improve
convergence in such cases.

When full conditionals are sampled exactly (conjugate Gibbs or the
auto-FC methods described below), every proposed move is accepted and
the chain satisfies the conditions of the Gibbs sampler theorem
directly. There is no approximation error in the stationary
distribution; the only concern is mixing speed.

  

## Automatic full conditionals in lucifer

### The problem

Deriving conjugate full conditionals from an arbitrary R Model function
is computationally intractable. The Model is a black box: it takes a
parameter vector and a data list, computes the unnormalized
log-posterior, and returns it alongside deviance, monitors, and
posterior predictive draws. No symbolic representation of the prior or
likelihood is available. This rules out the approach taken by JAGS/BUGS,
which relies on parsing a domain-specific language to identify conjugate
pairs.

### Slice sampling for continuous parameters

The solution for continuous parameters is **univariate slice sampling**
[\[4\]](#ref4). Given the current state \\\theta\_{-j}\\ and the
unnormalized density \\f(\theta_j) \propto p(\theta_j \mid \theta\_{-j},
y)\\, slice sampling introduces an auxiliary variable \\u \sim
\text{Uniform}(0, f(\theta_j^{(t-1)}))\\ and samples \\\theta_j^{(t)}\\
uniformly from the “slice” \\\\x : f(x) \> u\\\\. The resulting draw is
exactly distributed according to \\f\\, with acceptance probability one.

In practice, the slice \\\\x : f(x) \> u\\\\ is not known analytically.
The **stepping-out** procedure [\[4\]](#ref4) approximates it: starting
from the current value, an interval \\\[L, R\]\\ is expanded in steps of
width \\w\\ until both endpoints fall outside the slice, then a point is
drawn uniformly from \\\[L, R\]\\ and the interval is shrunk toward the
current value if the draw falls outside the slice. The procedure is
exact regardless of the choice of \\w\\; the width only affects
efficiency (too small means many expansions, too large means many
shrinkage steps). In auto-FC mode, \\w\\ is taken from the tuning
vector, which defaults to \\2.381^2 / K\\ (the optimal Metropolis
scaling, repurposed as a reasonable initial slice width).

**Correctness.** Slice sampling is an auxiliary-variable method. The
joint distribution of \\(\theta_j, u)\\ is uniform on the region under
the density curve; marginalizing over \\u\\ recovers \\f(\theta_j)\\
exactly. The stepping-out and shrinkage procedures preserve this joint
distribution by construction [\[4\]](#ref4).

### Enumeration for discrete parameters

For a discrete parameter \\\theta_j\\ with finite support
\\\mathcal{S}\_j = \\s_1, \ldots, s_M\\\\, the full conditional is a
categorical distribution:

\\p(\theta_j = s_k \mid \theta\_{-j}, y) = \frac{\exp\[\log p(s_k,
\theta\_{-j} \mid y)\]}{\sum\_{m=1}^M \exp\[\log p(s_m, \theta\_{-j}
\mid y)\]}.\\

This requires \\M\\ evaluations of the log-posterior, one for each
element of \\\mathcal{S}\_j\\. For Bernoulli indicators (\\M = 2\\), the
cost is exactly two Model evaluations per parameter per iteration. The
probabilities are computed in log space using the log-sum-exp trick for
numerical stability.

This approach is standard in Gibbs samplers for variable selection
models such as stochastic search variable selection (SSVS)
[\[3\]](#ref3), where binary inclusion indicators \\\gamma_j \in \\0,
1\\\\ select which predictors enter the model.

### The combined strategy

Auto-FC in lucifer applies the appropriate method to each parameter
based on its type:

- **Continuous parameters** (default): stepping-out slice sampling.
- **Discrete parameters** (declared via `Data$dparm`): exhaustive
  enumeration.

The classification is specified by the user through two optional
elements of the `Data` list:

- `Data$dparm`: an integer vector of parameter indices that are
  discrete.
- `Data$dsupport`: a list of support vectors, parallel to `dparm`. Each
  element is a numeric vector of the possible values for the
  corresponding discrete parameter. If `NULL`, all discrete parameters
  default to \\\\0, 1\\\\ (Bernoulli).

At each Gibbs iteration, the auto-FC closure performs a random-scan
sweep over all \\K\\ parameters, dispatching to slice sampling or
enumeration as appropriate. The current log-posterior value is passed
between updates to avoid redundant Model evaluations.

``` r

# All continuous (default): just omit Specs
fit <- lucifer(Model, Data, IV, Algorithm = "Gibbs")

# Mixed continuous + discrete (e.g., SSVS indicators)
Data$dparm <- pos.Gamma             # which parameters are discrete
Data$dsupport <- NULL               # default: {0, 1} (Bernoulli)
fit <- lucifer(Model, Data, IV, Algorithm = "Gibbs")

# Discrete with non-binary support
Data$dparm <- c(5, 6)
Data$dsupport <- list(0:3, 0:5)     # parameter 5 takes values 0-3, etc.
fit <- lucifer(Model, Data, IV, Algorithm = "Gibbs")
```

  

## Missing data as latent variables

A natural extension of Gibbs sampling is the treatment of missing
observations as latent variables. JAGS, for example, detects `NA` values
in the data and automatically includes them in the sampling scheme,
drawing imputed values from their full conditional \\p(y_i^{\text{miss}}
\mid \theta, y^{\text{obs}})\\ at each iteration. This is a special case
of the data augmentation framework of Tanner and Wong [\[5\]](#ref5),
which treats the complete-data posterior as the augmented target and
alternates between sampling parameters conditional on imputed data and
imputing data conditional on parameters.

In lucifer, the same approach works without any special infrastructure.
The user augments the parameter vector to include the missing values,
and the Model function fills in the missing entries before computing the
log-posterior. Auto-FC then samples the imputed values via slice
sampling, exactly as it does for any other continuous parameter. The
full conditional of an imputed observation \\y_i^{\text{miss}}\\ is
determined by the likelihood term for observation \\i\\ and any prior on
\\y_i\\; since the Model function evaluates both, the slice sampler
targets the correct conditional automatically.

The pattern is:

1.  Record which observations are missing and their positions.
2.  Append one parameter per missing value to `parm.names`.
3.  In the Model function, replace `NA`s with the corresponding
    parameter values before computing the likelihood.

The resulting Gibbs sampler performs proper Bayesian imputation: each
posterior draw of \\\theta\\ is accompanied by a draw of
\\y^{\text{miss}}\\ from its predictive distribution, and inferences
about \\\theta\\ marginalize over the uncertainty in the missing data.

  

## Example 1: linear regression

A simple linear regression with simulated data and known parameters
demonstrates the basic usage of auto-FC.

``` r

library(lucifer)
library(ggplot2)
set.seed(42)

## Simulate data
n <- 200
true_beta <- c(2.0, -1.5, 0.8)
true_sigma <- 1.5
X <- cbind(1, matrix(rnorm(n * 2), n, 2))
y <- rnorm(n, X %*% true_beta, true_sigma)

## Data list
Data <- list(
    y = y, X = X, n = n,
    parm.names = c("beta0", "beta1", "beta2", "sigma"),
    mon.names = "LP",
    pos.beta = 1:3,
    pos.sigma = 4
)

## Model
Model <- function(parm, Data) {
    beta  <- parm[Data$pos.beta]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    mu <- Data$X %*% beta
    LL <- sum(dnorm(Data$y, mu, sigma, log = TRUE))
    LP <- LL + sum(dnorm(beta, 0, 10, log = TRUE)) +
        dhalfcauchy(sigma, 25, log = TRUE)
    list(LP = LP, Dev = -2 * LL, Monitor = LP,
         yhat = rnorm(Data$n, mu, sigma), parm = parm)
}

## Fit with auto-FC Gibbs (no Specs)
fit <- lucifer(Model, Data,
    Initial.Values = c(rep(0, 3), 1),
    Iterations = 50000, Status = 10000,
    Thinning = 50, Algorithm = "Gibbs",
    Chains = 3)

## Compare posteriors to truth
post <- fit$Posterior1
true_vals <- c(true_beta, true_sigma)
pnames <- Data$parm.names
results <- data.frame(
    Parameter = pnames,
    True = true_vals,
    Mean = sapply(seq_along(pnames), function(i) mean(post[, i])),
    Lower = sapply(seq_along(pnames), function(i) quantile(post[, i], 0.025)),
    Upper = sapply(seq_along(pnames), function(i) quantile(post[, i], 0.975)),
    row.names = NULL
)
results$Coverage <- ifelse(true_vals >= results$Lower & true_vals <= results$Upper,
                           "Yes", "No")
knitr::kable(results, digits = 3, align = "lrrrrr",
             col.names = c("Parameter", "True", "Mean", "2.5%", "97.5%", "Covered"),
             caption = "Posterior summaries vs. true values (linear regression)")
```

All true values should fall within the 95% credible intervals. The
acceptance rate is exactly 1.

``` r

## Posterior density vs. true values via lucifer's native ground-truth grid
gt_lin <- setNames(true_vals, pnames)
plot(fit, ground_truth = gt_lin)
```

  

## Example 2: linear regression with missing data

We now introduce 20% missing values in the response variable and treat
them as latent parameters.

``` r

set.seed(42)

## Same data as Example 1, but introduce missing values
n <- 200
true_beta <- c(2.0, -1.5, 0.8)
true_sigma <- 1.5
X <- cbind(1, matrix(rnorm(n * 2), n, 2))
y_complete <- rnorm(n, X %*% true_beta, true_sigma)

## Randomly mask 20% of observations
missing_idx <- sort(sample(n, 40))
y_obs <- y_complete
y_obs[missing_idx] <- NA
n_miss <- length(missing_idx)

## Augment parameter vector: original params + imputed y values
parm.names <- c("beta0", "beta1", "beta2", "sigma",
                paste0("y_imp_", missing_idx))
n_parm <- length(parm.names)
pos.beta  <- 1:3
pos.sigma <- 4
pos.imp   <- 5:n_parm

Data <- list(
    y_obs = y_obs, X = X, n = n,
    missing_idx = missing_idx,
    n_miss = n_miss,
    parm.names = parm.names,
    mon.names = "LP",
    pos.beta = pos.beta,
    pos.sigma = pos.sigma,
    pos.imp = pos.imp
)

## Model: fill in missing y from parameters before computing likelihood
Model <- function(parm, Data) {
    beta  <- parm[Data$pos.beta]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma

    ## Reconstruct complete y
    y <- Data$y_obs
    y[Data$missing_idx] <- parm[Data$pos.imp]

    mu <- Data$X %*% beta
    LL <- sum(dnorm(y, mu, sigma, log = TRUE))
    LP <- LL + sum(dnorm(beta, 0, 10, log = TRUE)) +
        dhalfcauchy(sigma, 25, log = TRUE)
    list(LP = LP, Dev = -2 * LL, Monitor = LP,
         yhat = rnorm(Data$n, mu, sigma), parm = parm)
}

## Initial values: 0 for betas, 1 for sigma, mean(y_obs) for imputed
iv <- c(rep(0, 3), 1, rep(mean(y_obs, na.rm = TRUE), n_miss))

fit_miss <- lucifer(Model, Data,
    Initial.Values = iv,
    Iterations = 1000, Status = 100,
    Thinning = 1, Algorithm = "Gibbs",
    Chains = 3)

## Check regression parameters (should match Example 1 closely)
post <- fit_miss$Posterior1
true_reg <- c(true_beta, true_sigma)
reg_results <- data.frame(
    Parameter = parm.names[1:4],
    True = true_reg,
    Mean = sapply(1:4, function(i) mean(post[, i])),
    SD = sapply(1:4, function(i) sd(post[, i])),
    Lower = sapply(1:4, function(i) quantile(post[, i], 0.025)),
    Upper = sapply(1:4, function(i) quantile(post[, i], 0.975)),
    row.names = NULL
)
knitr::kable(reg_results, digits = 3, align = "lrrrrr",
             col.names = c("Parameter", "True", "Mean", "SD", "2.5%", "97.5%"),
             caption = "Posterior summaries with 20% missing data")

## Check imputed values against truth
imp_means <- colMeans(post[, pos.imp])
imp_rmse <- sqrt(mean((imp_means - y_complete[missing_idx])^2))
knitr::kable(data.frame(Metric = "RMSE (imputed vs. true)", Value = round(imp_rmse, 3)),
             align = "lr", col.names = c("Metric", "Value"),
             caption = "Imputation accuracy")
```

The regression coefficients recover the true values with slightly wider
credible intervals than the complete-data case (reflecting the
additional uncertainty from missing data). The imputed values cluster
around their true values, with posterior SDs reflecting the predictive
uncertainty.

The dedicated
[`plot_imputed()`](https://robustecologies.github.io/lucifer/reference/plot_imputed.md)
function visualises the imputation result directly in the time-series
space: observed values are drawn as a connected line with filled
circles, posterior median imputations as hollow diamonds at the missing
positions, 95% credible intervals as vertical bars, and ground truth
(when supplied) as small X markers. The function auto-detects parameters
whose names contain `imp` or `miss` in the posterior column names, so
for the conventional `y_imp_*` naming used here no explicit indicator
argument is required.

``` r

plot_imputed(fit_miss, data = y_obs, ground_truth = y_complete)
```

Two additional `Style` options are available when the missingness
pattern is structured. `Style = "ribbon"` converts the per-position
credible intervals into a continuous shaded band, which is the most
readable choice when the missing fraction is concentrated in contiguous
gaps:

``` r

plot_imputed(fit_miss, data = y_obs, Style = "ribbon",
             ground_truth = y_complete)
```

`Style = "draws"` overlays a spaghetti of posterior samples threading
through the missing positions, exposing the joint correlation structure
between nearby imputations that the marginal credible intervals collapse
away:

``` r

plot_imputed(fit_miss, data = y_obs, Style = "draws", n_draws = 80,
             ground_truth = y_complete)
```

### Example 2b: multivariate time series with missing observations

The same data-augmentation strategy extends naturally to multivariate
time series. Suppose we observe three independent AR(1) processes with
their own autoregressive coefficients and innovation variances, and that
two of the three series carry scattered missing observations while the
third remains complete. The augmented parameter vector includes one
latent variable per `NA` cell of the response matrix; auto-FC Gibbs
samples those latent variables jointly with the dynamic parameters via
slice sampling, exactly as in the univariate case.
[`plot_imputed()`](https://robustecologies.github.io/lucifer/reference/plot_imputed.md)
recognises that the input is a matrix, infers series labels from
[`colnames()`](https://rdrr.io/r/base/colnames.html), and produces one
panel per series so the imputation result remains visually
self-contained regardless of dimensionality.

``` r

set.seed(13)
Tt <- 60; J <- 3
true_phi   <- c(0.7, -0.4, 0.85)
true_sigma <- c(0.4, 0.6, 0.3)

## Simulate three independent AR(1) processes
Y_complete <- matrix(0, Tt, J)
for (j in 1:J) {
    Y_complete[1, j] <- rnorm(1, 0, true_sigma[j])
    for (t in 2:Tt) {
        Y_complete[t, j] <- true_phi[j] * Y_complete[t - 1, j] +
            rnorm(1, 0, true_sigma[j])
    }
}
colnames(Y_complete) <- c("series A", "series B", "series C")

## Mask: 8 NAs in series A, 5 NAs in series C, series B complete
Y_obs <- Y_complete
Y_obs[sort(sample(2:Tt,  8)), 1] <- NA
Y_obs[sort(sample(2:Tt,  5)), 3] <- NA

## NA positions in column-major order (lucifer convention)
na_flat <- which(is.na(as.vector(Y_obs)))
n_imp   <- length(na_flat)
n_imp
```

``` r

parm.names <- c(paste0("phi[",   1:J, "]"),
                paste0("sigma[", 1:J, "]"),
                paste0("y_imp[", na_flat, "]"))

Data_mv <- list(
    Y          = Y_obs, Tt = Tt, J = J,
    na_flat    = na_flat,
    parm.names = parm.names,
    mon.names  = "LP",
    pos.phi    = 1:J,
    pos.sigma  = (J + 1):(2 * J),
    pos.imp    = (2 * J + 1):(2 * J + n_imp)
)

Model_mv <- function(parm, Data) {
    phi   <- parm[Data$pos.phi]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma

    ## Reconstruct the complete matrix by filling NAs from the
    ## augmented parameter vector
    Yfull <- Data$Y
    Yfull[Data$na_flat] <- parm[Data$pos.imp]

    LL <- 0
    for (j in 1:Data$J) {
        LL <- LL + sum(dnorm(Yfull[-1, j],
                             phi[j] * Yfull[-Data$Tt, j],
                             sigma[j], log = TRUE))
    }
    LP <- LL +
        sum(dnorm(phi, 0, 1, log = TRUE)) +
        sum(dhalfcauchy(sigma, 5, log = TRUE))
    list(LP = LP, Dev = -2 * LL, Monitor = LP,
         yhat = as.vector(Yfull), parm = parm)
}

iv_mv  <- c(rep(0, J), rep(0.5, J), rep(0, n_imp))
fit_mv <- lucifer(Model_mv, Data_mv, iv_mv,
                  Iterations = 1000, Status = 1000,
                  Thinning = 5, Algorithm = "Gibbs", Chains = 1)
```

The faceted layout makes it immediately clear which series carry missing
observations, while series B (complete) is rendered without any
imputation overlay. The subtitle reports the total counts of observed
and imputed cells across all series, and the matching to `NA` positions
in column-major order is handled automatically by
[`plot_imputed()`](https://robustecologies.github.io/lucifer/reference/plot_imputed.md).

``` r

plot_imputed(fit_mv, data = Y_obs, ground_truth = Y_complete)
```

For the multivariate case the `Style = "ribbon"` variant is particularly
informative when the missingness is sparse: each panel carries a
continuous shaded band only at the imputed positions, leaving the
observed sections as a clean polyline.

``` r

plot_imputed(fit_mv, data = Y_obs, Style = "ribbon",
             ground_truth = Y_complete)
```

  

## Example 3: VAR(p) with SSVS

This example adapts the VAR(p)-SSVS model from the Examples vignette.
Stochastic search variable selection [\[3\]](#ref3) places
spike-and-slab priors on VAR autoregressive coefficients, using binary
inclusion indicators \\\Gamma\_{i,k,p} \in \\0, 1\\\\ to determine which
cross-series dependencies are active. The indicator variables are
discrete, making this model a natural test case for the mixed
continuous/discrete auto-FC strategy.

### Model specification

For \\J\\ time series of length \\T\\ with \\P\\ autoregressive lags at
positions \\L_1, \ldots, L_P\\:

\\Y\_{t,j} \sim \mathcal{N}(\mu\_{t,j},\\ \sigma_j^2), \qquad t = L_P +
1, \ldots, T,\quad j = 1, \ldots, J\\ \\\mu\_{t,j} = \alpha_j +
\sum\_{p=1}^{P} \sum\_{i=1}^{J} \Gamma\_{i,j,p}\\ \Phi\_{i,j,p}\\
Y\_{t - L_p,\\ i}\\ \\\alpha_j \sim \mathcal{N}(0, 1000)\\
\\\Gamma\_{i,k,p} \sim \text{Bern}(0.5)\\ \\\Phi\_{i,k,p} \mid
\Gamma\_{i,k,p} \sim (1 - \Gamma\_{i,k,p})\\\mathcal{N}(0,\\ \tau_0^2) +
\Gamma\_{i,k,p}\\\mathcal{N}(0,\\ \tau_1^2)\\ \\\sigma_j \sim
\text{Half-Cauchy}(25)\\

The spike variance \\\tau_0^2 = 0.01\\ concentrates inactive
coefficients near zero; the slab variance \\\tau_1^2 = 10\\ allows
active coefficients to take substantial values. The product
\\\Gamma\_{i,k,p} \cdot \Phi\_{i,k,p}\\ in the mean function ensures
that when an indicator is zero, the corresponding coefficient has no
effect on the likelihood, regardless of its value.

### Simulation

``` r

library(lucifer)
library(ggplot2)
set.seed(42)

J <- 3
TT <- 300
P <- 1
L <- 1

## True parameters
true_alpha <- c(0.3, -0.2, 0.1)
true_Gamma <- matrix(c(
    1, 0, 0,
    0, 1, 0,
    0, 0, 1), J, J, byrow = TRUE)   # only diagonal (own-lag) active
true_Phi <- matrix(c(
    0.6,  0.0, 0.0,
    0.0, -0.4, 0.0,
    0.0,  0.0, 0.5), J, J, byrow = TRUE)
true_sigma <- c(0.8, 1.0, 0.6)

## Check stationarity
ev <- eigen(true_Gamma * true_Phi)$values
stopifnot(all(Mod(ev) < 1))

## Simulate VAR(1)
Y <- matrix(0, TT, J)
Y[1, ] <- true_alpha
for (t in 2:TT) {
    Y[t, ] <- true_alpha +
        (true_Gamma * true_Phi) %*% Y[t - 1, ] +
        rnorm(J, 0, true_sigma)
}
Y <- Y[51:TT, ]    # discard burn-in
TT_eff <- nrow(Y)
```

### Data list with discrete parameter declaration

``` r

parm.names <- as.parm.names(list(
    alpha = rep(0, J),
    Gamma = matrix(0, J, J),
    Phi   = matrix(0, J, J),
    sigma = rep(0, J)))
pos.alpha <- grep("alpha", parm.names)
pos.Gamma <- grep("Gamma", parm.names)
pos.Phi   <- grep("Phi",   parm.names)
pos.sigma <- grep("sigma", parm.names)

MyData <- list(
    J = J, L = L, P = P, TT = TT_eff, Y = Y,
    parm.names = parm.names,
    mon.names  = "LP",
    pos.alpha = pos.alpha,
    pos.Gamma = pos.Gamma,
    pos.Phi   = pos.Phi,
    pos.sigma = pos.sigma,
    ## Declare Gamma indicators as discrete Bernoulli
    dparm = pos.Gamma,
    dsupport = NULL       # default: {0, 1}
)
```

### Model function

``` r

Model <- function(parm, Data) {
    alpha <- parm[Data$pos.alpha]
    Gamma <- matrix(parm[Data$pos.Gamma], Data$J, Data$J)
    Phi   <- matrix(parm[Data$pos.Phi],   Data$J, Data$J)
    Phi.Sigma <- ifelse(Gamma == 1, sqrt(10), 0.1)
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma

    ## Log-prior
    alpha.prior <- sum(dnormv(alpha, 0, 1000, log = TRUE))
    Gamma.prior <- sum(dbern(Gamma, 0.5, log = TRUE))
    Phi.prior   <- sum(dnorm(Phi, 0, Phi.Sigma, log = TRUE))
    sigma.prior <- sum(dhalfcauchy(sigma, 25, log = TRUE))

    ## Log-likelihood
    mu <- matrix(alpha, Data$TT - 1, Data$J, byrow = TRUE) +
        Data$Y[1:(Data$TT - 1), ] %*% t(Gamma * Phi)
    Sigma_mat <- matrix(sigma, Data$TT - 1, Data$J, byrow = TRUE)
    LL <- sum(dnorm(Data$Y[2:Data$TT, ], mu, Sigma_mat, log = TRUE))

    LP <- LL + alpha.prior + Gamma.prior + Phi.prior + sigma.prior
    yhat <- rnorm(length(mu), mu, Sigma_mat)
    list(LP = LP, Dev = -2 * LL, Monitor = LP,
         yhat = yhat, parm = parm)
}
```

### Fitting

``` r

## Initial values: all Gamma = 1, Phi near zero, sigma = 1
IV <- c(rep(0, J), rep(1, J * J), rep(0.1, J * J), rep(1, J))

fit <- lucifer(Model, MyData,
    Initial.Values = IV,
    Iterations = 1000, Status = 100,
    Thinning = 1, Algorithm = "Gibbs",
    Chains = 3)
```

### Results

The posterior inclusion probabilities for the \\\Gamma\\ indicators
reveal which cross-series dependencies are active:

``` r

post <- fit$Posterior1
## Posterior inclusion probabilities
gamma_means <- colMeans(post[, pos.Gamma])
gamma_mat <- matrix(gamma_means, J, J)
gamma_df <- expand.grid(i = 1:J, j = 1:J)
gamma_df$Posterior <- as.vector(gamma_mat)
gamma_df$True <- as.vector(true_Gamma)
gamma_table <- data.frame(
    From = gamma_df$i, To = gamma_df$j,
    `True` = gamma_df$True,
    `P(Gamma=1)` = gamma_df$Posterior,
    check.names = FALSE
)
knitr::kable(gamma_table, digits = 3, align = "rrrr",
             caption = "Posterior inclusion probabilities vs. true indicators")
```

Diagonal entries should have inclusion probabilities near 1 (strong
own-lag signal), while off-diagonal entries should be near 0 (no
cross-series coupling in the data-generating process). The continuous
parameters \\\alpha\\, \\\Phi\\, and \\\sigma\\ should recover their
true values within posterior credible intervals.

### Specialized SSVS visualization with `ssvs_summary()`

The dedicated
[`ssvs_summary()`](https://robustecologies.github.io/lucifer/reference/ssvs_summary.md)
infrastructure in lucifer wraps any spike-and-slab posterior into an S3
object that exposes seven specialized plot types covering inclusion
probabilities, the spike-vs-slab mixture structure, the conditional
posterior of each coefficient given inclusion, the matrix layout for
VAR-style models, the directed-network view, and Bayesian FDR/FOR curves
whenever ground truth is available. Using it on the auto-FC fit removes
the need for ad-hoc plotting code and provides a consistent visual
language across SSVS workflows.

``` r

gt_var <- setNames(as.vector(true_Gamma),
                   paste0("Gamma[",
                          rep(1:J, J), ",",
                          rep(1:J, each = J), "]"))

ss <- ssvs_summary(
    fit,
    indicators   = "Gamma",
    coefficients = paste0("Phi[",
                          rep(1:J, J), ",",
                          rep(1:J, each = J), "]"),
    spike_var    = 0.01,
    slab_var     = 10,
    prior_prob   = 0.5,
    ground_truth = gt_var,
    dims         = c(J, J)
)
print(ss)
```

The PIP bar chart is the canonical summary: each indicator becomes a
bar, the dotted line marks the prior probability
\\P(\Gamma\_{ij}=1)=0.5\\, and the dashed line marks the
median-probability decision threshold. Bars whose PIP is reinforced by
the data and crosses the threshold are shaded as selected, while bars
below remain grey. Ground truth status is overlaid as a separate symbol
layer so the visual sweep simultaneously reports both the posterior
decision and the truth.

``` r

plot(ss, type = "pip")
```

The heatmap reshapes the same PIP vector into the original \\J \times
J\\ adjacency layout via `dims = c(J, J)`. Tile color encodes the
inclusion probability on a diverging palette centered at the decision
boundary, and each tile is annotated with the numeric PIP and the
corresponding entry of the ground-truth indicator matrix. For the
simulated VAR(1) process the diagonal should be saturated blue and the
off-diagonal cells should remain near white.

``` r

plot(ss, type = "heatmap")
```

The network plot promotes the heatmap to a directed graph: time series
become circular nodes, edges become directed arrows whose opacity and
width are proportional to \\P(\Gamma\_{ij}=1)\\, and edges crossing the
threshold are coloured. Self-loops are intentionally absorbed into the
node markers so the graph emphasises Granger-causal cross-series
structure rather than own-lag persistence.

``` r

plot(ss, type = "network")
```

The conditional density plot decomposes the posterior of each
\\\Phi\_{ij}\\ coefficient by inclusion status, showing the distribution
conditional on \\\gamma=1\\ and the distribution conditional on
\\\gamma=0\\ as overlaid densities. This decomposition makes the
spike-and-slab mechanism visible: when the indicator is on, the
coefficient is free to take substantial values; when it is off, the
coefficient is constrained to the spike near zero. Diagonal
\\\Phi\_{jj}\\ coefficients should show essentially all of their mass in
the included branch, while off-diagonal entries should sit almost
entirely in the excluded branch.

``` r

plot(ss, type = "conditional")
```

Because
[`ssvs_summary()`](https://robustecologies.github.io/lucifer/reference/ssvs_summary.md)
was supplied with `ground_truth`, the Bayesian FDR (proportion of
selected indicators that are truly null) and Bayesian FOR (proportion of
excluded indicators that are truly active) can be computed exactly as
functions of the inclusion threshold. The resulting curves provide a
principled way to choose the threshold beyond the default 0.5: the
optimal value is whichever balances the two error rates given the cost
asymmetry of the application.

``` r

plot(ss, type = "fdr")
```

  

## Manual full-conditional derivation for VAR(p)-SSVS

While auto-FC eliminates the need for manual derivations, understanding
the full conditionals is instructive and can guide more efficient
implementations. We derive them here for the VAR(p)-SSVS model defined
above.

### Full conditional for \\\alpha_j\\

Holding all other parameters fixed, the log-posterior as a function of
\\\alpha_j\\ is quadratic, yielding a normal full conditional. Define
the residuals excluding \\\alpha_j\\:

\\r\_{t,j} = Y\_{t,j} - \sum\_{p=1}^{P} \sum\_{i=1}^{J}
\Gamma\_{i,j,p}\\ \Phi\_{i,j,p}\\ Y\_{t - L_p,\\ i}, \qquad t = L_P + 1,
\ldots, T.\\

Then \\r\_{t,j} \sim \mathcal{N}(\alpha_j, \sigma_j^2)\\ independently,
and with the prior \\\alpha_j \sim \mathcal{N}(0, 1000)\\:

\\\alpha_j \mid \text{rest} \sim \mathcal{N}(m_j, v_j),\\

where \\v_j = (T\_{\text{eff}} / \sigma_j^2 + 1/1000)^{-1}\\ and \\m_j =
v_j \cdot \sum_t r\_{t,j} / \sigma_j^2\\, with \\T\_{\text{eff}} = T -
L_P\\.

### Full conditional for \\\Gamma\_{i,k,p}\\

The full conditional of a binary indicator is a Bernoulli distribution
obtained by evaluating the log-posterior at both possible values. Let
\\\text{LP}\_g = \log p(\theta \mid y)\\ with \\\Gamma\_{i,k,p} = g\\
for \\g \in \\0, 1\\\\. Then:

\\p(\Gamma\_{i,k,p} = 1 \mid \text{rest}) = \frac{1}{1 +
\exp(\text{LP}\_0 - \text{LP}\_1)}.\\

The log-posteriors differ in two components: the prior on
\\\Phi\_{i,k,p}\\ (which switches between \\\mathcal{N}(0, \tau_0^2)\\
and \\\mathcal{N}(0, \tau_1^2)\\) and the likelihood (which includes or
excludes the \\\Phi\_{i,k,p} \cdot Y\_{t-L_p, i}\\ contribution). This
is exactly what auto-FC’s discrete enumeration computes.

### Full conditional for \\\Phi\_{i,k,p}\\

Conditional on \\\Gamma\_{i,k,p}\\:

- If \\\Gamma\_{i,k,p} = 0\\: the term \\\Gamma \cdot \Phi\\ vanishes
  from the mean function, so \\\Phi\_{i,k,p}\\ affects only the prior,
  giving \\\Phi\_{i,k,p} \mid \text{rest} \sim \mathcal{N}(0,
  \tau_0^2)\\.

- If \\\Gamma\_{i,k,p} = 1\\: define \\x_t = Y\_{t - L_p, i}\\ and let
  \\s\_{t,k}\\ be the residual of \\Y\_{t,k}\\ with all terms except
  \\\Phi\_{i,k,p} \cdot x_t\\ removed. Then \\s\_{t,k} \sim
  \mathcal{N}(\Phi\_{i,k,p} \cdot x_t, \sigma_k^2)\\, yielding a
  normal-normal conjugate update:

\\\Phi\_{i,k,p} \mid \Gamma\_{i,k,p} = 1, \text{rest} \sim
\mathcal{N}(\hat{m}, \hat{v}),\\

where \\\hat{v} = (1/\tau_1^2 + \sum_t x_t^2 / \sigma_k^2)^{-1}\\ and
\\\hat{m} = \hat{v} \cdot \sum_t s\_{t,k} \cdot x_t / \sigma_k^2\\.

### Full conditional for \\\sigma_j\\

The half-Cauchy prior on \\\sigma_j\\ is **not** conjugate with the
normal likelihood. No closed-form full conditional exists. Options
include:

- **Slice sampling** (what auto-FC uses): exact, automatic, ~5-10 Model
  evaluations.
- **Metropolis-within-Gibbs** (MWG): the approach taken in the original
  lucifer example, where `Specs = list(FC = myFC, MWG = pos.sigma)` uses
  a random-walk proposal for \\\sigma_j\\.
- **Reparameterization**: replace the half-Cauchy with an inverse-gamma
  prior to obtain conjugacy, at the cost of a different prior
  specification.

### The manual FC function

``` r

FC_varssvs <- function(parm, Data) {
    alpha <- parm[Data$pos.alpha]
    Gamma <- matrix(parm[Data$pos.Gamma], Data$J, Data$J)
    Phi   <- matrix(parm[Data$pos.Phi],   Data$J, Data$J)
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    T_eff <- Data$TT - Data$L[Data$P]
    t_idx <- (Data$L[Data$P] + 1):Data$TT

    ## --- alpha_j ---
    for (j in 1:Data$J) {
        resid_j <- Data$Y[t_idx, j]
        t_lag <- (Data$L[Data$P] - Data$L[1] + 1):(Data$TT - Data$L[1])
        resid_j <- resid_j -
            Data$Y[t_lag, ] %*% (Gamma[, j] * Phi[, j])
        prec <- T_eff / sigma[j]^2 + 1 / 1000
        m_j  <- sum(resid_j) / sigma[j]^2 / prec
        alpha[j] <- rnorm(1, m_j, 1 / sqrt(prec))
    }

    ## --- Gamma_{i,k} and Phi_{i,k} ---
    for (k in 1:Data$J) {
        for (i in 1:Data$J) {
            t_lag <- (Data$L[Data$P] - Data$L[1] + 1):(Data$TT - Data$L[1])
            x_t <- Data$Y[t_lag, i]
            ## Residual for Y_{t,k} excluding this term
            resid_ik <- Data$Y[t_idx, k] - alpha[k]
            for (ii in 1:Data$J) {
                if (ii == i) next
                resid_ik <- resid_ik -
                    Gamma[ii, k] * Phi[ii, k] * Data$Y[t_lag, ii]
            }

            ## Gamma_{i,k}: compare LP at 0 and 1
            lp1 <- log(0.5) + dnorm(Phi[i, k], 0, sqrt(10), log = TRUE) +
                sum(dnorm(resid_ik, Phi[i, k] * x_t, sigma[k], log = TRUE))
            lp0 <- log(0.5) + dnorm(Phi[i, k], 0, 0.1, log = TRUE) +
                sum(dnorm(resid_ik, 0, sigma[k], log = TRUE))
            prob1 <- 1 / (1 + exp(lp0 - lp1))
            Gamma[i, k] <- rbinom(1, 1, prob1)

            ## Phi_{i,k} | Gamma_{i,k}
            if (Gamma[i, k] == 0) {
                Phi[i, k] <- rnorm(1, 0, 0.1)
            } else {
                prec_prior <- 1 / 10
                prec_lik   <- sum(x_t^2) / sigma[k]^2
                v <- 1 / (prec_prior + prec_lik)
                m <- v * sum(resid_ik * x_t) / sigma[k]^2
                Phi[i, k] <- rnorm(1, m, sqrt(v))
            }
        }
    }

    ## sigma unchanged (sampled by MWG)
    c(alpha, as.vector(Gamma), as.vector(Phi), sigma)
}
```

To use this manual FC with MWG for sigma:

``` r

fit_manual <- lucifer(Model, MyData,
    Initial.Values = IV,
    Iterations = 10000, Status = 1000,
    Thinning = 10, Algorithm = "Gibbs",
    Specs = list(FC = FC_varssvs, MWG = pos.sigma),
    Chains = 3)
```

### Results (manual FC)

The manual-FC fit can be summarised through exactly the same
[`ssvs_summary()`](https://robustecologies.github.io/lucifer/reference/ssvs_summary.md)
machinery used for the auto-FC fit, which gives both the
inclusion-probability heatmap on the original \\J \times J\\ adjacency
layout and the directed-network view of the recovered VAR connectivity.

``` r

ss_m <- ssvs_summary(
    fit_manual,
    indicators   = "Gamma",
    coefficients = paste0("Phi[",
                          rep(1:J, J), ",",
                          rep(1:J, each = J), "]"),
    spike_var    = 0.01,
    slab_var     = 10,
    prior_prob   = 0.5,
    ground_truth = gt_var,
    dims         = c(J, J)
)
print(ss_m)
```

``` r

plot(ss_m, type = "heatmap")
```

``` r

plot(ss_m, type = "network")
```

The continuous parameters \\\alpha\\, \\\Phi\\ and \\\sigma\\ are
inspected with lucifer’s native ground-truth density grid, restricting
to the non-discrete part of the parameter vector (the binary \\\Gamma\\
indicators are already covered by the SSVS plots above).

``` r

cont_names <- parm.names[c(pos.alpha, pos.Phi, pos.sigma)]
true_cont  <- c(true_alpha, as.vector(true_Phi), true_sigma)
gt_cont    <- setNames(true_cont, cont_names)
plot(fit_manual, Parms = c("alpha", "Phi", "sigma"), ground_truth = gt_cont)
```

  

## Performance considerations

### Computational cost

The cost per Gibbs iteration scales linearly with the number of
parameters \\K\\. For a model with \\K_c\\ continuous and \\K_d\\
discrete parameters (with maximum support size \\M\\):

- Auto-FC: approximately \\(3\text{--}10) \cdot K_c + M \cdot K_d\\
  Model evaluations per iteration.
- Manual conjugate FC: \\1\\ evaluation (or \\0\\ if FC does not call
  Model).
- MWG: \\K\\ evaluations per iteration (one per parameter).

Auto-FC is therefore slower per iteration than conjugate Gibbs but
comparable to MWG, with the advantage of exact sampling (no
acceptance/rejection step) and no tuning parameters.

### Slice width

The stepping-out slice width \\w\\ affects efficiency but not
correctness. Too narrow a width requires many stepping-out expansions;
too wide causes unnecessary Model evaluations during shrinkage. The
default, derived from the Covar argument or the Laplace approximation,
is usually adequate. For parameters with dramatically different scales,
providing an informative Covar matrix or running a short adaptive
algorithm first can improve efficiency.

### When to prefer other algorithms

- **Highly correlated posteriors**: the componentwise updates of Gibbs
  (whether manual or auto-FC) mix slowly when parameters are strongly
  correlated. Consider NUTS, HMC, or adaptive multivariate algorithms
  (AM, RAM).
- **High-dimensional models** (\\K \> 100\\): the \\O(K)\\ Model
  evaluations per iteration become expensive. Gradient-based algorithms
  (HMC, NUTS, MALA) make \\O(1)\\ joint proposals.
- **Models with only continuous parameters and no conjugacy**: the Slice
  algorithm (not Gibbs) applies the same slice sampling strategy without
  the Gibbs framework overhead.

  

## References

**\[1\]** Gelfand, A.E. and Smith, A.F.M. (1990). *Sampling-based
approaches to calculating marginal densities.* Journal of the American
Statistical Association, 85(410), 398–409.
[doi:10.1080/01621459.1990.10476213](https://doi.org/10.1080/01621459.1990.10476213)

**\[2\]** Geman, S. and Geman, D. (1984). *Stochastic relaxation, Gibbs
distributions, and the Bayesian restoration of images.* IEEE
Transactions on Pattern Analysis and Machine Intelligence, 6(6),
721–741.
[doi:10.1109/TPAMI.1984.4767596](https://doi.org/10.1109/TPAMI.1984.4767596)

**\[3\]** George, E.I. and McCulloch, R.E. (1993). *Variable selection
via Gibbs sampling.* Journal of the American Statistical Association,
88(423), 881–889.
[doi:10.1080/01621459.1993.10476353](https://doi.org/10.1080/01621459.1993.10476353)

**\[4\]** Neal, R.M. (2003). *Slice sampling.* Annals of Statistics,
31(3), 705–767.
[doi:10.1214/aos/1056562461](https://doi.org/10.1214/aos/1056562461)

**\[5\]** Tanner, M.A. and Wong, W.H. (1987). *The calculation of
posterior distributions by data augmentation.* Journal of the American
Statistical Association, 82(398), 528–540.
[doi:10.1080/01621459.1987.10478458](https://doi.org/10.1080/01621459.1987.10478458)

**\[6\]** Turchin, V.F. (1971). *On the computation of multidimensional
integrals by the Monte Carlo method.* Theory of Probability and its
Applications, 16(4), 720–724.
[doi:10.1137/1116083](https://doi.org/10.1137/1116083)
