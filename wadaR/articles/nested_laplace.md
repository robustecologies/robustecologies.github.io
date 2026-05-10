# Nested Laplace approximation: INLA-like inference in lucifer

## Introduction

Integrated nested Laplace approximation (INLA) is a deterministic
strategy for approximate Bayesian inference on latent Gaussian models,
introduced by Rue, Martino, and Chopin [\[1\]](#ref1). Rather than
sampling from the posterior via MCMC, INLA decomposes the problem into a
sequence of Laplace approximations and low-dimensional numerical
integrations, achieving accuracy comparable to long MCMC runs in a
fraction of the time. The R-INLA software implements this strategy with
highly optimized sparse matrix routines and a specialized model
specification language, making it the natural tool for spatial,
temporal, and hierarchical models with Gaussian Markov random field
(GMRF) structure [\[2\]](#ref2).

Reimplementing R-INLA inside **lucifer** would be neither feasible nor
desirable: the sparse Cholesky infrastructure, GMRF theory, and
structured model specification that underpin R-INLA represent over a
decade of dedicated engineering. However, the *mathematical strategy*
behind INLA, namely the nested Laplace decomposition, is not tied to any
particular software architecture. **lucifer** already provides all the
building blocks needed to execute this strategy manually:
[`LaplaceApproximation()`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md)
for mode-finding and Hessian estimation,
[`VariationalBayes()`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md)
with Pathfinder or ADVI for fast approximate posteriors, and
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
with 55 MCMC algorithms for ground-truth comparison. This vignette
demonstrates how to compose these existing tools into an INLA-like
workflow for moderate-dimensional latent Gaussian models, where
“moderate” means problems with tens to low hundreds of parameters rather
than the thousands that require sparse matrix exploitation.

  

## The INLA decomposition

### Latent Gaussian models

INLA targets a specific class of Bayesian hierarchical models called
latent Gaussian models (LGMs). An LGM has three layers. The observations
\\y_i\\ follow some likelihood from the exponential family, conditional
on a latent field \\\boldsymbol{x}\\ and hyperparameters
\\\boldsymbol{\theta}\\:

\\ y_i \mid x_i, \boldsymbol{\theta} \sim p(y_i \mid \eta_i,
\boldsymbol{\theta}), \qquad \eta_i = g^{-1}(x_i) \\

where \\\eta_i\\ is the linear predictor and \\g\\ is a link function.
The latent field, which includes fixed effects, random effects, and any
other linear predictor components, is Gaussian conditional on the
hyperparameters:

\\ \boldsymbol{x} \mid \boldsymbol{\theta} \sim
\mathcal{N}(\boldsymbol{\mu}(\boldsymbol{\theta}),
\boldsymbol{Q}(\boldsymbol{\theta})^{-1}) \\

where \\\boldsymbol{Q}(\boldsymbol{\theta})\\ is the precision matrix.
The hyperparameters \\\boldsymbol{\theta}\\ (variance components,
correlation parameters, smoothing parameters) have their own prior
\\p(\boldsymbol{\theta})\\. The key structural assumption is that
\\\boldsymbol{\theta}\\ is low-dimensional (typically 2–6 components),
while \\\boldsymbol{x}\\ can be high-dimensional. This separation is
what makes the nested approximation efficient.

### Three-layer approximation

The posterior marginals of interest are \\\pi(x_i \mid \boldsymbol{y})\\
for each latent field component and \\\pi(\boldsymbol{\theta} \mid
\boldsymbol{y})\\ for the hyperparameters. INLA obtains these through
three nested approximations.

**Layer 1: marginal of the hyperparameters.** The marginal posterior
\\\pi(\boldsymbol{\theta} \mid \boldsymbol{y})\\ is approximated via the
Laplace formula of Tierney and Kadane [\[3\]](#ref3):

\\ \tilde{\pi}(\boldsymbol{\theta} \mid \boldsymbol{y}) \propto \left.
\frac{p(\boldsymbol{y}, \boldsymbol{x},
\boldsymbol{\theta})}{\tilde{\pi}\_G(\boldsymbol{x} \mid
\boldsymbol{\theta}, \boldsymbol{y})} \right\|\_{\boldsymbol{x} =
\boldsymbol{x}^\*(\boldsymbol{\theta})} \\

where \\\tilde{\pi}\_G(\boldsymbol{x} \mid \boldsymbol{\theta},
\boldsymbol{y})\\ is the Gaussian approximation to \\\pi(\boldsymbol{x}
\mid \boldsymbol{\theta}, \boldsymbol{y})\\ and
\\\boldsymbol{x}^\*(\boldsymbol{\theta})\\ is its mode. In practice,
this reduces to:

\\ \log \tilde{\pi}(\boldsymbol{\theta} \mid \boldsymbol{y}) = \log
p(\boldsymbol{y}, \boldsymbol{x}^\*(\boldsymbol{\theta}),
\boldsymbol{\theta}) + \frac{d_x}{2} \log(2\pi) + \frac{1}{2} \log
\|\boldsymbol{\Sigma}\_x(\boldsymbol{\theta})\| + \text{const} \\

where \\d_x = \dim(\boldsymbol{x})\\ and
\\\boldsymbol{\Sigma}\_x(\boldsymbol{\theta})\\ is the covariance matrix
of the Gaussian approximation (the inverse of the negative Hessian of
\\\log p(\boldsymbol{y}, \boldsymbol{x}, \boldsymbol{\theta})\\ with
respect to \\\boldsymbol{x}\\, evaluated at the conditional mode).

**Layer 2: conditional marginals of the latent field.** For each
hyperparameter configuration \\\boldsymbol{\theta}\_k\\, the conditional
posterior of the latent field component \\x_i\\ is approximated as
Gaussian:

\\ \tilde{\pi}(x_i \mid \boldsymbol{\theta}\_k, \boldsymbol{y}) =
\mathcal{N}\\\left(x_i;\\ \mu_i(\boldsymbol{\theta}\_k),\\
\sigma_i^2(\boldsymbol{\theta}\_k)\right) \\

where \\\mu_i(\boldsymbol{\theta}\_k)\\ is the \\i\\-th component of the
conditional mode \\\boldsymbol{x}^\*(\boldsymbol{\theta}\_k)\\ and
\\\sigma_i^2(\boldsymbol{\theta}\_k)\\ is the \\i\\-th diagonal element
of \\\boldsymbol{\Sigma}\_x(\boldsymbol{\theta}\_k)\\. R-INLA offers
more refined corrections (simplified Laplace, full Laplace) that improve
accuracy for skewed conditionals, but the Gaussian approximation is the
fastest and often sufficient.

**Layer 3: numerical integration.** The unconditional marginals are
obtained by integrating out \\\boldsymbol{\theta}\\ numerically over a
grid \\\\\boldsymbol{\theta}\_1, \ldots, \boldsymbol{\theta}\_K\\\\:

\\ \tilde{\pi}(x_i \mid \boldsymbol{y}) = \sum\_{k=1}^{K}
\tilde{\pi}(x_i \mid \boldsymbol{\theta}\_k, \boldsymbol{y}) \times
\tilde{\pi}(\boldsymbol{\theta}\_k \mid \boldsymbol{y}) \times \Delta_k
\\

The result is a Gaussian mixture, which can represent asymmetry, heavy
tails, and other non-Gaussian features that a single Gaussian
approximation would miss.

### Mapping to lucifer’s toolkit

Each layer maps directly to an existing **lucifer** function:

| INLA layer | Tool | What it provides |
|:---|:---|:---|
| Layer 1: explore \\\boldsymbol{\theta}\\-space | [`LaplaceApproximation()`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md) | Joint mode and marginal variance of \\\boldsymbol{\theta}\\ for grid construction |
| Layer 1: find \\\boldsymbol{x}^\*(\boldsymbol{\theta})\\ and \\\boldsymbol{\Sigma}\_x(\boldsymbol{\theta})\\ | [`optim()`](https://rdrr.io/r/stats/optim.html) + [`Hessian()`](https://robustecologies.github.io/lucifer/reference/Matrices.md) | Conditional mode and accurate Hessian-based covariance |
| Layer 3: numerical integration | Manual grid summation | Weighted Gaussian mixture |
| Validation | [`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md) with NUTS/HMCDA | Ground-truth MCMC posterior |
| Fast alternative | [`VariationalBayes()`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md) with Pathfinder or ADVI | Gaussian approximation without nesting |

The workflow requires the standard `Model` function (all parameters in
`parm`) for the joint Laplace approximation and MCMC, plus a plain R
function for the conditional log-posterior (latent field as argument,
hyperparameters as fixed parameters) used by
[`optim()`](https://rdrr.io/r/stats/optim.html) in the nested loop. The
Hessian at each conditional mode is computed via **lucifer**’s
[`Hessian()`](https://robustecologies.github.io/lucifer/reference/Matrices.md),
which uses Richardson extrapolation with step sizes calibrated for
second derivatives; this accuracy is essential for the Laplace marginal
formula, where the log-determinant of the conditional covariance drives
the hyperparameter weights.

  

## Example: Poisson GLMM with random intercepts

### Data simulation

Consider a Poisson generalized linear mixed model with group-level
random intercepts, a standard example in the INLA literature. Counts
\\y\_{ij}\\ from group \\j\\ follow a Poisson distribution with
log-linear predictor:

\\ y\_{ij} \sim \text{Poisson}(\exp(\eta\_{ij})), \qquad \eta\_{ij} =
\beta_0 + \beta_1 \\ x\_{ij} + u_j \\

where \\u_j \sim \mathcal{N}(0, \sigma_u^2)\\ are group-level random
intercepts. The hyperparameter is \\\theta = \log \sigma_u\\, which
controls the between-group variability. The latent field comprises
\\\boldsymbol{x} = (\beta_0, \beta_1, u_1, \ldots, u_J)\\, giving \\d_x
= J + 2\\ parameters to optimize at each grid point.

``` r

library(lucifer)
set.seed(42)

# Design
J   <- 10       # groups
n_j <- 20       # observations per group
N   <- J * n_j  # total

# True parameter values
beta0_true   <- 1.0
beta1_true   <- 0.5
sigma_u_true <- 0.8

# Simulate
group <- rep(1:J, each = n_j)
x_cov <- rnorm(N, 0, 1)
u_true <- rnorm(J, 0, sigma_u_true)
eta <- beta0_true + beta1_true * x_cov + u_true[group]
y <- rpois(N, exp(eta))
```

### Model specification

The joint model places weakly informative priors on the fixed effects
and a half-Cauchy prior on \\\sigma_u\\. Since **lucifer** operates on a
continuous, unconstrained parameter space, we parameterize the variance
component as \\\theta = \log \sigma_u\\ and include the Jacobian
correction \\\|\mathrm{d}\sigma_u / \mathrm{d}\theta\| = \sigma_u\\ in
the log-posterior.

``` r

parm.names <- c("beta[1]", "beta[2]", "log.sigma.u",
                paste0("u[", 1:J, "]"))
mon.names  <- c("sigma.u", "LP")

Data <- list(
  N = N, J = J, y = y, x = x_cov, group = group,
  parm.names = parm.names, mon.names = mon.names
)

Model <- function(parm, Data) {
  # --- extract parameters ---
  beta  <- parm[1:2]
  log_sigma_u <- parm[3]
  sigma_u <- exp(log_sigma_u)
  u <- parm[4:(3 + Data$J)]

  # --- priors ---
  beta_prior    <- sum(dnorm(beta, 0, 10, log = TRUE))
  sigma_u_prior <- dhalfcauchy(sigma_u, 5, log = TRUE) + log_sigma_u
  u_prior       <- sum(dnorm(u, 0, sigma_u, log = TRUE))

  # --- likelihood ---
  eta <- beta[1] + beta[2] * Data$x + u[Data$group]
  mu  <- exp(eta)
  LL  <- sum(dpois(Data$y, mu, log = TRUE))

  # --- assemble ---
  LP  <- LL + beta_prior + sigma_u_prior + u_prior
  Dev <- -2 * LL
  yhat <- rpois(Data$N, mu)

  return(list(LP = LP, Dev = Dev,
              Monitor = c(sigma_u, LP),
              yhat = yhat, parm = parm))
}

# Initial values
parm0 <- c(0, 0, log(1), rep(0, J))
```

  

## Approach 1: joint Laplace approximation

The simplest approximation runs
[`LaplaceApproximation()`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md)
on the full model, treating all parameters (fixed effects, log-variance,
random effects) jointly. This finds the joint posterior mode and
estimates a Gaussian posterior via the Hessian. It corresponds to INLA’s
Gaussian approximation *without* nesting: a single multivariate Gaussian
centered at the joint mode.

``` r

fit_joint <- LaplaceApproximation(
  Model, parm0, Data,
  Method     = "TR",
  Iterations = 300,
  sir        = FALSE
)

# Joint mode and covariance
joint_mode  <- fit_joint$Summary1[, "Mode"]
joint_covar <- fit_joint$Covar
joint_sd    <- fit_joint$Summary1[, "SD"]

# Marginal for each parameter: N(mode_i, sd_i^2)
fit_joint$Summary1
```

This gives reasonable point estimates but assumes the joint posterior is
Gaussian, which is often inadequate for variance components. The
marginal of \\\log \sigma_u\\ tends to be skewed, especially with few
groups, and a single Gaussian cannot capture this asymmetry.

  

## Approach 2: nested Laplace approximation

The INLA-like approach exploits the fact that \\\theta = \log \sigma_u\\
is one-dimensional while the latent field \\\boldsymbol{x}\\ has \\d_x =
12\\ components. We evaluate the conditional Laplace approximation at a
grid of \\\theta\\ values, weight each configuration by
\\\tilde{\pi}(\theta \mid \boldsymbol{y})\\, and combine the conditional
Gaussian marginals into a mixture.

### The conditional log-posterior

We define the conditional log-posterior as a plain R function of the
latent field \\\boldsymbol{x}\\, with \\\theta\\ passed as a fixed
argument. This function will be maximized with
[`optim()`](https://rdrr.io/r/stats/optim.html) at each grid point, and
its Hessian computed with **lucifer**’s
[`Hessian()`](https://robustecologies.github.io/lucifer/reference/Matrices.md)
for accurate covariance estimation. The hyperparameter prior is included
so that the log-posterior at the mode gives the correct Laplace marginal
weight.

``` r

neg_LP_cond <- function(parm_x, log_sigma_u, Data) {
  beta <- parm_x[1:2]
  u    <- parm_x[3:(2 + Data$J)]
  sigma_u <- exp(log_sigma_u)

  beta_prior    <- sum(dnorm(beta, 0, 10, log = TRUE))
  sigma_u_prior <- dhalfcauchy(sigma_u, 5, log = TRUE) + log_sigma_u
  u_prior       <- sum(dnorm(u, 0, sigma_u, log = TRUE))

  eta <- beta[1] + beta[2] * Data$x + u[Data$group]
  mu  <- exp(eta)
  LL  <- sum(dpois(Data$y, mu, log = TRUE))

  -(LL + beta_prior + sigma_u_prior + u_prior)
}
```

### Grid construction

The grid over \\\theta\\ is centered at the joint mode with spread
determined by the joint Hessian. We use \\\pm 4\\ standard deviations to
ensure the grid captures the full marginal, including any rightward skew
that is typical of variance components with moderate group counts.

``` r

theta_mode <- joint_mode["log.sigma.u"]
theta_sd   <- joint_sd["log.sigma.u"]

K <- 31
theta_grid <- seq(theta_mode - 4 * theta_sd,
                  theta_mode + 4 * theta_sd,
                  length.out = K)
```

### Conditional optimization loop

For each \\\theta_k\\, we maximize the conditional log-posterior with
[`optim()`](https://rdrr.io/r/stats/optim.html) (BFGS) and compute the
Hessian at the mode with
[`Hessian()`](https://robustecologies.github.io/lucifer/reference/Matrices.md).
To interface with
[`Hessian()`](https://robustecologies.github.io/lucifer/reference/Matrices.md),
we wrap the conditional log-posterior in a thin `Model`-compatible
function that returns the standard `list(LP, Dev, Monitor, yhat, parm)`.
The Hessian is inverted to obtain the conditional covariance matrix
\\\boldsymbol{\Sigma}\_x(\theta_k)\\. Warm-starting from the previous
grid point’s mode ensures fast convergence across the grid.

``` r

d_x <- 2 + J  # dimension of latent field

# Initialize from joint mode (excluding log.sigma.u)
parm_x_init <- joint_mode[c(1:2, 4:(3 + J))]

# Wrapper for Hessian(): given a fixed theta_k, returns a Model-compatible
# function whose LP is the conditional log-posterior of the latent field
make_cond_model <- function(log_sigma_u, Data) {
  function(parm, Data) {
    LP <- -neg_LP_cond(parm, log_sigma_u = log_sigma_u, Data = Data)
    list(LP = LP, Dev = -2 * LP, Monitor = 1, yhat = 1, parm = parm)
  }
}

results <- vector("list", K)

for (k in seq_len(K)) {
  theta_k <- theta_grid[k]

  # Find conditional mode
  opt <- optim(parm_x_init, neg_LP_cond,
               log_sigma_u = theta_k, Data = Data,
               method = "BFGS", control = list(maxit = 500))

  # Hessian at mode via lucifer's Hessian() (Richardson, d = 0.1)
  cond_model <- make_cond_model(theta_k, Data)
  H <- -Hessian(cond_model, opt$par, Data)
  Sigma <- solve(H)
  log_det <- as.numeric(determinant(Sigma, logarithm = TRUE)$modulus)

  results[[k]] <- list(
    mode    = opt$par,
    covar   = Sigma,
    LP      = -opt$value,
    log_det = log_det
  )

  # Warm-start next grid point
  parm_x_init <- opt$par
}
```

### Marginal posteriors via Gaussian mixture

The unnormalized log-weight for each grid point is the Laplace
approximation to \\\log \pi(\theta_k \mid \boldsymbol{y})\\, computed
from the conditional log-posterior at the mode and the log-determinant
of the conditional covariance:

``` r

log_w <- sapply(results, function(r) {
  r$LP + 0.5 * d_x * log(2 * pi) + 0.5 * r$log_det
})

# Normalize (log-sum-exp for stability)
log_w <- log_w - max(log_w)
w <- exp(log_w)
w <- w / sum(w)
```

The marginal posterior of any latent field component \\x_i\\ is now a
Gaussian mixture:

\\ \tilde{\pi}(x_i \mid \boldsymbol{y}) = \sum\_{k=1}^{K} w_k \\
\mathcal{N}\\\left(x_i;\\ \mu_i^{(k)},\\ \sigma_i^{2(k)}\right) \\

where \\\mu_i^{(k)}\\ and \\\sigma_i^{2(k)}\\ are the conditional mode
and variance from the \\k\\-th grid point. We define a helper function
to evaluate this mixture:

``` r

gmix_density <- function(x_vals, param_idx, results, weights) {
  dens <- numeric(length(x_vals))
  for (k in seq_along(results)) {
    mu_k <- results[[k]]$mode[param_idx]
    sd_k <- sqrt(results[[k]]$covar[param_idx, param_idx])
    dens <- dens + weights[k] * dnorm(x_vals, mu_k, sd_k)
  }
  dens
}
```

The marginal of \\\theta\\ itself (and hence \\\sigma_u = e^\theta\\) is
given directly by the normalized weights on the grid. We can obtain a
smooth density estimate via interpolation:

``` r

# Marginal density of log(sigma_u) on the grid
theta_dens <- w / diff(theta_grid)[1]  # rescale to density units

# Smooth interpolation
theta_smooth <- splinefun(theta_grid, theta_dens)
```

  

## Ground truth: MCMC

Full MCMC with the No-U-Turn Sampler provides the reference posterior
against which all approximations are compared.

``` r

fit_mcmc <- lucifer(
  Model, Data, parm0,
  Iterations = 20000,
  Algorithm  = "NUTS",
  Specs      = NULL,
  Status     = 5000,
  Thinning   = 5
)

# Stationary posterior samples
posterior <- fit_mcmc$Posterior2
```

  

## Comparing marginal posteriors

The payoff of the nested Laplace approach is visible when we overlay
marginal density estimates from all methods. The comparison highlights
three regimes: parameters where all methods agree (well-identified fixed
effects), parameters where nesting improves on joint Laplace (variance
components with skewed marginals), and random effects whose marginal
shape depends on the hyperparameter uncertainty.

``` r

library(ggplot2)

# Parameters to compare
params <- list(
  list(name = expression(beta[0]),      joint_idx = 1, cond_idx = 1, mcmc_col = 1),
  list(name = expression(beta[1]),      joint_idx = 2, cond_idx = 2, mcmc_col = 2),
  list(name = expression(log~sigma[u]), joint_idx = 3, cond_idx = NA, mcmc_col = 3),
  list(name = expression(u[1]),         joint_idx = 4, cond_idx = 3, mcmc_col = 4)
)

plot_list <- list()

for (i in seq_along(params)) {
  p <- params[[i]]
  mcmc_samp <- posterior[, p$mcmc_col]
  xlim <- range(mcmc_samp) + c(-0.3, 0.3) * diff(range(mcmc_samp))
  x_seq <- seq(xlim[1], xlim[2], length.out = 300)

  # MCMC kernel density
  kd <- density(mcmc_samp, from = xlim[1], to = xlim[2], n = 300)
  d_mcmc <- approx(kd$x, kd$y, xout = x_seq)$y

  # Joint LA: single Gaussian
  d_joint <- dnorm(x_seq, joint_mode[p$joint_idx], joint_sd[p$joint_idx])

  # Nested LA: Gaussian mixture (or theta marginal)
  if (!is.na(p$cond_idx)) {
    d_nested <- gmix_density(x_seq, p$cond_idx, results, w)
  } else {
    d_nested <- pmax(0, theta_smooth(x_seq))
  }

  df <- data.frame(
    value   = rep(x_seq, 3),
    density = c(d_mcmc, d_joint, d_nested),
    method  = factor(
      rep(c("MCMC (NUTS)", "Joint Laplace", "Nested Laplace"), each = 300),
      levels = c("MCMC (NUTS)", "Nested Laplace", "Joint Laplace")
    )
  )

  plot_list[[i]] <- ggplot(df, aes(x = value, y = density,
                                    colour = method, linetype = method)) +
    geom_line(linewidth = 0.8) +
    labs(x = p$name, y = "Density") +
    scale_colour_manual(values = c("black", "#E41A1C", "#377EB8")) +
    scale_linetype_manual(values = c("solid", "solid", "dashed")) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none")
}

# Add legend to last panel and arrange
plot_list[[4]] <- plot_list[[4]] +
  theme(legend.position = "bottom", legend.title = element_blank())

library(gridExtra)
gridExtra::grid.arrange(
  plot_list[[1]], plot_list[[2]],
  plot_list[[3]], plot_list[[4]],
  nrow = 2
)
```

The comparison reveals the expected pattern. For \\\beta_1\\, all three
methods produce nearly identical marginals because the slope is
well-identified by the data and its posterior is close to Gaussian. For
\\\beta_0\\ and the random effect \\u_1\\, the nested Laplace marginal
is wider than the joint Laplace marginal, because it integrates over
uncertainty in \\\sigma_u\\ that the joint approximation absorbs into a
single mode. The most informative comparison appears in \\\log
\sigma_u\\: the MCMC marginal (black) is right-skewed, typical for
variance components with moderate group counts; the nested Laplace
mixture (red) captures this skewness through the varying weights across
the \\\theta\\-grid; and the joint Laplace marginal (blue dashed)
imposes symmetry by construction, underestimating both the mean and the
variance.

  

## Quantifying the approximation

Beyond visual comparison, we can assess the quality of each
approximation numerically. Two useful diagnostics are the marginal mean
and standard deviation (which should match MCMC), and the coverage of
credible intervals.

``` r

# Nested Laplace: mixture mean and variance for parameter i
gmix_moments <- function(param_idx, results, weights) {
  mu <- sapply(results, function(r) r$mode[param_idx])
  v  <- sapply(results, function(r) r$covar[param_idx, param_idx])
  mean_mix <- sum(weights * mu)
  var_mix  <- sum(weights * (v + mu^2)) - mean_mix^2
  c(mean = mean_mix, sd = sqrt(var_mix))
}

# Compare moments for beta[1], beta[2], log.sigma.u, u[1]
idx <- c(1, 2, 3, 4)
rows <- lapply(idx, function(i) {
  mcmc_m <- mean(posterior[, i])
  mcmc_s <- sd(posterior[, i])
  jla_m  <- joint_mode[i]
  jla_s  <- joint_sd[i]
  ci     <- if (i <= 2) i else if (i == 3) NA else i - 1
  if (!is.na(ci)) {
    nla <- gmix_moments(ci, results, w)
  } else {
    nla <- c(mean = sum(w * theta_grid),
             sd = sqrt(sum(w * theta_grid^2) - sum(w * theta_grid)^2))
  }
  data.frame(
    Parameter = Data$parm.names[i],
    MCMC      = sprintf("%.3f (%.3f)", mcmc_m, mcmc_s),
    Joint.LA  = sprintf("%.3f (%.3f)", jla_m, jla_s),
    Nested.LA = sprintf("%.3f (%.3f)", nla[1], nla[2])
  )
})
diag_df <- do.call(rbind, rows)
knitr::kable(diag_df, align = "lccc",
             col.names = c("Parameter", "MCMC (NUTS)", "Joint Laplace",
                           "Nested Laplace"),
             caption = "Posterior mean (sd) by approximation method")
```

  

## Generalizing to higher-dimensional hyperparameters

The one-dimensional grid above is computationally trivial, but the
nested Laplace strategy extends naturally to low-dimensional
\\\boldsymbol{\theta}\\. For \\\dim(\boldsymbol{\theta}) = 2\\, the grid
becomes a two-dimensional lattice; for \\\dim(\boldsymbol{\theta}) =
3\\, a three-dimensional lattice. The cost grows as \\O(K^{d\_\theta}
\times C\_{\text{LA}})\\ where \\C\_{\text{LA}}\\ is the cost of each
conditional Laplace approximation, making the approach practical for
\\d\_\theta \leq 4\\ or so.

For higher-dimensional hyperparameter spaces, R-INLA uses a central
composite design (CCD) or similar space-filling grid to keep the number
of evaluation points manageable. In **lucifer**, one could combine the
nested Laplace loop with
[`IterativeQuadrature()`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md)
using the CAGH (componentwise adaptive Gauss-Hermite) algorithm, which
adaptively places quadrature points along each hyperparameter dimension.
The conditional `Model_cond` function would be called at each quadrature
node, and the integration handled by the IQ machinery rather than a
manual grid.

  

## Discussion

### When does nesting help?

The nested Laplace approximation provides the most benefit when two
conditions hold simultaneously. First, the marginal posterior of the
hyperparameters is non-Gaussian: skewed, heavy-tailed, or multimodal.
This happens most often with variance components when the number of
groups is small (say, \\J \< 20\\), making the likelihood information
about \\\sigma_u\\ limited relative to the prior. Second, the latent
field marginals depend meaningfully on the hyperparameter value, so that
integrating over \\\boldsymbol{\theta}\\ uncertainty changes their
shape. When the posterior is close to Gaussian everywhere (large data,
well-identified parameters), the joint Laplace approximation suffices
and nesting adds computational cost without improving accuracy.

### Limitations relative to R-INLA

Three structural limitations separate this manual workflow from R-INLA.
First, **lucifer** uses dense matrix algebra, so the cost of each
conditional Laplace approximation scales as \\O(d_x^3)\\ in the Hessian
computation. R-INLA exploits the sparse precision matrix
\\\boldsymbol{Q}(\boldsymbol{\theta})\\ inherent in GMRFs, reducing this
to \\O(d_x)\\ or \\O(d_x^{3/2})\\ depending on the sparsity pattern. For
models with thousands of latent parameters (spatial fields, spline
bases), this difference is decisive. Second, R-INLA offers the
simplified Laplace correction [\[1\]](#ref1) for the conditional
marginals, which improves accuracy for skewed posteriors without the
full cost of a per-component Laplace approximation. The approach
presented here uses only the Gaussian approximation at Layer 2. Third,
R-INLA includes the SPDE approach for continuous spatial models
[\[4\]](#ref4), which approximates Gaussian random fields by GMRFs on
triangulated meshes; this has no analogue in **lucifer**.

### Practical recommendations in lucifer

For models where the nested approach is appropriate (moderate dimension,
few hyperparameters, non-Gaussian hyperparameter marginals), the
workflow presented here provides a practical middle ground between joint
Laplace approximation (fast but potentially inaccurate) and full MCMC
(accurate but slow). The following guidelines help decide which tool to
use:

- **Joint Laplace**
  ([`LaplaceApproximation()`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md)):
  start here. If the mode is well-defined and the Hessian is positive
  definite, the Gaussian approximation may be sufficient. Check by
  comparing against a short MCMC run.
- **Nested Laplace** (manual loop): use when the joint Laplace marginals
  for variance components are clearly wrong (e.g., substantial
  probability mass below zero for \\\log \sigma_u\\ indicates the
  Gaussian is too wide, or the marginal appears symmetric when it should
  be skewed).
- **Pathfinder** (`VariationalBayes(Method = "Pathfinder")`): a fast
  alternative that explores multiple optimization trajectories. It
  provides Gaussian approximations like joint Laplace but can sometimes
  find better modes.
- **Full MCMC**
  ([`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
  with NUTS): the gold standard. Use when approximation quality is
  uncertain or when the posterior is expected to be strongly
  non-Gaussian.

  

## References

**\[1\]** Rue, H., Martino, S. & Chopin, N. (2009). Approximate Bayesian
inference for latent Gaussian models by using integrated nested Laplace
approximations. *Journal of the Royal Statistical Society: Series B*,
**71**(2), 319–392.
[doi:10.1111/j.1467-9868.2008.00700.x](https://doi.org/10.1111/j.1467-9868.2008.00700.x)

**\[2\]** Rue, H. & Held, L. (2005). *Gaussian Markov random fields:
theory and applications*. Chapman & Hall/CRC. ISBN: 978-1584884323.

**\[3\]** Tierney, L. & Kadane, J.B. (1986). Accurate approximations for
posterior moments and marginal densities. *Journal of the American
Statistical Association*, **81**(393), 82–86.
[doi:10.1080/01621459.1986.10478240](https://doi.org/10.1080/01621459.1986.10478240)

**\[4\]** Lindgren, F., Rue, H. & Lindstrom, J. (2011). An explicit link
between Gaussian fields and Gaussian Markov random fields: the
stochastic partial differential equation approach. *Journal of the Royal
Statistical Society: Series B*, **73**(4), 423–498.
[doi:10.1111/j.1467-9868.2011.00777.x](https://doi.org/10.1111/j.1467-9868.2011.00777.x)
