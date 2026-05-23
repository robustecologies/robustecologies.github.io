# Stochastic differential equation model

Constructs a structured model specification for stochastic dynamical
systems that automatically generates a
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)-compatible
Model function. The user specifies drift, diffusion, observation model,
and priors; the package handles discretization, likelihood computation,
and (optionally) particle filtering. The generated Model function plugs
directly into the existing MCMC, SMC, and ABC infrastructure with zero
modifications.

## Usage

``` r
SDE(
  drift = NULL,
  diffusion = NULL,
  data,
  times,
  x0 = NULL,
  obs.model = "gaussian",
  obs.noise = NULL,
  obs.loglik = NULL,
  obs.link = NULL,
  obs.params = NULL,
  prior = NULL,
  parm.names = NULL,
  state.names = NULL,
  method = "ekf",
  family = NULL,
  dt = NULL,
  N.particles = 200L,
  pf.type = "bootstrap",
  ess.threshold = 0.5,
  constraints = NULL,
  Initial.Values = NULL,
  jump.rate = NULL,
  jump.size = NULL,
  compile = FALSE,
  vectorized = FALSE,
  interpretation = "ito"
)
```

## Arguments

- drift:

  A function with signature `function(x, theta, t)` that returns the
  drift vector (length `dim_x`).

- diffusion:

  A function with signature `function(x, theta, t)` that returns the
  diffusion coefficient. May return a scalar (isotropic noise), a vector
  of length `dim_x` (diagonal noise), or a `dim_x x dim_x` matrix
  (correlated noise).

- data:

  Numeric vector or matrix of observations. For multivariate
  observations, rows are time points and columns are variables.

- times:

  Numeric vector of observation times with length equal to the number of
  observations.

- x0:

  Initial state specification. Either a numeric vector of fixed initial
  values, or a function with signature `function(theta)` that returns
  the initial state vector. Functions are useful when the initial state
  depends on parameters.

- obs.model:

  Character string specifying the observation model: `"gaussian"`
  (default), `"poisson"`, `"binomial"`, `"custom"`, or any family
  registered in the SSM observation registry (`"negbin"`, `"student_t"`,
  `"zero_inflated_poisson"`). When a registry family is selected and
  `obs.params` is provided, the observation log-likelihood is built
  automatically from the registry without requiring a custom
  `obs.loglik` function.

- obs.noise:

  A function with signature `function(theta)` that returns the
  observation noise specification. Required when
  `obs.model = "gaussian"`. Three return types are supported: a scalar
  (isotropic noise, standard deviation applied to all observation
  dimensions); a numeric vector of length `dim_y` (diagonal noise, one
  standard deviation per observation dimension); or a `dim_y x dim_y`
  matrix (full observation noise covariance matrix, used directly).

- obs.loglik:

  A function with signature `function(y, x, theta)` that returns the
  log-likelihood of observation `y` given latent state `x` and
  parameters `theta`. Required when `obs.model = "custom"`.

- obs.link:

  A function with signature `function(x, theta)` that maps the latent
  state to the observation mean. Required for `"poisson"` and
  `"binomial"` observation models; optional for `"gaussian"` (defaults
  to identity).

- obs.params:

  Optional named list of observation-model parameters for registry-based
  families. Each element is either a fixed scalar/vector or a function
  of `theta` for parameters estimated jointly with the latent states.
  For example, `list(sigma = function(theta) exp(theta[4]))` makes the
  Gaussian observation noise a function of the fourth parameter.
  Required scale parameters depend on the chosen family: `"gaussian"`
  needs `sigma`; `"negbin"` needs `size`; `"student_t"` needs `sigma`
  and `df`; `"zero_inflated_poisson"` needs `p_zero`; `"binomial"` needs
  `n_trials`. When `obs.params` is provided with a registered family,
  the observation log-likelihood is constructed automatically and
  `obs.loglik`, `obs.noise`, and `obs.link` are not required (though a
  custom `obs.link` overrides the registry default if supplied).

- prior:

  A function with signature `function(theta)` that returns the log-prior
  density for the parameter vector.

- parm.names:

  Character vector of parameter names.

- state.names:

  Optional character vector of state variable names. Defaults to `x1`,
  `x2`, etc.

- method:

  Character string specifying the discretization method: `"exact"`
  (analytical transition densities for known families), `"ekf"`
  (Extended Kalman Filter, deterministic approximate likelihood,
  recommended for Gaussian observations), `"euler"` (stochastic
  Euler-Maruyama single-trajectory approximation), `"milstein"`
  (stochastic Milstein scheme), or `"particle"` (bootstrap particle
  filter).

- family:

  Character string identifying the SDE family for exact methods: `"ou"`,
  `"gbm"`, `"cir"`, or `"vasicek"`. Required when `method = "exact"`.

- dt:

  Numeric discretization step for Euler-Maruyama and Milstein schemes.
  Defaults to `min(diff(times)) / 10`, clamped to
  `[1e-6, min(diff(times))]`.

- N.particles:

  Integer number of particles for the bootstrap particle filter. Default
  200.

- pf.type:

  Character string specifying the particle filter type: `"bootstrap"`
  (default), `"auxiliary"` (look-ahead weighting, better ESS for
  informative observations), or `"bridge"` (guided proposals toward next
  observation, best for sparse data). Only used when
  `method = "particle"`.

- ess.threshold:

  Numeric fraction in `(0, 1)` controlling the ESS-based adaptive
  resampling threshold for the particle filter. The filter resamples
  when ESS \< `ess.threshold * N.particles`. Default 0.5.

- constraints:

  Optional parameter constraints. Either a named list where each element
  is a numeric vector `c(lower, upper)` keyed by parameter name, or a
  matrix with `dim_theta` rows and 2 columns (lower bound, upper bound).
  Parameters not listed in a named list are left unconstrained. The
  generated Model function applies
  [`interval`](https://robustecologies.github.io/lucifer/reference/interval.md)
  to enforce these bounds at every evaluation, preventing the sampler
  from exploring degenerate regions (e.g., negative scale parameters).
  When `Initial.Values` is `NULL` and constraints are provided, sensible
  starting values are generated automatically: 1.0 for `(0, Inf)`,
  midpoint for finite intervals, and 0 for unconstrained parameters.

- Initial.Values:

  Optional numeric vector of starting parameter values. If `NULL` and
  `constraints` are provided, smart defaults are generated from the
  constraint bounds. If both are `NULL`, defaults to
  `rep(0.1, dim_theta)`.

- jump.rate:

  Optional function with signature `function(x, theta, t)` returning the
  Poisson jump rate for jump-diffusion models.

- jump.size:

  Optional function with signature `function(x, theta, t)` returning a
  random jump size vector for jump-diffusion models.

- compile:

  Logical. If `TRUE` and `drift`/`diffusion` were supplied as character
  strings, expressions, or formulas, the drift/diffusion expressions are
  compiled to native C++ via
  [`compile.sde_model`](https://robustecologies.github.io/lucifer/reference/compile.sde_model.md)
  for zero-overhead dispatch. Default `FALSE`.

- vectorized:

  Logical. If `TRUE`, the user-supplied `drift` and `diffusion` closures
  accept and return matrices (one row per particle) so that
  particle-filter and EKF engines can evaluate all particles in a single
  call. Default `FALSE` (scalar per-particle calls).

- interpretation:

  Character; either `"ito"` (default) or `"stratonovich"`. Selects the
  stochastic-calculus interpretation used to derive transition densities
  and discretisation schemes.

## Value

An object of class `sde_model`, a list containing:

- Model:

  A function with signature `function(parm, Data)` returning the
  standard five-component list (LP, Dev, Monitor, yhat, parm) for use
  with
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

- Data:

  A list containing `parm.names`, `mon.names`, and internal SDE
  specification fields.

- Initial.Values:

  Numeric vector of starting parameter values.

- drift:

  The drift function.

- diffusion:

  The diffusion function.

- dim_x:

  State dimension.

- dim_theta:

  Parameter dimension.

- method:

  Discretization method used.

- family:

  SDE family (if exact).

- obs.model:

  Observation model type.

- state.names:

  State variable names.

- parm.names:

  Parameter names.

- times:

  Observation times.

- data:

  Observed data.

## Details

Stochastic differential equation model

The core design is a bridge pattern: `SDE()` constructs an S3 object of
class `sde_model` that encapsulates the SDE specification and produces a
standard `Model(parm, Data)` closure. This closure internally dispatches
to C++ likelihood engines depending on the model structure. The user
then passes `sde$Model`, `sde$Data`, and `sde$Initial.Values` to
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`SMC`](https://robustecologies.github.io/lucifer/reference/SMC.md), or
[`ABC`](https://robustecologies.github.io/lucifer/reference/ABC.md)
unchanged.

Four likelihood engines are available. The `"exact"` method uses
closed-form Gaussian (OU, GBM, Vasicek) or non-central chi-squared (CIR)
transition densities for maximum speed and accuracy. The `"euler"`
method approximates the transition density via Euler-Maruyama
discretization with sub-stepping, treating each transition as Gaussian
with mean \\x + f(x)\Delta t\\ and variance \\g(x)^2 \Delta t\\. The
`"milstein"` method adds the higher-order correction \\0.5 g g' (\Delta
W^2 - \Delta t)\\ to reduce discretization bias. The `"particle"` method
uses a bootstrap particle filter to produce an unbiased
log-marginal-likelihood estimate, suitable for pseudo-marginal MCMC via
non-Gaussian or highly nonlinear observation models.

Jump-diffusion processes of the form \\dX = f(X)\\dt + g(X)\\dW +
J\\dN(\lambda)\\ are supported by specifying `jump.rate` and
`jump.size`. Jumps are integrated into the Euler-Maruyama, Milstein, and
particle filter engines.

Parameter constraints should be specified for any parameter that must
remain within a specific range (e.g., positive scale parameters). The
generated Model function calls
[`interval`](https://robustecologies.github.io/lucifer/reference/interval.md)
on each constrained parameter before evaluating the likelihood, and
returns the constrained values in `parm` so the sampler tracks the
feasible region. Without constraints, scale parameters initialized near
zero can produce degenerate likelihoods that trap the sampler.

**Algorithm choice.** Euler, Milstein, and particle-filter likelihoods
call user-supplied R functions (drift, diffusion) at every sub-step,
making each evaluation expensive. Gradient-based algorithms (NUTS,
HMCDA, HMC, MALA) multiply this cost by \\2d+1\\ per gradient evaluation
and by the number of leapfrog steps per iteration. For most SDE models,
non-gradient samplers such as HARM, AMWG, twalk, or RAM produce
comparable effective sample sizes at a fraction of the wall-clock time.
Reserve NUTS for models with user-supplied analytical gradients via the
`Grad` argument to
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
or for low-dimensional models where the per-iteration cost is
acceptable.

## References

Iacus, S.M. (2008). *Simulation and Inference for Stochastic
Differential Equations*. Springer.
[doi:10.1007/978-0-387-75839-8](https://doi.org/10.1007/978-0-387-75839-8)

Kloeden, P.E. and Platen, E. (1992). *Numerical Solution of Stochastic
Differential Equations*. Springer.
[doi:10.1007/978-3-662-12616-5](https://doi.org/10.1007/978-3-662-12616-5)

Andrieu, C., Doucet, A., and Holenstein, R. (2010). "Particle Markov
chain Monte Carlo methods." *JRSS-B*, 72(3), p. 269–342.
[doi:10.1111/j.1467-9868.2009.00736.x](https://doi.org/10.1111/j.1467-9868.2009.00736.x)

## See also

[`compile.sde_model`](https://robustecologies.github.io/lucifer/reference/compile.sde_model.md),
[`plot.sde_model`](https://robustecologies.github.io/lucifer/reference/plot.sde_model.md),
[`print.sde_model`](https://robustecologies.github.io/lucifer/reference/print.sde_model.md),
[`simulate.sde_model`](https://robustecologies.github.io/lucifer/reference/simulate.sde_model.md),
[`summary.sde_model`](https://robustecologies.github.io/lucifer/reference/summary.sde_model.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate OU data with known parameters
set.seed(42)
n <- 200; dt <- 0.1
kappa <- 2; mu <- 5; sigma <- 0.8; obs_sd <- 0.3
x <- numeric(n); x[1] <- mu
for (i in 2:n) {
    x[i] <- x[i-1] + kappa * (mu - x[i-1]) * dt +
        sigma * sqrt(dt) * rnorm(1)
}
y <- x + rnorm(n, 0, obs_sd)
times <- seq(0, by = dt, length.out = n)

# Build model using the OU family
sde <- SDE(data = y, times = times, family = "ou")

# Inspect
print(sde)
summary(sde)

# Fit via MCMC (AMWG recommended for SDE models)
fit <- SDE.fit(sde, Algorithm = "AMWG",
    Iterations = 50000, Thinning = 10)
print(fit)
summary(fit)

# Posterior trajectory (filtered states tracking data)
plot(fit, type = "trajectory")

# Posterior predictive paths
plot(fit, type = "predictive")

# Forecast 30 steps ahead
pred <- predict(fit, n.ahead = 30)
plot(pred)
} # }
```
