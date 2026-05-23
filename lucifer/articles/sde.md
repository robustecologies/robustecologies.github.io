# Bayesian inference for stochastic differential equations

## Introduction

Stochastic differential equations (SDEs) describe continuous-time
dynamical systems perturbed by noise. They appear whenever a
deterministic model, expressed as an ordinary differential equation,
must be augmented with random forcing to capture environmental
variability, measurement error, demographic stochasticity, or model
misspecification. The canonical form is

\\dX_t = f(X_t, \theta)\\dt + g(X_t, \theta)\\dW_t\\

where \\X_t \in \mathbb{R}^d\\ is the state vector, \\f: \mathbb{R}^d
\times \Theta \to \mathbb{R}^d\\ is the drift (the deterministic
skeleton of the dynamics), \\g: \mathbb{R}^d \times \Theta \to
\mathbb{R}^{d \times m}\\ is the diffusion coefficient (the noise
amplitude, which may depend on the state), \\\theta \in \Theta\\ is the
parameter vector, and \\W_t \in \mathbb{R}^m\\ is a standard Wiener
process [\[1\]](#ref1). Bayesian inference requires evaluating or
approximating the likelihood \\p(y\_{1:T} \mid \theta) =
\prod\_{k=1}^{T} p(y_k \mid y\_{1:k-1}, \theta)\\, which in turn
requires the transition density \\p(X\_{t\_{k+1}} \mid X\_{t_k},
\theta)\\. This transition density is analytically known only for a
handful of linear processes (Ornstein-Uhlenbeck, geometric Brownian
motion, Cox-Ingersoll-Ross); for general nonlinear SDEs, it must be
approximated [\[2\]](#ref2).

The R ecosystem offers several packages for SDE inference. The `sde`
package [\[3\]](#ref3) provides simulation and maximum-likelihood
estimation for one-dimensional diffusions with exact transition
densities but does not include MCMC machinery. The `yuima` package
[\[4\]](#ref4) handles quasi-likelihood estimation, simulation of
multidimensional diffusions, and some Bayesian methods but has a complex
object-oriented interface. The `pomp` package [\[5\]](#ref5) provides a
general partially-observed Markov process framework with iterated
filtering (IF2) and particle MCMC but does not integrate with a
general-purpose Bayesian inference engine. lucifer’s SDE module fills
the gap by offering multiple likelihood approximation strategies (exact
transition densities, Extended Kalman Filter, three particle filter
variants), plugging them directly into 55 MCMC algorithms, 8 variational
methods, SMC, and ABC with zero code changes, and providing SDE-aware
model comparison via LOO-PSIS, DIC, and robust Bayesian sensitivity
analysis [\[6\]](#ref6).

## Mathematical foundations

### Ito stochastic calculus

The mathematical foundation of SDEs rests on Ito’s stochastic calculus
[\[7\]](#ref7). An Ito SDE \\dX = f(X)dt + g(X)dW\\ defines a diffusion
process whose sample paths are continuous (by the continuity of the
Wiener process) but nowhere differentiable (by the unbounded variation
of Brownian motion). The Ito integral \\\int_0^t g(X_s)\\dW_s\\ is a
martingale with variance \\\int_0^t g(X_s)^2\\ds\\, and the fundamental
tool for manipulating SDEs is the Ito formula: for a twice continuously
differentiable function \\h(x, t)\\,

\\dh(X_t, t) = \left\[\frac{\partial h}{\partial t} + f\frac{\partial
h}{\partial x} + \frac{1}{2}g^2\frac{\partial^2 h}{\partial
x^2}\right\]dt + g\frac{\partial h}{\partial x}\\dW_t\\

The additional \\\frac{1}{2}g^2 h''\\ term (the Ito correction) has no
counterpart in ordinary calculus and arises from the quadratic variation
of Brownian motion: \\(dW)^2 = dt\\ in the Ito sense. This correction is
what makes, for example, the log-transform of geometric Brownian motion
yield a drift-corrected linear SDE: if \\dS = \mu S\\dt + \sigma
S\\dW\\, then \\d\log S = (\mu - \sigma^2/2)\\dt + \sigma\\dW\\
[\[1\]](#ref1).

For multidimensional systems \\dX = f(X)\\dt + G(X)\\dW\\ with \\X \in
\mathbb{R}^d\\ and \\G \in \mathbb{R}^{d \times m}\\, the Ito formula
generalizes to

\\dh = \left\[\frac{\partial h}{\partial t} + \sum_i f_i \frac{\partial
h}{\partial x_i} + \frac{1}{2}\sum\_{i,j} (GG^\top)\_{ij}
\frac{\partial^2 h}{\partial x_i \partial x_j}\right\]dt + \sum\_{i,k}
G\_{ik}\frac{\partial h}{\partial x_i}\\dW_k\\

The diffusion matrix \\\Sigma(x) = G(x)G(x)^\top\\ determines the local
covariance structure of the process and appears in the Fokker-Planck
equation governing the evolution of the marginal density [\[8\]](#ref8).

### Transition densities and the Fokker-Planck equation

The transition density \\p(x, t \mid x_0, t_0)\\ satisfies the forward
Kolmogorov (Fokker-Planck) equation

\\\frac{\partial p}{\partial t} = -\sum_i \frac{\partial}{\partial
x_i}\[f_i\\ p\] + \frac{1}{2}\sum\_{i,j}\frac{\partial^2}{\partial x_i
\partial x_j}\[\Sigma\_{ij}\\ p\]\\

with initial condition \\p(x, t_0 \mid x_0, t_0) = \delta(x - x_0)\\.
This PDE is analytically solvable only for a few processes: linear drift
with constant or state-proportional diffusion. For the
Ornstein-Uhlenbeck process \\dX = \kappa(\mu - X)dt + \sigma\\dW\\, the
transition density is Gaussian with mean \\\mu + (x_0 -
\mu)e^{-\kappa\Delta t}\\ and variance \\\frac{\sigma^2}{2\kappa}(1 -
e^{-2\kappa\Delta t})\\ [\[9\]](#ref9). For the Cox-Ingersoll-Ross
process \\dX = \kappa(\mu - X)dt + \sigma\sqrt{X}\\dW\\, the transition
density involves a non-central chi-squared distribution
[\[10\]](#ref10). For geometric Brownian motion \\dX = \mu X\\dt +
\sigma X\\dW\\, the log-price follows a Gaussian random walk.

### Numerical discretization

When the transition density is not available in closed form, it must be
approximated numerically. The two standard approaches are Euler-Maruyama
and Milstein [\[2\]](#ref2).

**Euler-Maruyama.** Given a step size \\\Delta t\\, the Euler-Maruyama
scheme approximates the transition as

\\X\_{k+1} = X_k + f(X_k)\Delta t + g(X_k)\sqrt{\Delta t}\\Z_k, \quad
Z_k \sim N(0, I_m)\\

The implied transition density is Gaussian: \\X\_{k+1} \mid X_k \sim
N(X_k + f(X_k)\Delta t,\\ G(X_k)G(X_k)^\top\Delta t)\\. This
approximation has weak order 1 and strong order 0.5 [\[2\]](#ref2). When
the inter-observation interval is large relative to the dynamics,
sub-stepping is necessary: the interval \\(t_k, t\_{k+1})\\ is divided
into \\n\_{sub} = \lceil (t\_{k+1} - t_k)/\delta t \rceil\\ sub-steps,
each of size \\\delta s = (t\_{k+1} - t_k)/n\_{sub}\\. lucifer
automatically computes the number of sub-steps from the user-specified
`dt` parameter.

**Milstein scheme.** The Milstein scheme adds a correction term
involving the derivative of the diffusion coefficient:

\\X\_{k+1} = X_k + f\Delta t + g\sqrt{\Delta t}\\Z +
\frac{1}{2}g\\\frac{\partial g}{\partial x}(\Delta t\\Z^2 - \Delta t)\\

This correction arises from the Ito-Taylor expansion of \\g(X)\\ and
achieves strong order 1 (vs. 0.5 for Euler-Maruyama), meaning the
pathwise approximation error scales as \\O(\Delta t)\\ rather than
\\O(\sqrt{\Delta t})\\ [\[2\]](#ref2). lucifer computes the diffusion
derivative \\\partial g / \partial x\\ via adaptive central finite
differences with step size \\h = \max(10^{-5}, 10^{-5}\|x\|)\\.

### Milstein scheme and strong convergence

The Milstein scheme achieves strong order 1 by including the Ito-Taylor
correction from the second-order term:

\\X\_{k+1} = X_k + f\Delta t + g\sqrt{\Delta t}\\Z +
\frac{1}{2}g\frac{\partial g}{\partial x}\left((\sqrt{\Delta t}\\Z)^2 -
\Delta t\right)\\

where \\Z \sim N(0,1)\\. The key quantity is \\g\\g'\\: for the CIR
process (\\g = \sigma\sqrt{x}\\), \\g' = \sigma/(2\sqrt{x})\\, so
\\g\\g' = \sigma^2/2\\, a constant that makes the Milstein correction
straightforward. For state-dependent diffusions like the FitzHugh-Nagumo
model, the correction involves nontrivial derivatives that lucifer
computes via adaptive central finite differences with step \\h =
\max(10^{-5}, 10^{-5}\|x\|)\\.

The distinction between weak and strong convergence is important for
practitioners: weak convergence (order 1 for both Euler and Milstein)
governs the accuracy of expectations \\E\[h(X_T)\]\\ and hence
transition density approximations. Strong convergence governs pathwise
accuracy \\E\[\|X_T - \hat{X}\_T\|\]\\. For likelihood-based inference,
weak order matters; for simulation and posterior predictive checks,
strong order matters. The Milstein scheme improves pathwise accuracy at
the cost of additional diffusion evaluations for \\g'\\.

### Extended Kalman Filter

For Gaussian observation models \\y_k \sim N(H x_k, R)\\, the Extended
Kalman Filter (EKF) provides a deterministic approximate likelihood by
propagating a Gaussian belief state through the nonlinear dynamics
[\[11\]](#ref11). Between observation times, the mean \\m\\ and
covariance \\P\\ of the filtering distribution evolve according to the
continuous-time Riccati-type equations

\\\dot{m} = f(m), \quad \dot{P} = J_f(m)\\P + P\\J_f(m)^\top +
G(m)G(m)^\top\\

where \\J_f = \partial f/\partial x \mid\_{x=m}\\ is the Jacobian of the
drift evaluated at the current mean estimate, and \\G\\ is the diffusion
matrix. lucifer discretizes these equations via Euler sub-stepping with
the same `dt` parameter used for the Euler-Maruyama scheme. At each
observation time \\t_k\\, the standard Kalman update applies:

\\S = H P H^\top + R, \quad K = P H^\top S^{-1}, \quad m \leftarrow m +
K(y_k - Hm), \quad P \leftarrow (I - KH)P(I - KH)^\top + KRK^\top\\

The covariance update uses the Joseph form, which is numerically more
stable than the simpler \\(I - KH)P\\ because it guarantees positive
semi-definiteness even under finite-precision arithmetic
[\[11\]](#ref11).

The log-likelihood increment is \\\ell_k = -\frac{1}{2}\left\[d_y
\log(2\pi) + \log\|S\| + v^\top S^{-1} v\right\]\\ where \\v = y_k -
Hm\\ is the innovation and \\d_y\\ is the observation dimension. The
total log-likelihood is \\\sum_k \ell_k\\.

The EKF has two key properties that make it attractive for MCMC. First,
it is deterministic: the same parameters always produce the same
log-likelihood, eliminating the noise that plagues pseudo-marginal
methods. Second, it is fast: each observation requires one Jacobian
evaluation (via \\2d\\ drift calls for central differences) plus matrix
operations, which is much cheaper than running \\N\\ particle
trajectories. The approximation is exact for linear SDEs (the OU
process) and deteriorates for strongly nonlinear systems where the
Gaussian assumption breaks down [\[11\]](#ref11).

### Kalman filtering for exact families with observation noise

When an exact family (OU, GBM, CIR) is observed with Gaussian noise, the
system is a linear-Gaussian state-space model and the correct marginal
likelihood \\p(y\_{1:T} \mid \theta)\\ requires integrating out the
latent states. lucifer uses an affine Kalman filter tailored to each
family. For the OU process, the state-space representation is

\\X\_{k+1} = e^{-\kappa\Delta t}\\X_k + \mu(1 - e^{-\kappa\Delta t}) +
w_k, \quad Y_k = X_k + v_k\\

where \\w_k \sim N(0, Q)\\ with \\Q = \frac{\sigma^2}{2\kappa}(1 -
e^{-2\kappa\Delta t})\\ and \\v_k \sim N(0, R)\\ with \\R =
\sigma\_{\text{obs}}^2\\. The Kalman filter recursion tracks the
filtered state estimate \\\hat{X}\_k = E\[X_k \mid Y\_{1:k}\]\\ and its
uncertainty \\P_k = \text{Var}(X_k \mid Y\_{1:k})\\, yielding the exact
marginal log-likelihood as the sum of innovation log-densities. For CIR,
whose transition variance depends on the state, a quasi-Kalman filter
plugs the filtered mean into the variance formula. For GBM with
observation noise, an EKF in log-space handles the nonlinear observation
map \\h(z) = \exp(z)\\.

When observation noise is zero or negligible (\< \\10^{-10}\\), lucifer
uses the direct closed-form transition densities, since observations are
then identical to the latent states and no filtering is needed.

### Particle filtering

When the observation model is non-Gaussian (Poisson counts, binomial
proportions, custom likelihoods) or the dynamics are highly nonlinear,
the particle filter provides an unbiased estimate of the marginal
likelihood [\[12\]](#ref12). The bootstrap particle filter (Gordon et
al. 1993 [\[13\]](#ref13)) represents the filtering distribution by a
weighted particle cloud \\\\x_k^{(i)}, w_k^{(i)}\\\_{i=1}^N\\. At each
observation, particles are propagated forward via Euler-Maruyama (the
prior proposal), weighted by the observation likelihood \\w_k^{(i)}
\propto p(y_k \mid x_k^{(i)})\\, and resampled when the effective sample
size \\\text{ESS} = 1/\sum_i (w_k^{(i)})^2\\ drops below a threshold
(default: \\0.5N\\). The log-marginal-likelihood estimate

\\\log \hat{p}(y\_{1:T} \mid \theta) = \sum\_{k=1}^{T}
\log\left(\frac{1}{N}\sum\_{i=1}^{N} w_k^{(i)}\right)\\

is unbiased in the sense that \\\mathbb{E}\[\hat{p}\] = p(y\_{1:T} \mid
\theta)\\, which ensures that the pseudo-marginal MCMC chain targets the
correct posterior [\[14\]](#ref14).

lucifer provides three particle filter variants accessible via the
`pf.type` parameter:

The **bootstrap** PF (`pf.type = "bootstrap"`) uses the prior dynamics
as the proposal distribution. It is robust and unbiased but can be
inefficient when observations are highly informative, since particles
propagated under the prior may drift far from the observed data.

The **auxiliary** PF (`pf.type = "auxiliary"`, Pitt and Shephard 1999
[\[15\]](#ref15)) adds a first-stage look-ahead step. Before
propagating, each particle receives a preliminary weight based on how
well a one-step deterministic prediction matches the next observation.
Particles with higher predictive likelihood are resampled more
frequently before propagation, improving ESS when the data are
informative.

The **bridge** PF (`pf.type = "bridge"`, Delyon and Hu 2006
[\[16\]](#ref16)) modifies the drift during particle propagation to
guide particles toward the next observation:

\\\tilde{f}(x, t) = f(x, t) + \frac{g(x,t)^2\\(y\_{k+1} - x)}{t\_{k+1} -
t}\\

This Brownian-bridge-like proposal dramatically improves efficiency for
sparse observations where the inter-observation interval is large
relative to the diffusion time scale. The importance weight corrects for
the modified proposal to maintain unbiasedness.

### Pseudo-marginal MCMC

The key insight of Andrieu et al. (2010) [\[14\]](#ref14) is that
replacing the exact likelihood in a Metropolis-Hastings algorithm with
an unbiased estimate preserves the correct stationary distribution. This
means any particle filter can be used within any MCMC algorithm, and the
resulting chain targets the exact posterior \\p(\theta \mid y\_{1:T})\\
despite the stochastic likelihood estimate. The efficiency depends on
the variance of the log-likelihood estimate: high variance leads to low
acceptance rates and poor mixing. lucifer’s
[`SDE.fit.pmmh()`](https://robustecologies.github.io/lucifer/reference/SDE.fit.pmmh.md)
wrapper sets appropriate defaults (HARM algorithm, 500 particles) and
allows users to tune the number of particles and choose between
bootstrap, auxiliary, and bridge proposals.

### Numerical stability considerations

Several numerical issues arise in SDE discretization that lucifer
handles automatically:

**Positivity preservation.** The CIR and Lotka-Volterra processes should
remain positive, but Euler-Maruyama can produce negative values. lucifer
clamps states at \\10^{-6}\\ in the built-in family drift/diffusion
functions.

**Sub-stepping.** When the inter-observation interval \\\Delta
t\_{obs}\\ is large relative to the dynamics, the Euler-Maruyama
approximation degrades. lucifer automatically sub-steps: \\n\_{sub} =
\lceil \Delta t\_{obs}/\delta t \rceil\\ sub-steps of size \\\delta s =
\Delta t\_{obs}/n\_{sub}\\, where \\\delta t\\ defaults to \\\min(\Delta
t\_{obs})/10\\.

**Log-sum-exp trick.** The particle filter computes
\\\log\left(\frac{1}{N}\sum_i \exp(w_i)\right)\\ for the incremental
log-normalizing constant. Direct computation overflows for large
\\w_i\\. lucifer uses \\\max(w) + \log\left(\sum_i \exp(w_i -
\max(w))\right) - \log N\\.

**Covariance symmetry.** The EKF covariance matrix \\P\\ accumulates
numerical asymmetry through the Riccati equation. lucifer enforces
symmetry via \\P \leftarrow \frac{1}{2}(P + P^\top)\\ after each
sub-step.

## The `SDE()` constructor

The primary entry point for SDE modelling in lucifer is the
[`SDE()`](https://robustecologies.github.io/lucifer/reference/SDE.md)
function. It constructs an object of class `sde_model` that encapsulates
the full model specification: dynamics, observation model, priors,
parameter constraints, and inference method. The generated
`Model(parm, Data)` closure plugs directly into
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`SMC()`](https://robustecologies.github.io/lucifer/reference/SMC.md),
[`ABC()`](https://robustecologies.github.io/lucifer/reference/ABC.md),
and all other lucifer inference backends.

### Using built-in families

The simplest way to specify an SDE is via the `family` argument, which
selects from 11 pre-built model templates:

| Family          | Equation                     | Method   | Parameters |
|:----------------|:-----------------------------|:---------|:----------:|
| ou              | dX = k(m-X)dt + s dW         | exact    |     4      |
| vasicek         | dX = k(m-X)dt + s dW         | exact    |     4      |
| gbm             | dX = mX dt + sX dW           | exact    |     2      |
| cir             | dX = k(m-X)dt + s sqrt(X) dW | exact    |     4      |
| lotka_volterra  | 2D predator-prey             | particle |     5      |
| sir             | 2D SIR epidemic              | particle |     3      |
| fitzhugh_nagumo | 2D neural oscillator         | particle |     6      |
| heston          | 2D stochastic volatility     | particle |     5      |
| double_well     | dX = a(X-X^3)dt + s dW       | ekf      |     3      |
| merton          | GBM + log-normal jumps       | euler    |     5      |
| kou             | GBM + double-exp jumps       | euler    |     6      |

Built-in SDE families {.table}

Each family supplies default drift, diffusion, prior, parameter names,
constraints, initial values, and the recommended inference method.
User-supplied arguments override the defaults. For example:

``` r

# Use the OU family but with a custom prior
sde <- SDE(data = y, times = times, family = "ou",
           prior = function(theta) {
               if (theta[1] <= 0 || theta[3] <= 0 || theta[4] <= 0)
                   return(-Inf)
               dnorm(theta[1], 1, 0.5, log = TRUE) +   # informative on kappa
                   dnorm(theta[2], 0, 5, log = TRUE) +
                   dexp(theta[3], 2, log = TRUE) +
                   dexp(theta[4], 10, log = TRUE)
           })
```

### Custom SDE specification

For models not covered by the built-in families, supply drift and
diffusion functions directly. Each function takes the state vector `x`,
the parameter vector `theta`, and the current time `t`, and returns a
vector of the same dimension as the state.

``` r

# Damped nonlinear oscillator with multiplicative noise
sde_osc <- SDE(
    drift = function(x, theta, t) {
        omega <- theta[1]; gamma <- theta[2]; beta <- theta[3]
        c(x[2],                              # dx1 = x2
          -omega^2 * x[1] - gamma * x[2] -   # dx2 = -w^2 x1 - g x2 - b x1^3
              beta * x[1]^3)
    },
    diffusion = function(x, theta, t) {
        c(0, theta[4])  # noise only on velocity component
    },
    data = y_osc, times = times_osc, x0 = c(1, 0),
    obs.model = "gaussian",
    obs.noise = function(theta) theta[5],
    obs.link = function(x, theta) x[1],  # observe position only
    prior = function(theta) {
        if (any(theta[c(1,2,4,5)] <= 0)) return(-Inf)
        sum(dnorm(theta, 0, 5, log = TRUE))
    },
    parm.names = c("omega", "gamma", "beta", "sigma", "obs_sd"),
    state.names = c("position", "velocity"),
    method = "ekf"
)
```

### Observation models

lucifer supports four observation models. The `"gaussian"` model \\y_k
\sim N(h(X_k, \theta), \sigma\_{obs}^2)\\ is the default and compatible
with all inference methods. The `"poisson"` model \\y_k \sim
\text{Pois}(\lambda_k)\\ where \\\lambda_k = h(X_k, \theta)\\ is the
rate via `obs.link`, and the `"binomial"` model \\y_k \sim \text{Bin}(n,
p_k)\\ where \\p_k = h(X_k, \theta)\\, both require the particle filter
since the EKF assumes Gaussian observations. The `"custom"` model allows
any likelihood via a user-supplied function `obs.loglik(y_k, x, theta)`
returning a scalar log-likelihood.

The `obs.link` function maps the latent state to the observation space.
When omitted for Gaussian observations, the identity is assumed (direct
observation of the first `dim_y` state components). When present, it
allows partial observation (e.g., observing voltage but not recovery in
the FitzHugh-Nagumo model) or nonlinear observation operators.

### Jump-diffusion processes

Jump-diffusion SDEs of the form \\dX = f(X)dt + g(X)dW +
J\\dN(\lambda)\\ extend continuous diffusions with discontinuous jumps
arriving at Poisson rate \\\lambda\\. lucifer supports jump-diffusion by
specifying `jump.rate` (a function returning the Poisson intensity) and
`jump.size` (a function returning a random jump vector). Jumps are
integrated into the Euler-Maruyama, Milstein, and particle filter
engines. The built-in Merton and Kou families are jump-diffusions.

``` r

# Mean-reverting process with occasional regime-switching jumps
sde_jump <- SDE(
    drift = function(x, theta, t) theta[1] * (theta[2] - x),
    diffusion = function(x, theta, t) theta[3],
    data = y, times = times, x0 = y[1],
    obs.model = "gaussian",
    obs.noise = function(theta) theta[6],
    prior = function(theta) {
        if (theta[1] <= 0 || theta[3] <= 0 || theta[4] <= 0 ||
            theta[5] <= 0 || theta[6] <= 0) return(-Inf)
        sum(dnorm(theta, 0, 10, log = TRUE))
    },
    parm.names = c("kappa", "mu", "sigma", "lambda",
                   "jump_sd", "obs_sd"),
    method = "euler",
    jump.rate = function(x, theta, t) theta[4],
    jump.size = function(x, theta, t) rnorm(1, 0, theta[5])
)
```

### Parameter constraints

SDE parameters often have natural constraints. The `constraints`
argument enforces them via
[`interval()`](https://robustecologies.github.io/lucifer/reference/interval.md)
at every likelihood evaluation:

``` r

sde_con <- SDE(
    drift = function(x, theta, t) theta[1] * (theta[2] - x),
    diffusion = function(x, theta, t) theta[3],
    data = y_ou, times = times, x0 = y_ou[1],
    obs.model = "gaussian", obs.noise = function(theta) theta[4],
    prior = function(theta) sum(dnorm(theta, 0, 10, log = TRUE)),
    parm.names = c("kappa", "mu", "sigma", "obs_sd"),
    method = "ekf",
    constraints = list(kappa = c(0, Inf), sigma = c(0, Inf),
                       obs_sd = c(0, Inf)))
cat("Auto initial values:", sde_con$Initial.Values, "\n")
```

### Compiled SDE expressions

lucifer can compile user-supplied drift and diffusion expressions from R
to native C++, eliminating the overhead of R callbacks in the particle
filter and EKF engines. This is achieved by specifying drift and
diffusion as character strings rather than R functions.

``` r

# Expression-based (auto-compiled to C++)
sde_compiled <- SDE(
    drift = "theta[1] * (theta[2] - x[1])",
    diffusion = "theta[3]",
    data = y_ou, times = times, x0 = y_ou[1],
    obs.model = "gaussian",
    obs.noise = function(theta) theta[4],
    prior = function(theta) {
        if (theta[1] <= 0 || theta[3] <= 0 || theta[4] <= 0)
            return(-Inf)
        sum(dnorm(theta, 0, 10, log = TRUE))
    },
    parm.names = c("kappa", "mu", "sigma", "obs_sd"),
    method = "ekf",
    compile = TRUE
)

# Verify compilation
cat("Compiled:", isTRUE(sde_compiled$Data$.sde_compiled))
```

The expression parser translates R syntax to C++ at model construction
time via
[`Rcpp::sourceCpp()`](https://rdrr.io/pkg/Rcpp/man/sourceCpp.html).
Supported operations include arithmetic, `exp`, `log`, `sqrt`, `sin`,
`cos`, `abs`, `pmax`, `pmin`, `^`, and subscripting (`x[i]`,
`theta[j]`). The compiled code is cached by hash, so repeated calls with
the same expression skip recompilation. For particle filter models, the
compiled drift and diffusion enable OpenMP-parallel particle propagation
using a thread-safe xoshiro256+ PRNG.

Alternatively, explicit compilation can be triggered on any
expression-based model:

``` r

sde <- SDE(drift = "theta[1] * (theta[2] - x[1])",
           diffusion = "theta[3]", ...)
sde <- compile(sde)  # S3 method
```

For built-in families (`family = "ou"`, etc.), native C++ dispatch is
automatic via pre-compiled drift/diffusion headers. No explicit
compilation step is needed.

### Stratonovich auto-conversion

lucifer’s SDE engines use Ito calculus internally, but models can be
specified in either Ito or Stratonovich form. The `interpretation`
parameter triggers automatic symbolic conversion when the Stratonovich
form is used. When drift and diffusion are specified as expressions
(character strings or
[`quote()`](https://rdrr.io/r/base/substitute.html)), lucifer’s symbolic
differentiation engine computes \\g'(x)\\ automatically and adds the
correction to the drift. The user specifies
`interpretation = "stratonovich"` and the rest is handled internally.
Full details and worked examples appear in the appendix on Ito
vs. Stratonovich interpretation.

## Practical considerations for MCMC convergence

SDE models typically present strong posterior correlations between
process parameters (e.g., the mean-reversion rate \\\kappa\\ and
diffusion coefficient \\\sigma\\ in the OU process are jointly
constrained by the process variance \\\sigma^2/(2\kappa)\\). For this
reason, the AMWG (Adaptive Metropolis-within-Gibbs) algorithm is
generally recommended over multivariate proposals like HARM, as its
componentwise adaptation handles these correlations more effectively. A
minimum of 50,000 iterations with thinning of 10 is advisable for most
SDE models. Always check convergence via trace plots, effective sample
size (ESS), and the Gelman-Rubin diagnostic when running multiple
chains.

For models using the Euler-Maruyama or Milstein discretization, the
log-likelihood is stochastic (a single simulated trajectory produces a
different value at every evaluation) and should be considered a rough
approximation. The EKF (for Gaussian observations) or particle filter
(for general observations) are strongly preferred for MCMC inference.

## Built-in SDE families

### Ornstein-Uhlenbeck process

\\dX = \kappa(\mu - X)\\dt + \sigma\\dW\\

The transition density is Gaussian:

\\X(t) \mid X(s) \sim N\left(\mu + (X(s) - \mu)e^{-\kappa\Delta t},\\
\frac{\sigma^2}{2\kappa}(1 - e^{-2\kappa\Delta t})\right)\\

The stationary distribution is \\X\_\infty \sim N(\mu,
\sigma^2/(2\kappa))\\. The autocorrelation function is \\\rho(\tau) =
e^{-\kappa\|\tau\|}\\, so the decorrelation time is \\1/\kappa\\. The
spectral density is Lorentzian: \\S(f) = \sigma^2/\[2\pi(\kappa^2 +
4\pi^2 f^2)\]\\.

This example demonstrates the complete workflow: simulate data with
known parameters, fit the model, inspect diagnostics, generate posterior
predictive checks, and forecast.

``` r

# Simulate from the OU process with known parameters
set.seed(42)
true_kappa <- 2; true_mu <- 3; true_sigma <- 0.5; true_obs_sd <- 0.1
times <- seq(0, 10, by = 0.1)
n <- length(times)

x_true <- numeric(n); x_true[1] <- true_mu
for (k in 2:n) {
    dt <- times[k] - times[k - 1]
    ekt <- exp(-true_kappa * dt)
    x_true[k] <- rnorm(1,
        true_mu + (x_true[k - 1] - true_mu) * ekt,
        true_sigma * sqrt((1 - exp(-2 * true_kappa * dt)) /
            (2 * true_kappa)))
}
y_ou <- x_true + rnorm(n, 0, true_obs_sd)
```

``` r

# Construct model and fit
sde_ou <- SDE(data = y_ou, times = times, family = "ou")
print(sde_ou)
summary(sde_ou)

# Plot observed data and prior predictive paths
plot(sde_ou, type = "data")
plot(sde_ou, type = "prior_predictive", n.paths = 100)

# Fit via MCMC (AMWG recommended for SDE models: better mixing with
# correlated parameters like kappa and sigma)
fit_ou <- SDE.fit(sde_ou, Algorithm = "AMWG",
                   Iterations = 50000, Status = 5000, Thinning = 10)
summary(fit_ou)
```

``` r

# MCMC diagnostics: trace plots, density, autocorrelation
plot(fit_ou, type = "posterior")

# Trajectory: posterior mean + 95% credible bands overlaid on data
plot(fit_ou, type = "trajectory")

# Posterior predictive paths
plot(fit_ou, type = "predictive", n.paths = 100)

# One-step-ahead prediction residuals
plot(fit_ou, type = "residuals")
```

The trajectory plot shows the posterior mean of the filtered state (blue
line) with 95% credible bands (shaded region) overlaid on the observed
data (black points). The predictive plot shows 30 forward-simulated
paths drawn from the posterior, providing a visual posterior predictive
check. The residual plot shows the difference between each observation
and the one-step-ahead prediction from the posterior mean parameters;
well-calibrated residuals should be roughly centered at zero with
constant variance.

``` r

# Forecast 30 steps ahead with prediction intervals
pred_ou <- predict(fit_ou, n.ahead = 30, n.paths = 200)
print(pred_ou)
plot(pred_ou)
```

The prediction intervals widen over time as process uncertainty
accumulates. For the OU process, the variance converges to the
stationary variance \\\sigma^2/(2\kappa)\\ as the forecast horizon
grows, so the intervals stabilize rather than growing without bound.

### Cox-Ingersoll-Ross process

\\dX = \kappa(\mu - X)\\dt + \sigma\sqrt{X}\\dW\\

The transition density involves a non-central chi-squared distribution.
Defining \\c = 2\kappa/\[\sigma^2(1-e^{-\kappa\Delta t})\]\\, \\q =
2cX(s)e^{-\kappa\Delta t}\\ (non-centrality), and \\\nu =
4\kappa\mu/\sigma^2\\ (degrees of freedom), the transition is \\2cX(t)
\sim \chi^2\_\nu(q)\\. The Feller condition \\2\kappa\mu \geq \sigma^2\\
ensures \\X(t) \> 0\\ a.s. The stationary distribution is Gamma:
\\X\_\infty \sim \text{Gamma}(2\kappa\mu/\sigma^2, 2\kappa/\sigma^2)\\.

The Cox-Ingersoll-Ross process is the canonical model for interest rates
and other positive-valued mean-reverting quantities [\[10\]](#ref10).
Its transition density involves a non-central chi-squared distribution,
which lucifer evaluates via R’s `dnchisq` function.

``` r

set.seed(123)
true_k <- 1.5; true_m <- 2; true_s <- 0.4; true_osd <- 0.05
times_cir <- seq(0, 15, by = 0.15)
n_cir <- length(times_cir)

x_cir <- numeric(n_cir); x_cir[1] <- true_m
for (k in 2:n_cir) {
    dt <- times_cir[k] - times_cir[k - 1]
    ekt <- exp(-true_k * dt)
    c_val <- 2 * true_k / (true_s^2 * (1 - ekt))
    q <- 2 * c_val * x_cir[k - 1] * ekt
    nu <- 4 * true_k * true_m / true_s^2
    x_cir[k] <- rchisq(1, df = nu, ncp = q) / (2 * c_val)
    x_cir[k] <- max(x_cir[k], 1e-6)
}
y_cir <- x_cir + rnorm(n_cir, 0, true_osd)
y_cir <- pmax(y_cir, 0.01)
```

``` r

sde_cir <- SDE(data = y_cir, times = times_cir, family = "cir")
fit_cir <- SDE.fit(sde_cir, Algorithm = "HARM",
                    Iterations = 20000, Status = 1000, Thinning = 5)
summary(fit_cir)
plot(fit_cir, type = "trajectory")
plot(sde_cir, type = "prior_predictive", n.paths = 100)
plot(fit_cir, type = "predictive", n.paths = 100)
```

The Feller condition \\2\kappa\mu \geq \sigma^2\\ determines whether the
CIR process stays strictly positive. When violated, the process can
reach zero with positive probability. lucifer’s CIR family constructor
issues a warning when the initial parameter values violate this
condition, helping the user diagnose potential boundary issues before
fitting.

### Double-well potential

The double-well potential \\dX = \alpha(X - X^3)dt + \sigma\\dW\\ has
two stable equilibria at \\X = \pm 1\\ separated by an unstable
equilibrium at \\X = 0\\. Noise-induced transitions between the wells
produce bistable dynamics that are characteristic of many physical and
biological systems [\[17\]](#ref17). The EKF handles this mildly
nonlinear system efficiently.

``` r

set.seed(99)
true_alpha <- 1; true_sig <- 0.6; true_osd2 <- 0.1
times_dw <- seq(0, 20, by = 0.05)
n_dw <- length(times_dw)

x_dw <- numeric(n_dw); x_dw[1] <- 0.5
for (k in 2:n_dw) {
    dt <- times_dw[k] - times_dw[k - 1]
    x_dw[k] <- x_dw[k - 1] +
        true_alpha * (x_dw[k - 1] - x_dw[k - 1]^3) * dt +
        true_sig * sqrt(dt) * rnorm(1)
}
y_dw <- x_dw + rnorm(n_dw, 0, true_osd2)
```

``` r

sde_dw <- SDE(data = y_dw, times = times_dw, family = "double_well")
cat("Method:", sde_dw$method, "\n")  # should be "ekf"

fit_dw <- SDE.fit(sde_dw, Algorithm = "HARM",
                   Iterations = 20000, Status = 1000, Thinning = 5)
summary(fit_dw)
plot(fit_dw, type = "trajectory")
plot(sde_dw, type = "prior_predictive", n.paths = 100)
plot(fit_dw, type = "predictive", n.paths = 100)
```

Because the EKF is deterministic, each MCMC iteration produces the same
log-likelihood for the same parameter values. This eliminates the noisy
likelihood estimates that degrade mixing in pseudo-marginal methods. The
cost per iteration is comparable to the Euler-Maruyama stochastic engine
but without the acceptance rate penalties caused by likelihood noise.

### Geometric Brownian motion

\\dS = \mu S\\dt + \sigma S\\dW\\

By Ito’s formula, \\d\log S = (\mu - \sigma^2/2)dt + \sigma\\dW\\. The
transition density in log-space is

\\\log S(t) \mid \log S(s) \sim N\left(\log S(s) + (\mu -
\tfrac{1}{2}\sigma^2)\Delta t,\\ \sigma^2\Delta t\right)\\

The expected value is \\E\[S(t)\] = S(0)e^{\mu t}\\ (exponential
growth), and the variance is \\\text{Var}\[S(t)\] = S(0)^2 e^{2\mu
t}(e^{\sigma^2 t} - 1)\\. GBM has no stationary distribution; it is
non-stationary.

Geometric Brownian motion is the standard model for asset prices
[\[19\]](#ref19). The log-price follows a random walk with drift:
\\d\log S = (\mu - \sigma^2/2)dt + \sigma\\dW\\, a consequence of Ito’s
formula applied to the exponential transformation.

``` r

set.seed(77)
true_mu_gbm <- 0.1; true_sig_gbm <- 0.25
times_gbm <- seq(0, 1, by = 1/252)
n_gbm <- length(times_gbm)
log_s <- numeric(n_gbm); log_s[1] <- log(100)
for (k in 2:n_gbm) {
    dt <- times_gbm[k] - times_gbm[k - 1]
    log_s[k] <- log_s[k - 1] +
        (true_mu_gbm - 0.5 * true_sig_gbm^2) * dt +
        true_sig_gbm * sqrt(dt) * rnorm(1)
}
prices <- exp(log_s)
```

``` r

sde_gbm <- SDE(data = prices, times = times_gbm, family = "gbm")
fit_gbm <- SDE.fit(sde_gbm, Algorithm = "HARM",
                    Iterations = 15000, Thinning = 3)
summary(fit_gbm)
plot(fit_gbm, type = "trajectory")
plot(sde_gbm, type = "prior_predictive", n.paths = 100)
plot(fit_gbm, type = "predictive", n.paths = 100)

# Forecast 60 trading days
pred_gbm <- predict(fit_gbm, n.ahead = 60, n.paths = 500)
plot(pred_gbm)
```

### Merton jump-diffusion

\\dS = \mu S\\dt + \sigma S\\dW + S(e^J - 1)\\dN(\lambda), \quad J \sim
N(\mu_J, \sigma_J^2)\\

The characteristic function of log-returns over interval \\\Delta t\\ is

\\\phi(u) = \exp\left\[iu(\mu - \tfrac{1}{2}\sigma^2 -
\lambda(e^{\mu_J + \sigma_J^2/2} - 1))\Delta t - \tfrac{1}{2}\sigma^2
u^2 \Delta t + \lambda\Delta t(e^{iu\mu_J - \sigma_J^2 u^2/2} -
1)\right\]\\

The excess kurtosis relative to GBM is \\\lambda\Delta t(\mu_J^4 +
6\mu_J^2\sigma_J^2 + 3\sigma_J^4)/(\sigma^2\Delta t + \lambda\Delta
t(\mu_J^2 + \sigma_J^2))^2\\, which is always positive, producing the
heavy tails observed in financial returns.

The Merton model [\[19\]](#ref19) extends geometric Brownian motion with
log-normal jumps arriving at Poisson rate \\\lambda\\. The jump
component generates heavy-tailed returns and discontinuous price paths.
The model has 5 parameters: drift \\\mu\\, diffusion \\\sigma\\, jump
intensity \\\lambda\\, jump mean \\\mu_J\\, and jump standard deviation
\\\sigma_J\\.

``` r

set.seed(123)
true_mu_m <- 0.05; true_sig_m <- 0.15; true_lam <- 2
true_jmu <- -0.02; true_jsd <- 0.08
times_m <- seq(0, 2, by = 1/252)
n_m <- length(times_m)
S <- numeric(n_m); S[1] <- 100
for (k in 2:n_m) {
    dt_m <- times_m[k] - times_m[k - 1]
    dW <- rnorm(1, 0, sqrt(dt_m))
    n_jumps <- rpois(1, true_lam * dt_m)
    J_sum <- if (n_jumps > 0)
        sum(rnorm(n_jumps, true_jmu, true_jsd)) else 0
    S[k] <- S[k - 1] * exp((true_mu_m - 0.5 * true_sig_m^2) * dt_m +
        true_sig_m * dW + J_sum)
}
```

``` r

# Log-returns distribution: heavier tails than normal
log_ret <- diff(log(S))
ret_df <- data.frame(returns = log_ret)
ggplot2::ggplot(ret_df, ggplot2::aes(x = .data$returns)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
        bins = 80, fill = "#2196F3", alpha = 0.6, color = "white") +
    ggplot2::stat_function(fun = dnorm,
        args = list(mean = mean(log_ret), sd = sd(log_ret)),
        color = "#FF5722", linewidth = 0.8, linetype = "dashed") +
    ggplot2::labs(x = "Log-return", y = "Density",
        title = "Merton jump-diffusion: heavy-tailed returns",
        subtitle = "Dashed line: fitted normal (fails to capture tails)") +
    theme_relab()
```

``` r

sde_mer <- SDE(data = S, times = times_m, family = "merton")
fit_mer <- SDE.fit(sde_mer, Algorithm = "HARM",
                    Iterations = 15000, Thinning = 5)
summary(fit_mer)
plot(fit_mer, type = "trajectory")
plot(fit_mer, type = "predictive", n.paths = 25)
```

### Heston stochastic volatility

The Heston model [\[21\]](#ref21) extends Black-Scholes [\[37\]](#ref37)
by making the volatility itself stochastic, following a CIR-type
mean-reverting process. The 2D system is

\\d\log S = \left(\mu - \frac{V}{2}\right)dt + \sqrt{V}\\dW_1, \quad dV
= \kappa(\theta_v - V)\\dt + \sigma_v\sqrt{V}\\dW_2\\

with \\\text{Cor}(dW_1, dW_2) = \rho\\. The parameters are: \\\mu\\
(drift), \\\kappa\\ (volatility mean-reversion speed), \\\theta_v\\
(long-run variance), \\\sigma_v\\ (volatility of volatility), and
\\\rho\\ (leverage effect, typically negative for equities). The Feller
condition \\2\kappa\theta_v \geq \sigma_v^2\\ ensures variance stays
positive.

``` r

set.seed(42)
true_mu_h <- 0.05; true_kappa_h <- 2; true_theta_v <- 0.04
true_sig_v <- 0.3; true_rho <- -0.7
times_h <- seq(0, 2, by = 1/252)
n_h <- length(times_h)
logS <- numeric(n_h); V_h <- numeric(n_h)
logS[1] <- log(100); V_h[1] <- true_theta_v
for (k in 2:n_h) {
    dt_h <- times_h[k] - times_h[k - 1]
    z1 <- rnorm(1); z2 <- true_rho * z1 + sqrt(1 - true_rho^2) * rnorm(1)
    v <- max(V_h[k - 1], 1e-8)
    logS[k] <- logS[k - 1] + (true_mu_h - v/2) * dt_h +
        sqrt(v * dt_h) * z1
    V_h[k] <- max(v + true_kappa_h * (true_theta_v - v) * dt_h +
        true_sig_v * sqrt(v * dt_h) * z2, 1e-8)
}
prices_h <- exp(logS)
```

``` r

# Price and volatility paths
heston_df <- data.frame(
    time = rep(times_h, 2),
    value = c(prices_h, sqrt(V_h) * sqrt(252) * 100),
    variable = rep(c("Price", "Annualized volatility (%)"),
        each = n_h)
)
ggplot2::ggplot(heston_df, ggplot2::aes(x = .data$time,
    y = .data$value)) +
    ggplot2::geom_line(color = "#2196F3", linewidth = 0.4) +
    ggplot2::facet_wrap(~variable, ncol = 1, scales = "free_y") +
    ggplot2::labs(x = "Time (years)", y = NULL,
        title = "Heston stochastic volatility: simulated paths") +
    theme_relab()
```

``` r

sde_heston <- SDE(data = as.matrix(prices_h), times = times_h,
                   family = "heston")
fit_heston <- SDE.fit.pmmh(sde_heston, N.particles = 100,
                             Iterations = 3000, Thinning = 2)
summary(fit_heston)
plot(fit_heston, type = "trajectory")
```

``` r

# QQ-plot: model vs. empirical return distribution
ret_h <- diff(log(prices_h))
post_h <- colMeans(if (!is.null(fit_heston$Posterior2)) fit_heston$Posterior2 else fit_heston$Posterior1)
# Simulate returns from fitted model
set.seed(1)
n_sim <- 5000
sim_ret <- rnorm(n_sim, post_h[1] * (1/252),
    sqrt(post_h[3]) * sqrt(1/252))

qq_df <- data.frame(
    empirical = sort(ret_h),
    theoretical = quantile(sim_ret,
        probs = seq_along(ret_h) / (length(ret_h) + 1))
)
ggplot2::ggplot(qq_df, ggplot2::aes(x = .data$theoretical,
    y = .data$empirical)) +
    ggplot2::geom_point(alpha = 0.3, size = 0.8, color = "#2196F3") +
    ggplot2::geom_abline(slope = 1, intercept = 0,
        linetype = "dashed", color = "#FF5722") +
    ggplot2::labs(x = "Fitted model quantiles",
        y = "Empirical quantiles",
        title = "QQ-plot: Heston vs. empirical returns") +
    theme_relab()
```

### FitzHugh-Nagumo neural oscillator

\\\frac{dV}{dt} = \frac{V - V^3/3 - W + I\_{ext}}{\varepsilon}, \quad
\frac{dW}{dt} = \varepsilon(V + a - bW)\\

with noise \\dV \mathrel{+}= (\sigma/\sqrt{\varepsilon})\\dW_1\\. The
nullclines are \\W = V - V^3/3 + I\_{ext}\\ (cubic, fast) and \\W = (V +
a)/b\\ (linear, slow). Their intersection determines the fixed point.
When the fixed point lies on the middle branch of the cubic nullcline,
the system exhibits oscillations (limit cycle). The parameter
\\\varepsilon \ll 1\\ enforces time-scale separation: \\V\\ is fast
(spikes) and \\W\\ is slow (recovery).

The FitzHugh-Nagumo model [\[22\]](#ref22) is a 2D simplification of the
Hodgkin-Huxley equations describing neural excitability. The system
exhibits limit cycles, excitable behavior, and noise-driven transitions
between resting and spiking states. The dynamics are

\\\frac{dV}{dt} = \frac{1}{\varepsilon}\left(V - \frac{V^3}{3} - W +
I\_{ext}\right), \quad \frac{dW}{dt} = \varepsilon(V + a - bW)\\

where \\V\\ is the membrane voltage (fast variable), \\W\\ is a recovery
variable (slow variable), \\\varepsilon\\ controls the time-scale
separation, \\a\\ and \\b\\ shape the nullclines, and \\I\_{ext}\\ is an
external current. Noise enters through the voltage equation with
amplitude \\\sigma/\sqrt{\varepsilon}\\. Only \\V\\ is observed with
Gaussian noise.

``` r

# Simulate FitzHugh-Nagumo with known parameters
set.seed(42)
true_eps <- 0.08; true_a <- 0.7; true_b <- 0.8
true_I <- 0.5; true_sig <- 0.3; true_osd <- 0.15
times_fhn <- seq(0, 30, by = 0.15)
n_fhn <- length(times_fhn)
V <- numeric(n_fhn); W <- numeric(n_fhn)
V[1] <- 0; W[1] <- 0
dt_fhn <- 0.001  # fine integration step
for (k in 2:n_fhn) {
    n_steps <- round((times_fhn[k] - times_fhn[k - 1]) / dt_fhn)
    v <- V[k - 1]; w <- W[k - 1]
    for (s in seq_len(n_steps)) {
        dv <- (v - v^3/3 - w + true_I) / true_eps
        dw <- true_eps * (v + true_a - true_b * w)
        v <- v + dv * dt_fhn +
            (true_sig / sqrt(true_eps)) * sqrt(dt_fhn) * rnorm(1)
        w <- w + dw * dt_fhn
    }
    V[k] <- v; W[k] <- w
}
y_fhn <- V + rnorm(n_fhn, 0, true_osd)
```

``` r

# Phase portrait: voltage vs. recovery (colored by time)
phase_df <- data.frame(V = V, W = W, time = times_fhn)
ggplot2::ggplot(phase_df, ggplot2::aes(x = .data$V, y = .data$W,
    color = .data$time)) +
    ggplot2::geom_path(linewidth = 0.4, alpha = 0.8) +
    ggplot2::scale_color_viridis_c(option = "plasma") +
    ggplot2::labs(x = "Voltage V", y = "Recovery W",
        color = "Time",
        title = "FitzHugh-Nagumo phase portrait") +
    theme_relab()
```

``` r

sde_fhn <- SDE(data = y_fhn, times = times_fhn,
               family = "fitzhugh_nagumo")
fit_fhn <- SDE.fit.pmmh(sde_fhn, N.particles = 150,
                          Iterations = 5000, Thinning = 3)
summary(fit_fhn)
```

``` r

# Full diagnostic suite
plot(fit_fhn, type = "trajectory")
plot(fit_fhn, type = "predictive", n.paths = 15)
plot(fit_fhn, type = "residuals")
plot(fit_fhn, type = "posterior")
```

The trajectory plot shows the posterior mean filtered voltage with 95%
credible bands overlaid on the noisy observations. The residuals should
be approximately centered at zero; systematic patterns indicate model
misspecification (e.g., wrong time-scale separation \\\varepsilon\\).

### Lotka-Volterra and N-species community

#### Lotka-Volterra predator-prey

For \\N\\ species with generalized interaction matrix:

\\\frac{dX_i}{dt} = X_i\left(r_i - \sum\_{j=1}^{N} A\_{ij} X_j\right) +
\sigma_i X_i\\dW_i\\

The interaction matrix \\A\\ encodes: \\A\_{ii} \> 0\\ (self-regulation
/ carrying capacity), \\A\_{ij} \> 0\\ for \\i \neq j\\ (competition or
predation), \\A\_{ij} \< 0\\ (mutualism or facilitation). The
deterministic fixed point satisfies \\A\mathbf{x}^\* = \mathbf{r}\\, so
\\\mathbf{x}^\* = A^{-1}\mathbf{r}\\ when \\A\\ is invertible. Stability
requires all eigenvalues of \\-\text{diag}(\mathbf{x}^\*) A\\ to have
negative real parts (community matrix criterion, May 1974
[\[30\]](#ref30)).

The stochastic Lotka-Volterra predator-prey model has nonlinear dynamics
with Poisson observations, requiring particle filtering for inference.
This example demonstrates the
[`SDE.fit.pmmh()`](https://robustecologies.github.io/lucifer/reference/SDE.fit.pmmh.md)
wrapper and the auxiliary particle filter.

``` r

set.seed(42)
times_lv <- seq(0, 20, by = 0.5)
n_lv <- length(times_lv)

# Simulate from the stochastic Lotka-Volterra system
alpha <- 0.5; beta <- 0.025; gamma <- 0.3
sigma1 <- 0.1; sigma2 <- 0.1
prey <- numeric(n_lv); pred_pop <- numeric(n_lv)
prey[1] <- 50; pred_pop[1] <- 30
dt_sim <- 0.01
for (k in 2:n_lv) {
    t_int <- times_lv[k] - times_lv[k - 1]
    n_steps <- round(t_int / dt_sim)
    x1 <- prey[k - 1]; x2 <- pred_pop[k - 1]
    for (s in seq_len(n_steps)) {
        x1 <- max(x1 + (alpha * x1 - beta * x1 * x2) * dt_sim +
            sigma1 * x1 * sqrt(dt_sim) * rnorm(1), 1e-3)
        x2 <- max(x2 + (beta * x1 * x2 - gamma * x2) * dt_sim +
            sigma2 * x2 * sqrt(dt_sim) * rnorm(1), 1e-3)
    }
    prey[k] <- rpois(1, max(x1, 0.1))
    pred_pop[k] <- rpois(1, max(x2, 0.1))
}
```

``` r

sde_lv <- SDE(data = cbind(prey, pred_pop), times = times_lv,
              family = "lotka_volterra", pf.type = "auxiliary")
print(sde_lv)

fit_lv <- SDE.fit.pmmh(sde_lv, N.particles = 300,
                         Iterations = 5000, pf.type = "auxiliary")
summary(fit_lv)

# Trajectory and phase portrait
plot(fit_lv, type = "trajectory")
plot(fit_lv, type = "phase")

# Filtered state estimates with ESS
plot(fit_lv, type = "filtered")
```

The filtered state plot shows the particle filter’s running estimate of
the latent prey and predator populations (top panels) and the effective
sample size over time (bottom panel). Low ESS regions indicate time
points where the particle cloud degenerated and resampling occurred.

#### N-species community ecology

The generalized Lotka-Volterra family supports arbitrary numbers of
interacting species. For \\N\\ species, the dynamics are

\\\frac{dX_i}{dt} = X_i\left(r_i - \sum\_{j=1}^{N} A\_{ij} X_j\right) +
\sigma_i X_i\\dW_i\\

where \\r_i\\ is the intrinsic growth rate of species \\i\\, \\A\_{ij}\\
is the interaction coefficient (competition, predation, mutualism), and
\\\sigma_i\\ is the noise amplitude. The parameter vector has \\N +
N^2 + N = N(N+2)\\ entries: growth rates, the interaction matrix
(row-major), and noise amplitudes.

``` r

# Simulate a 3-species system
set.seed(42)
times_c <- seq(0, 15, by = 0.5)
sp1 <- rpois(length(times_c), 50 + 20 * sin(0.5 * times_c))
sp2 <- rpois(length(times_c), 30 + 10 * cos(0.3 * times_c))
sp3 <- rpois(length(times_c), 20 + 15 * sin(0.7 * times_c + 1))

# Using the lotka_volterra family with 3 columns
sde_3sp <- SDE(data = cbind(sp1, sp2, sp3), times = times_c,
               family = "lotka_volterra")
cat("Parameters:", length(sde_3sp$parm.names), "(3+9+3=15)")
print(sde_3sp)
```

For named species, the `community` family reads column names from the
data matrix:

``` r

community_data <- cbind(Fox = sp1, Rabbit = sp2, Grass = sp3)
sde_eco <- SDE(data = community_data, times = times_c,
               family = "community")
cat("Species:", paste(sde_eco$state.names, collapse = ", "))
```

The `community` family uses hierarchical priors: log-normal on growth
rates, normal on interaction coefficients (allowing competition and
mutualism), and half-normal on noise amplitudes. For large \\N\\, the
\\N^2\\ interaction parameters may require regularizing priors or
structured parameterizations.

#### 3-species competitive community

This example demonstrates lucifer’s N-species community model on a
simulated 3-species competitive system. The interaction matrix \\A\\
determines the community structure: diagonal elements represent
self-regulation (carrying capacity effects), while off-diagonal elements
encode interspecific competition or facilitation [\[30\]](#ref30),
[\[31\]](#ref31).

``` r

# Simulate a 3-species competitive community
set.seed(77)
N_sp <- 3
true_r <- c(0.5, 0.3, 0.4)      # intrinsic growth rates
true_A <- matrix(c(
    0.01,  0.005, 0.002,    # species 1: self-reg, competes with 2 and 3
    0.003, 0.008, 0.004,    # species 2: weaker self-reg
    0.001, 0.003, 0.012     # species 3: strongest self-reg
), N_sp, N_sp, byrow = TRUE)
true_sig <- c(0.05, 0.04, 0.06)  # process noise

times_eco <- seq(0, 30, by = 0.5)
n_eco <- length(times_eco)
X_eco <- matrix(NA, n_eco, N_sp)
X_eco[1, ] <- c(50, 40, 30)  # initial populations
dt_eco <- 0.01

for (k in 2:n_eco) {
    n_sub <- round((times_eco[k] - times_eco[k - 1]) / dt_eco)
    x <- X_eco[k - 1, ]
    for (s in seq_len(n_sub)) {
        x <- pmax(x, 1e-3)
        drift <- x * (true_r - as.numeric(true_A %*% x))
        diff_v <- true_sig * x
        x <- x + drift * dt_eco + diff_v * sqrt(dt_eco) * rnorm(N_sp)
    }
    X_eco[k, ] <- pmax(x, 0.1)
}
# Poisson observations
Y_eco <- matrix(rpois(n_eco * N_sp, pmax(X_eco, 0.1)), n_eco, N_sp)
colnames(Y_eco) <- c("Prey", "Mesopredator", "Apex")
```

``` r

# Plot observed counts
eco_df <- data.frame(
    time = rep(times_eco, N_sp),
    count = as.vector(Y_eco),
    species = rep(colnames(Y_eco), each = n_eco)
)
ggplot2::ggplot(eco_df, ggplot2::aes(x = .data$time, y = .data$count,
    color = .data$species)) +
    ggplot2::geom_line(alpha = 0.7) +
    ggplot2::geom_point(size = 0.8, alpha = 0.5) +
    ggplot2::scale_color_manual(values = c(
        Prey = "#4CAF50", Mesopredator = "#FF9800", Apex = "#F44336")) +
    ggplot2::labs(x = "Time", y = "Count",
        title = "3-species community: observed Poisson counts",
        color = "Species") +
    theme_relab()
```

``` r

# Fit the community model
sde_eco <- SDE(data = Y_eco, times = times_eco, family = "community")
cat("Parameters:", sde_eco$dim_theta,
    "(3 growth + 9 interaction + 3 noise = 15)\n")
cat("State names:", paste(sde_eco$state.names, collapse = ", "), "\n")

fit_eco <- SDE.fit.pmmh(sde_eco, N.particles = 100,
                          Iterations = 3000, Thinning = 2)
summary(fit_eco)
```

``` r

# Filtered trajectories for all species
plot(fit_eco, type = "trajectory")
```

``` r

# Reconstruct posterior interaction matrix
posterior_eco <- fit_eco$Posterior2
if (is.null(posterior_eco)) posterior_eco <- fit_eco$Posterior1
A_post <- colMeans(posterior_eco[, 4:12])  # columns 4-12 are A entries
A_mat <- matrix(A_post, N_sp, N_sp, byrow = TRUE)
rownames(A_mat) <- colnames(A_mat) <- colnames(Y_eco)

# Heatmap of posterior mean interaction matrix
a_df <- expand.grid(
    To = colnames(Y_eco),
    From = colnames(Y_eco)
)
a_df$value <- as.vector(t(A_mat))
ggplot2::ggplot(a_df, ggplot2::aes(x = .data$From, y = .data$To,
    fill = .data$value)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(
        label = sprintf("%.4f", .data$value)),
        color = "white", size = 3.5) +
    ggplot2::scale_fill_gradient2(low = "#2196F3", mid = "grey90",
        high = "#F44336", midpoint = 0) +
    ggplot2::labs(x = "Effect of", y = "On",
        fill = "Interaction\nstrength",
        title = "Posterior mean interaction matrix") +
    theme_relab()
```

``` r

# Compare posterior estimates with true values
true_theta <- c(true_r, as.vector(t(true_A)), true_sig)
post_mean <- colMeans(posterior_eco)
post_sd <- apply(posterior_eco, 2, sd)

comp_df <- data.frame(
    Parameter = sde_eco$parm.names,
    True = round(true_theta, 4),
    `Posterior mean` = round(post_mean, 4),
    `Posterior SD` = round(post_sd, 4),
    `Covered` = ifelse(
        true_theta >= post_mean - 1.96 * post_sd &
        true_theta <= post_mean + 1.96 * post_sd,
        "Yes", "No"),
    check.names = FALSE
)
knitr::kable(comp_df, caption = "3-species community: parameter recovery")
```

The interaction matrix heatmap reveals the community structure: diagonal
elements (self-regulation) should be largest, and off-diagonal elements
show interspecific effects. Red cells indicate competition (positive
\\A\_{ij}\\, which reduces growth), while blue cells would indicate
facilitation (negative \\A\_{ij}\\).

## Inference strategies

### Comparison of EKF vs. particle filter

For the OU process (a linear SDE), the EKF and bootstrap particle filter
should produce similar results, but the EKF is faster and deterministic.
This section compares them:

``` r

# EKF
sde_ekf <- SDE(
    drift = function(x, theta, t) theta[1] * (theta[2] - x),
    diffusion = function(x, theta, t) theta[3],
    data = y_ou, times = times, x0 = y_ou[1],
    obs.model = "gaussian",
    obs.noise = function(theta) theta[4],
    prior = function(theta) {
        if (theta[1] <= 0 || theta[3] <= 0 || theta[4] <= 0)
            return(-Inf)
        sum(dnorm(theta, 0, 10, log = TRUE))
    },
    parm.names = c("kappa", "mu", "sigma", "obs_sd"),
    method = "ekf"
)

# Particle filter (same model)
sde_pf <- SDE(
    drift = function(x, theta, t) theta[1] * (theta[2] - x),
    diffusion = function(x, theta, t) theta[3],
    data = y_ou, times = times, x0 = y_ou[1],
    obs.model = "gaussian",
    obs.noise = function(theta) theta[4],
    obs.loglik = function(y_k, x, theta) {
        dnorm(y_k, x, theta[4], log = TRUE)
    },
    prior = function(theta) {
        if (theta[1] <= 0 || theta[3] <= 0 || theta[4] <= 0)
            return(-Inf)
        sum(dnorm(theta, 0, 10, log = TRUE))
    },
    parm.names = c("kappa", "mu", "sigma", "obs_sd"),
    method = "particle", N.particles = 200
)

# Compare likelihoods at true parameters
theta_true <- c(true_kappa, true_mu, true_sigma, true_obs_sd)
ekf_ll <- sde_ekf$Model(theta_true, sde_ekf$Data)$Monitor[1]
set.seed(1)
pf_ll <- sde_pf$Model(theta_true, sde_pf$Data)$Monitor[1]

# EKF is deterministic (same theta -> same LL)
r1 <- sde_ekf$Model(theta_true, sde_ekf$Data)$Monitor[1]
r2 <- sde_ekf$Model(theta_true, sde_ekf$Data)$Monitor[1]

knitr::kable(data.frame(
    Method = c("EKF (deterministic)", "Particle filter (N=200)"),
    `Log-likelihood` = round(c(ekf_ll, pf_ll), 4),
    Deterministic = c("Yes", "No"),
    check.names = FALSE
), caption = "EKF vs. particle filter likelihood comparison at true parameters")
```

### Algorithm recommendations

The following table summarizes the recommended MCMC algorithms for SDE
models based on the inference method:

| Scenario | Recommended | Avoid |
|:---|:---|:---|
| Exact likelihood, \< 10 parameters | HARM or AMWG | NUTS (wasteful) |
| EKF, \< 20 parameters | HARM or twalk | NUTS (expensive gradients) |
| Particle filter, any dimension | HARM via SDE.fit.pmmh() | Any gradient-based |
| Jump-diffusion, any method | HARM or RAM | NUTS |
| Quick exploration | RWM with short run | Gradient-based |

Algorithm recommendations for SDE models {.table}

### Multi-algorithm comparison

Different MCMC algorithms have different strengths for SDE models. This
example compares HARM, twalk, and AMWG on the same OU dataset using
[`SDE.fit()`](https://robustecologies.github.io/lucifer/reference/SDE.fit.md),
evaluating acceptance rate, effective sample size, and wall-clock time.

``` r

algos <- c("HARM", "twalk", "AMWG")
fits_algo <- list()
for (algo in algos) {
    t0 <- proc.time()[3]
    fits_algo[[algo]] <- SDE.fit(sde_ou, Algorithm = algo,
        Iterations = 10000, Thinning = 5)
    fits_algo[[algo]]$wall_time <- proc.time()[3] - t0
}
```

``` r

algo_results <- data.frame(
    Algorithm = algos,
    `Acceptance rate` = sapply(fits_algo,
        function(f) round(f$Acceptance.Rate, 3)),
    `Wall time (s)` = sapply(fits_algo,
        function(f) round(f$wall_time, 1)),
    `ESS (kappa)` = sapply(fits_algo, function(f) {
        post <- f$Posterior2
        if (is.null(post)) post <- f$Posterior1
        round(coda::effectiveSize(coda::as.mcmc(post[, 1])), 0)
    }),
    check.names = FALSE
)
knitr::kable(algo_results,
    caption = "Algorithm comparison on OU process")
```

``` r

# Trace plot comparison
trace_dfs <- lapply(algos, function(algo) {
    post <- fits_algo[[algo]]$Posterior2
    if (is.null(post)) post <- fits_algo[[algo]]$Posterior1
    data.frame(
        iteration = seq_len(nrow(post)),
        kappa = post[, 1],
        mu = post[, 2],
        algorithm = algo
    )
})
trace_df <- do.call(rbind, trace_dfs)

ggplot2::ggplot(trace_df, ggplot2::aes(x = .data$iteration,
    y = .data$kappa, color = .data$algorithm)) +
    ggplot2::geom_line(alpha = 0.5, linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = true_kappa, linetype = "dashed") +
    ggplot2::facet_wrap(~algorithm, ncol = 1) +
    ggplot2::scale_color_manual(values = c(
        HARM = "#2196F3", twalk = "#FF9800", AMWG = "#4CAF50")) +
    ggplot2::labs(x = "Iteration", y = expression(kappa),
        title = expression(paste("Trace plots of ", kappa,
            " across algorithms"))) +
    ggplot2::guides(color = "none") +
    theme_relab()
```

## Posterior diagnostics and model checking

### MCMC convergence diagnostics

SDE models present unique convergence challenges: the log-likelihood may
be noisy (particle filter), the parameter space may have ridges (e.g.,
\\\kappa\\ and \\\sigma\\ are correlated in the OU process), and the
observation noise and process noise may be confounded. This section
demonstrates diagnostics that go beyond standard trace plots.

``` r

# Joint posterior: kappa vs. sigma (typically correlated)
post_ou <- fit_ou$Posterior2
if (is.null(post_ou)) post_ou <- fit_ou$Posterior1
pairs_df <- data.frame(
    kappa = post_ou[, 1], mu = post_ou[, 2],
    sigma = post_ou[, 3], obs_sd = post_ou[, 4]
)

# Pairs plot with contours
ggplot2::ggplot(pairs_df, ggplot2::aes(x = .data$kappa,
    y = .data$sigma)) +
    ggplot2::geom_point(alpha = 0.1, size = 0.5, color = "#2196F3") +
    ggplot2::geom_density_2d(color = "#1565C0", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = true_kappa, linetype = "dashed") +
    ggplot2::geom_hline(yintercept = true_sigma, linetype = "dashed") +
    ggplot2::labs(x = expression(kappa), y = expression(sigma),
        title = expression(paste("Joint posterior: ", kappa,
            " vs. ", sigma)),
        subtitle = "Dashed: true values") +
    theme_relab()
```

``` r

# Running mean of each parameter (should stabilize)
n_post <- nrow(post_ou)
running_df <- data.frame(
    iteration = rep(seq_len(n_post), 4),
    running_mean = c(
        cumsum(post_ou[, 1]) / seq_len(n_post),
        cumsum(post_ou[, 2]) / seq_len(n_post),
        cumsum(post_ou[, 3]) / seq_len(n_post),
        cumsum(post_ou[, 4]) / seq_len(n_post)
    ),
    parameter = rep(c("kappa", "mu", "sigma", "obs_sd"),
        each = n_post),
    true_value = rep(c(true_kappa, true_mu, true_sigma, true_obs_sd),
        each = n_post)
)

ggplot2::ggplot(running_df, ggplot2::aes(x = .data$iteration,
    y = .data$running_mean)) +
    ggplot2::geom_line(color = "#2196F3", linewidth = 0.4) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = .data$true_value),
        linetype = "dashed", color = "black") +
    ggplot2::facet_wrap(~parameter, scales = "free_y") +
    ggplot2::labs(x = "Iteration", y = "Running mean",
        title = "MCMC convergence: running means",
        subtitle = "Dashed: true values") +
    theme_relab()
```

``` r

# Geweke z-scores per parameter
if (requireNamespace("coda", quietly = TRUE)) {
    geweke_z <- coda::geweke.diag(
        coda::as.mcmc(post_ou))$z
    geweke_df <- data.frame(
        parameter = names(geweke_z),
        z_score = as.numeric(geweke_z),
        significant = abs(as.numeric(geweke_z)) > 1.96
    )
    ggplot2::ggplot(geweke_df, ggplot2::aes(x = .data$parameter,
        y = .data$z_score, fill = .data$significant)) +
        ggplot2::geom_col(alpha = 0.8) +
        ggplot2::geom_hline(yintercept = c(-1.96, 1.96),
            linetype = "dashed", color = "#FF5722") +
        ggplot2::scale_fill_manual(values = c("FALSE" = "#4CAF50",
            "TRUE" = "#F44336")) +
        ggplot2::labs(x = NULL, y = "Geweke z-score",
            title = "Geweke convergence diagnostic",
            subtitle = "Values outside dashed lines indicate non-convergence",
            fill = "|z| > 1.96") +
        ggplot2::coord_flip() +
        theme_relab()
}
```

The Geweke diagnostic compares the mean of the first 10% of the chain
with the mean of the last 50%. A z-score exceeding \\\pm 1.96\\
indicates the chain has not converged at the 5% significance level. For
SDE models, the observation noise parameter \\\sigma\_{obs}\\ typically
converges fastest (it is directly identified by the residual variance),
while the reversion speed \\\kappa\\ converges slowest (it requires long
series to distinguish from unit-root behavior).

### Posterior predictive checks

Posterior predictive checking (PPC) is essential for assessing whether a
fitted SDE captures the key features of the observed data
[\[25\]](#ref25). For time series, relevant summaries include the
autocorrelation function, the marginal distribution, the distribution of
increments, and the running variance.

``` r

set.seed(42)
posterior <- fit_ou$Posterior2
if (is.null(posterior)) posterior <- fit_ou$Posterior1
n_rep <- 50
idx <- sample.int(nrow(posterior), n_rep, replace = TRUE)
reps <- list()
for (r in seq_len(n_rep)) {
    sim <- simulate(sde_ou, theta = posterior[idx[r], ], n.paths = 1)
    reps[[r]] <- sim$trajectories[, 1, 1]
}
```

``` r

require_package("ggplot2")
acf_obs <- acf(y_ou, lag.max = 20, plot = FALSE)$acf[, , 1]
acf_reps <- t(sapply(reps, function(r)
    acf(r, lag.max = 20, plot = FALSE)$acf[, , 1]))
acf_df <- data.frame(
    lag = rep(0:20, n_rep + 1),
    acf = c(acf_obs, as.vector(t(acf_reps))),
    group = c(rep("Observed", 21),
              rep(paste0("r", seq_len(n_rep)), each = 21)),
    type = c(rep("observed", 21), rep("replicate", 21 * n_rep)))
ggplot2::ggplot(acf_df, ggplot2::aes(x = .data$lag, y = .data$acf,
    group = .data$group, color = .data$type, alpha = .data$type)) +
    ggplot2::geom_line() +
    ggplot2::scale_color_manual(values = c(observed = "black",
        replicate = "#2196F3")) +
    ggplot2::scale_alpha_manual(values = c(observed = 1, replicate = 0.15)) +
    ggplot2::labs(x = "Lag", y = "ACF",
        title = "PPC: autocorrelation function") +
    ggplot2::guides(alpha = "none") + theme_relab()
```

``` r

dens_df <- data.frame(
    value = c(y_ou, unlist(reps)),
    type = c(rep("Observed", length(y_ou)),
             rep("Posterior predictive", length(unlist(reps)))))
ggplot2::ggplot(dens_df, ggplot2::aes(x = .data$value, fill = .data$type)) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::scale_fill_manual(values = c("Observed" = "grey30",
        "Posterior predictive" = "#2196F3")) +
    ggplot2::labs(x = "Value", y = "Density",
        title = "PPC: marginal density") + theme_relab()
```

``` r

dy_obs <- diff(y_ou)
dy_reps <- unlist(lapply(reps, diff))
inc_df <- data.frame(
    increment = c(dy_obs, dy_reps),
    type = c(rep("Observed", length(dy_obs)),
             rep("Posterior predictive", length(dy_reps))))
ggplot2::ggplot(inc_df, ggplot2::aes(x = .data$increment, fill = .data$type)) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::scale_fill_manual(values = c("Observed" = "grey30",
        "Posterior predictive" = "#2196F3")) +
    ggplot2::labs(x = expression(Delta*y), y = "Density",
        title = "PPC: increment distribution") + theme_relab()
```

``` r

running_var <- function(x, w = 10) {
    n <- length(x); v <- numeric(n - w + 1)
    for (i in seq_along(v)) v[i] <- var(x[i:(i + w - 1)]); v
}
rv_obs <- running_var(y_ou)
rv_q <- apply(do.call(rbind, lapply(reps, running_var)), 2,
              quantile, probs = c(0.025, 0.5, 0.975))
rv_df <- data.frame(time = times[10:length(times)], observed = rv_obs,
    median = rv_q[2, ], lo = rv_q[1, ], hi = rv_q[3, ])
ggplot2::ggplot(rv_df) +
    ggplot2::geom_ribbon(ggplot2::aes(x = .data$time, ymin = .data$lo,
        ymax = .data$hi), fill = "#2196F3", alpha = 0.3) +
    ggplot2::geom_line(ggplot2::aes(x = .data$time, y = .data$median),
        color = "#1565C0") +
    ggplot2::geom_line(ggplot2::aes(x = .data$time, y = .data$observed),
        color = "black") +
    ggplot2::labs(x = "Time", y = "Running variance (window = 10)",
        title = "PPC: running variance") + theme_relab()
```

These four PPC plots assess different aspects of the model fit. The ACF
check verifies temporal correlation structure. The marginal density
check ensures the overall distribution is correct. The increment
distribution is sensitive to diffusion misspecification. The running
variance check detects non-stationarity or heteroscedasticity.

#### Comprehensive posterior predictive checking

This section demonstrates a systematic PPC workflow for time series,
going beyond the basic checks above. The approach compares multiple
summary statistics of replicated datasets against the observed data to
assess calibration, dynamics, and distributional assumptions.

``` r

# Use the OU fit from above
set.seed(42)
posterior_ou <- fit_ou$Posterior2
if (is.null(posterior_ou)) posterior_ou <- fit_ou$Posterior1
n_rep <- 100
idx_rep <- sample.int(nrow(posterior_ou), n_rep, replace = TRUE)
y_reps <- list()
for (r in seq_len(n_rep)) {
    sim <- simulate(sde_ou, theta = posterior_ou[idx_rep[r], ], n.paths = 1)
    y_reps[[r]] <- sim$trajectories[, 1, 1]
}
```

#### Distribution of extremes

Extreme values are sensitive to tail behavior and noise amplitude. If
the model underestimates \\\sigma\\, the replicated maxima will be
systematically lower than the observed maximum.

``` r

obs_max <- max(y_ou)
obs_min <- min(y_ou)
rep_max <- sapply(y_reps, max)
rep_min <- sapply(y_reps, min)

ext_df <- data.frame(
    statistic = c(rep_max, rep_min),
    type = rep(c("Maximum", "Minimum"), each = n_rep)
)
ggplot2::ggplot(ext_df, ggplot2::aes(x = .data$statistic)) +
    ggplot2::geom_histogram(bins = 30, fill = "#2196F3", alpha = 0.6,
        color = "white") +
    ggplot2::geom_vline(data = data.frame(
        type = c("Maximum", "Minimum"),
        obs = c(obs_max, obs_min)),
        ggplot2::aes(xintercept = .data$obs),
        color = "black", linetype = "dashed", linewidth = 0.8) +
    ggplot2::facet_wrap(~type, scales = "free_x") +
    ggplot2::labs(x = "Value", y = "Count",
        title = "PPC: distribution of extremes",
        subtitle = "Dashed line: observed statistic") +
    theme_relab()
```

#### Lag-1 scatterplot

The lag-1 scatterplot reveals the autoregressive structure. For the OU
process, \\X\_{k+1} \mid X_k\\ is linear, so the lag-1 cloud should be
elliptical. Deviations from ellipticity suggest nonlinear dynamics that
the OU model misses.

``` r

# Observed lag-1
lag_obs <- data.frame(x = y_ou[-length(y_ou)], y = y_ou[-1],
    type = "Observed")
# One replicate for comparison
lag_rep <- data.frame(x = y_reps[[1]][-length(y_reps[[1]])],
    y = y_reps[[1]][-1], type = "Replicate")
lag_df <- rbind(lag_obs, lag_rep)
ggplot2::ggplot(lag_df, ggplot2::aes(x = .data$x, y = .data$y,
    color = .data$type)) +
    ggplot2::geom_point(alpha = 0.4, size = 1) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dotted",
        color = "grey50") +
    ggplot2::scale_color_manual(values = c("Observed" = "black",
        "Replicate" = "#2196F3")) +
    ggplot2::labs(x = expression(y[k]), y = expression(y[k+1]),
        title = "PPC: lag-1 scatterplot", color = NULL) +
    theme_relab()
```

#### Spectral density

The power spectral density reveals periodic structure that the model
should reproduce. For the OU process, the theoretical spectrum is
Lorentzian: \\S(f) = \frac{\sigma^2}{2\pi(\kappa^2 + 4\pi^2 f^2)}\\.

``` r

# Observed spectrum
spec_obs <- spectrum(y_ou, plot = FALSE)
spec_reps <- lapply(y_reps[1:30], function(r)
    spectrum(r, plot = FALSE)$spec[, 1])
spec_mat <- do.call(cbind, spec_reps)
spec_q <- apply(spec_mat, 1, quantile, probs = c(0.025, 0.5, 0.975))

spec_df <- data.frame(
    freq = spec_obs$freq,
    obs = log10(spec_obs$spec[, 1]),
    median = log10(spec_q[2, ]),
    lo = log10(spec_q[1, ]),
    hi = log10(spec_q[3, ])
)
ggplot2::ggplot(spec_df) +
    ggplot2::geom_ribbon(ggplot2::aes(x = .data$freq, ymin = .data$lo,
        ymax = .data$hi), fill = "#2196F3", alpha = 0.3) +
    ggplot2::geom_line(ggplot2::aes(x = .data$freq, y = .data$median),
        color = "#1565C0") +
    ggplot2::geom_line(ggplot2::aes(x = .data$freq, y = .data$obs),
        color = "black", linewidth = 0.5) +
    ggplot2::labs(x = "Frequency", y = expression(log[10](S(f))),
        title = "PPC: power spectral density") +
    theme_relab()
```

#### Bayesian p-values

The Bayesian p-value for a test statistic \\T\\ is \\\Pr(T(y\_{rep})
\geq T(y\_{obs}) \mid y\_{obs})\\. Values near 0.5 indicate good
calibration; values near 0 or 1 indicate systematic discrepancy.

``` r

# Compute Bayesian p-values for multiple statistics
T_obs <- c(
    mean = mean(y_ou),
    sd = sd(y_ou),
    skew = moments::skewness(y_ou),
    kurt = moments::kurtosis(y_ou),
    acf1 = acf(y_ou, lag.max = 1, plot = FALSE)$acf[2, , 1],
    range = diff(range(y_ou))
)
T_rep <- sapply(y_reps, function(r) {
    c(mean = mean(r), sd = sd(r),
      skew = moments::skewness(r),
      kurt = moments::kurtosis(r),
      acf1 = acf(r, lag.max = 1, plot = FALSE)$acf[2, , 1],
      range = diff(range(r)))
})
p_values <- rowMeans(T_rep >= T_obs)
knitr::kable(data.frame(
    Statistic = names(T_obs),
    Observed = round(T_obs, 4),
    `Bayesian p-value` = round(p_values, 3),
    Calibrated = ifelse(p_values > 0.05 & p_values < 0.95,
        "Yes", "No"),
    check.names = FALSE
), caption = "Bayesian p-values for summary statistics")
```

#### Cumulative distribution function (probability integral transform)

The probability integral transform (PIT) checks whether the observed
data are consistent with the predictive distribution. For a
well-calibrated model, the PIT values \\u_k = F(y_k \mid y\_{1:k-1},
\hat{\theta})\\ should be uniformly distributed on \\\[0, 1\]\\.

``` r

# Compute PIT values using posterior mean parameters
theta_hat <- colMeans(if (!is.null(fit_ou$Posterior2)) fit_ou$Posterior2 else fit_ou$Posterior1)
res_hat <- sde_ou$Model(theta_hat, sde_ou$Data)
pred_mean <- res_hat$yhat
obs_sd_hat <- theta_hat[4]
pit_values <- pnorm(y_ou, mean = pred_mean, sd = obs_sd_hat)

pit_df <- data.frame(u = pit_values)
ggplot2::ggplot(pit_df, ggplot2::aes(x = .data$u)) +
    ggplot2::geom_histogram(ggplot2::aes(
        y = ggplot2::after_stat(density)),
        bins = 20, fill = "#2196F3", alpha = 0.6, color = "white") +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed",
        color = "#FF5722") +
    ggplot2::labs(x = "PIT value u", y = "Density",
        title = "PPC: probability integral transform (should be uniform)",
        subtitle = "Dashed line: uniform reference") +
    theme_relab()
```

#### Cross-correlation structure (multivariate)

For multivariate SDEs, the cross-correlation between state components
carries information about the interaction structure. This check is
especially relevant for community ecology models.

``` r

# Applicable to the 3-species community model
# Compare observed cross-correlations with posterior predictive
if (exists("fit_eco")) {
    obs_cor <- cor(Y_eco)
    n_rep_eco <- 50
    post_eco <- fit_eco$Posterior2
    if (is.null(post_eco)) post_eco <- fit_eco$Posterior1
    idx_eco <- sample.int(nrow(post_eco), n_rep_eco, replace = TRUE)
    rep_cors <- array(NA, dim = c(N_sp, N_sp, n_rep_eco))
    for (r in seq_len(n_rep_eco)) {
        sim_eco <- simulate(sde_eco,
            theta = post_eco[idx_eco[r], ], n.paths = 1)
        rep_cors[, , r] <- cor(sim_eco$trajectories[, , 1])
    }
    cor_names <- c("Prey-Meso", "Prey-Apex", "Meso-Apex")
    cor_obs <- c(obs_cor[1,2], obs_cor[1,3], obs_cor[2,3])
    cor_reps <- rbind(
        rep_cors[1, 2, ], rep_cors[1, 3, ], rep_cors[2, 3, ])
    cor_df <- data.frame(
        correlation = as.vector(t(cor_reps)),
        pair = rep(cor_names, each = n_rep_eco)
    )
    obs_cor_df <- data.frame(pair = cor_names, obs = cor_obs)

    ggplot2::ggplot(cor_df, ggplot2::aes(x = .data$correlation)) +
        ggplot2::geom_histogram(bins = 20, fill = "#2196F3",
            alpha = 0.6, color = "white") +
        ggplot2::geom_vline(data = obs_cor_df,
            ggplot2::aes(xintercept = .data$obs),
            color = "black", linetype = "dashed", linewidth = 0.8) +
        ggplot2::facet_wrap(~pair, scales = "free_x") +
        ggplot2::labs(x = "Cross-correlation", y = "Count",
            title = "PPC: interspecific cross-correlations",
            subtitle = "Dashed: observed; histogram: posterior predictive") +
        theme_relab()
}
```

### Model comparison and selection

#### DIC-based model comparison

A central question in SDE modelling is which process best describes the
data. lucifer supports model comparison via DIC (computed automatically
during fitting) and LOO-PSIS (via the
[`LOO()`](https://robustecologies.github.io/lucifer/reference/LOO.md)
function). The following example compares three candidate models for
mean-reverting data: the OU process, the CIR process, and the
double-well potential.

``` r

# Use the OU data from above
sde_m1 <- SDE(data = y_ou, times = times, family = "ou")
sde_m2 <- SDE(data = pmax(y_ou, 0.01), times = times, family = "cir")
sde_m3 <- SDE(data = y_ou, times = times, family = "double_well")

fit_m1 <- SDE.fit(sde_m1, Algorithm = "HARM",
                   Iterations = 20000, Thinning = 5)
fit_m2 <- SDE.fit(sde_m2, Algorithm = "HARM",
                   Iterations = 20000, Thinning = 5)
fit_m3 <- SDE.fit(sde_m3, Algorithm = "HARM",
                   Iterations = 20000, Thinning = 5)
```

``` r

# Extract DIC from each fit
dic_table <- data.frame(
    Model = c("OU (exact)", "CIR (exact)", "Double-well (EKF)"),
    DIC = c(
        extract_dic(fit_m1)$DIC,
        extract_dic(fit_m2)$DIC,
        extract_dic(fit_m3)$DIC
    )
)
dic_table$Delta_DIC <- dic_table$DIC - min(dic_table$DIC)
knitr::kable(dic_table, digits = 1,
             caption = "DIC comparison of three SDE models")
```

Since the data were generated from an OU process, we expect the OU model
to have the lowest DIC. The CIR model is a close competitor because it
reduces to the OU when the state is far from zero, but the non-central
chi-squared likelihood adds unnecessary complexity. The double-well
model penalizes because the cubic drift does not match the linear
data-generating process.

#### Posterior predictive comparison

A visual comparison of posterior predictive trajectories reveals how
each model captures the dynamics:

``` r

par(mfrow = c(3, 1))
plot(fit_m1, type = "predictive", n.paths = 20)
title(main = "OU process", line = 0.5)
plot(fit_m2, type = "predictive", n.paths = 20)
title(main = "CIR process", line = 0.5)
plot(fit_m3, type = "predictive", n.paths = 20)
title(main = "Double-well potential", line = 0.5)
```

#### Comprehensive model selection workflow

This section demonstrates a systematic workflow for selecting between
competing SDE models using DIC, posterior predictive loss, and visual
diagnostics. The workflow applies the principle that no single criterion
is sufficient; multiple lines of evidence should converge.

``` r

# Four candidate models for the OU data:
# M1: OU (correct)
# M2: GBM (wrong: non-mean-reverting)
# M3: Double-well (wrong: cubic nonlinearity)
# M4: CIR (wrong: sqrt-diffusion)
models <- list(
    OU = SDE(data = y_ou, times = times, family = "ou"),
    GBM = SDE(data = abs(y_ou), times = times, family = "gbm"),
    `Double-well` = SDE(data = y_ou, times = times,
        family = "double_well"),
    CIR = SDE(data = abs(y_ou) + 0.1, times = times, family = "cir")
)

fits <- lapply(models, function(m)
    SDE.fit(m, Algorithm = "HARM", Iterations = 15000, Thinning = 5))
```

``` r

dic_results <- data.frame(
    Model = names(fits),
    DIC = sapply(fits, function(f) round(extract_dic(f)$DIC, 1)),
    pD = sapply(fits, function(f) round(extract_dic(f)$pD, 1)),
    `Mean deviance` = sapply(fits, function(f) round(extract_dic(f)$Dbar, 1)),
    check.names = FALSE
)
dic_results$Delta_DIC <- dic_results$DIC - min(dic_results$DIC)
knitr::kable(dic_results,
    caption = "Model selection via DIC")
```

``` r

# PPC comparison: posterior predictive density overlay for each model
ppc_list <- list()
for (nm in names(fits)) {
    post <- fits[[nm]]$Posterior2
    if (is.null(post)) post <- fits[[nm]]$Posterior1
    n_r <- min(50, nrow(post))
    idx_r <- sample.int(nrow(post), n_r)
    reps <- numeric(0)
    for (r in seq_len(n_r)) {
        sim <- simulate(models[[nm]], theta = post[idx_r[r], ],
            n.paths = 1)
        reps <- c(reps, sim$trajectories[, 1, 1])
    }
    ppc_list[[nm]] <- data.frame(
        value = c(y_ou, reps),
        type = c(rep("Observed", length(y_ou)),
                 rep("Predictive", length(reps))),
        model = nm
    )
}
ppc_all <- do.call(rbind, ppc_list)

ggplot2::ggplot(ppc_all, ggplot2::aes(x = .data$value,
    fill = .data$type)) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::scale_fill_manual(values = c("Observed" = "grey30",
        "Predictive" = "#2196F3")) +
    ggplot2::facet_wrap(~model, scales = "free_y", ncol = 2) +
    ggplot2::labs(x = "Value", y = "Density",
        title = "PPC: marginal predictive density by model",
        fill = NULL) +
    theme_relab()
```

``` r

# Residual comparison across models
res_list <- list()
for (nm in names(fits)) {
    post_nm <- if (!is.null(fits[[nm]]$Posterior2)) fits[[nm]]$Posterior2 else fits[[nm]]$Posterior1
    post_mean <- colMeans(post_nm)
    res_out <- models[[nm]]$Model(post_mean, models[[nm]]$Data)
    yhat_m <- res_out$yhat
    obs_m <- if (nm %in% c("GBM", "CIR")) abs(y_ou) + 0.1 else y_ou
    residuals_m <- obs_m[seq_along(yhat_m)] - yhat_m
    res_list[[nm]] <- data.frame(
        time = times[seq_along(residuals_m)],
        residual = residuals_m,
        model = nm
    )
}
res_all <- do.call(rbind, res_list)

ggplot2::ggplot(res_all, ggplot2::aes(x = .data$time,
    y = .data$residual)) +
    ggplot2::geom_point(alpha = 0.3, size = 0.8, color = "#2196F3") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
        color = "grey50") +
    ggplot2::geom_smooth(method = "loess", se = FALSE,
        color = "#F44336", linewidth = 0.6) +
    ggplot2::facet_wrap(~model, ncol = 2) +
    ggplot2::labs(x = "Time", y = "Residual",
        title = "One-step-ahead residuals by model",
        subtitle = "Red: LOESS smooth (should be flat at zero)") +
    theme_relab()
```

The correct model (OU) should show: lowest DIC, predictive density
matching the observed distribution, and flat residuals with constant
variance. Misspecified models typically show systematic residual
patterns: GBM residuals trend because it cannot mean-revert; double-well
residuals oscillate because the cubic drift overreacts; CIR residuals
show heteroscedasticity because the sqrt-diffusion mismatches the
constant noise.

#### Observation model comparison

The same latent dynamics can be observed through different observation
models. This example compares fitting an OU process with observation
noise (correct) versus without (misspecified):

``` r

sde_correct <- SDE(data = y_ou, times = times, family = "ou")
sde_misspec <- SDE(
    drift = function(x, theta, t) theta[1] * (theta[2] - x),
    diffusion = function(x, theta, t) theta[3],
    data = y_ou, times = times, x0 = y_ou[1],
    obs.model = "gaussian",
    obs.noise = function(theta) 1e-10,
    prior = function(theta) {
        if (theta[1] <= 0 || theta[3] <= 0) return(-Inf)
        sum(dnorm(theta, 0, 10, log = TRUE))
    },
    parm.names = c("kappa", "mu", "sigma"),
    method = "exact", family = "ou")

fit_c <- SDE.fit(sde_correct, Algorithm = "HARM", Iterations = 15000, Thinning = 3)
fit_m <- SDE.fit(sde_misspec, Algorithm = "HARM", Iterations = 15000, Thinning = 3)

dic_c <- extract_dic(fit_c)$DIC
dic_m <- extract_dic(fit_m)$DIC
knitr::kable(data.frame(
    Model = c("With observation noise (correct)", "Without observation noise (misspecified)"),
    DIC = round(c(dic_c, dic_m), 1),
    Delta_DIC = round(c(dic_c, dic_m) - min(dic_c, dic_m), 1)
), caption = "Observation model comparison via DIC")
```

Omitting observation noise forces \\\sigma\\ to absorb the measurement
variability, producing biased parameter estimates and a worse DIC.

### Detecting model misspecification

This example demonstrates how posterior predictive checks systematically
reveal specific types of model misspecification. We fit two wrong models
to OU data and show how different PPC statistics diagnose each failure
mode.

``` r

# OU data: the correct model
# Misspecification 1: GBM (no mean reversion, exponential growth)
# Misspecification 2: Double-well (cubic nonlinearity)

sde_right <- SDE(data = y_ou, times = times, family = "ou")
sde_wrong1 <- SDE(data = abs(y_ou) + 0.1, times = times, family = "gbm")
sde_wrong2 <- SDE(data = y_ou, times = times, family = "double_well")

fit_right <- SDE.fit(sde_right, Algorithm = "HARM",
    Iterations = 15000, Thinning = 5)
fit_wrong1 <- SDE.fit(sde_wrong1, Algorithm = "HARM",
    Iterations = 15000, Thinning = 5)
fit_wrong2 <- SDE.fit(sde_wrong2, Algorithm = "HARM",
    Iterations = 15000, Thinning = 5)
```

``` r

# Side-by-side trajectory comparison
par_list <- list(
    list(fit = fit_right, name = "OU (correct)"),
    list(fit = fit_wrong1, name = "GBM (no mean-reversion)"),
    list(fit = fit_wrong2, name = "Double-well (cubic drift)")
)

traj_dfs <- lapply(par_list, function(p) {
    post <- p$fit$Posterior2
    if (is.null(post)) post <- p$fit$Posterior1
    theta_hat <- colMeans(post)
    yhat <- p$fit$sde$Model(theta_hat, p$fit$sde$Data)$yhat
    data.frame(
        time = times[seq_along(yhat)],
        observed = y_ou[seq_along(yhat)],
        fitted = yhat,
        model = p$name
    )
})
traj_all <- do.call(rbind, traj_dfs)

ggplot2::ggplot(traj_all) +
    ggplot2::geom_point(ggplot2::aes(x = .data$time,
        y = .data$observed), color = "grey40", size = 0.5, alpha = 0.5) +
    ggplot2::geom_line(ggplot2::aes(x = .data$time,
        y = .data$fitted), color = "#2196F3", linewidth = 0.6) +
    ggplot2::facet_wrap(~model, ncol = 1) +
    ggplot2::labs(x = "Time", y = "Value",
        title = "Filtered trajectories: correct vs. misspecified models",
        subtitle = "Grey points: observations; blue line: posterior mean") +
    theme_relab()
```

``` r

# PPC comparison: autocorrelation structure
models_list <- list(
    "OU (correct)" = list(sde = sde_right, fit = fit_right),
    "GBM (wrong)" = list(sde = sde_wrong1, fit = fit_wrong1),
    "Double-well (wrong)" = list(sde = sde_wrong2, fit = fit_wrong2)
)

acf_comp_dfs <- list()
for (nm in names(models_list)) {
    post <- models_list[[nm]]$fit$Posterior2
    if (is.null(post)) post <- models_list[[nm]]$fit$Posterior1
    n_r <- min(30, nrow(post))
    idx_r <- sample.int(nrow(post), n_r)
    acf_reps <- matrix(NA, n_r, 21)
    for (r in seq_len(n_r)) {
        sim <- simulate(models_list[[nm]]$sde,
            theta = post[idx_r[r], ], n.paths = 1)
        acf_reps[r, ] <- acf(sim$trajectories[, 1, 1],
            lag.max = 20, plot = FALSE)$acf[, , 1]
    }
    acf_q <- apply(acf_reps, 2, quantile, probs = c(0.025, 0.5, 0.975))
    acf_obs <- acf(y_ou, lag.max = 20, plot = FALSE)$acf[, , 1]
    acf_comp_dfs[[nm]] <- data.frame(
        lag = 0:20,
        obs = acf_obs,
        median = acf_q[2, ],
        lo = acf_q[1, ],
        hi = acf_q[3, ],
        model = nm
    )
}
acf_comp <- do.call(rbind, acf_comp_dfs)

ggplot2::ggplot(acf_comp) +
    ggplot2::geom_ribbon(ggplot2::aes(x = .data$lag, ymin = .data$lo,
        ymax = .data$hi), fill = "#2196F3", alpha = 0.3) +
    ggplot2::geom_line(ggplot2::aes(x = .data$lag, y = .data$median),
        color = "#1565C0") +
    ggplot2::geom_line(ggplot2::aes(x = .data$lag, y = .data$obs),
        color = "black", linewidth = 0.6) +
    ggplot2::facet_wrap(~model, ncol = 1) +
    ggplot2::labs(x = "Lag", y = "ACF",
        title = "PPC: autocorrelation under correct and misspecified models",
        subtitle = "Black: observed; blue: posterior predictive 95% band") +
    theme_relab()
```

The diagnostic pattern is revealing. GBM fails the ACF check because it
cannot mean-revert: the simulated ACF decays too slowly (or not at all)
compared to the observed exponential decay \\\rho(\tau) =
e^{-\kappa\tau}\\. The double-well model may pass the ACF check at short
lags (the cubic drift produces local mean-reversion around each well)
but generates occasional large jumps between wells that inflate the
variance and produce excess kurtosis not present in the data.

### Simulation-based calibration

Simulation-based calibration (SBC) verifies that the inference procedure
is well-calibrated: for data simulated from the prior predictive, the
posterior rank statistics should be uniformly distributed. This is a
stringent self-consistency check that detects bugs in the likelihood,
incorrect constraint handling, or poor MCMC mixing [\[32\]](#ref32).

``` r

# SBC for the OU model with a small number of replications
set.seed(42)
n_sbc <- 20  # increase to 200+ for production
n_obs_sbc <- 50
times_sbc <- seq(0, 5, by = 0.1)
ranks <- matrix(NA, n_sbc, 4)
colnames(ranks) <- c("kappa", "mu", "sigma", "obs_sd")

for (rep in seq_len(n_sbc)) {
    # Draw parameters from prior
    kappa_p <- abs(rnorm(1, 0, 3))
    mu_p <- rnorm(1, 0, 5)
    sigma_p <- abs(rnorm(1, 0, 2))
    obs_sd_p <- abs(rnorm(1, 0, 1))
    theta_p <- c(max(kappa_p, 0.1), mu_p,
                 max(sigma_p, 0.01), max(obs_sd_p, 0.01))

    # Simulate data from prior predictive
    x_sbc <- numeric(length(times_sbc))
    x_sbc[1] <- mu_p
    for (k in 2:length(times_sbc)) {
        dt_s <- times_sbc[k] - times_sbc[k - 1]
        ekt <- exp(-theta_p[1] * dt_s)
        x_sbc[k] <- rnorm(1,
            theta_p[2] + (x_sbc[k-1] - theta_p[2]) * ekt,
            theta_p[3] * sqrt((1 - exp(-2 * theta_p[1] * dt_s)) /
                (2 * theta_p[1])))
    }
    y_sbc <- x_sbc + rnorm(length(times_sbc), 0, theta_p[4])

    # Fit
    sde_sbc <- SDE(data = y_sbc, times = times_sbc, family = "ou")
    fit_sbc <- SDE.fit(sde_sbc, Algorithm = "RWM",
                        Iterations = 2000, Thinning = 2)
    post_sbc <- fit_sbc$Posterior1

    # Compute ranks
    for (j in 1:4)
        ranks[rep, j] <- mean(post_sbc[, j] < theta_p[j])
}
```

``` r

# SBC rank histogram: should be uniform
rank_df <- data.frame(
    rank = as.vector(ranks),
    parameter = rep(colnames(ranks), each = n_sbc)
)
ggplot2::ggplot(rank_df, ggplot2::aes(x = .data$rank)) +
    ggplot2::geom_histogram(bins = 10, fill = "#2196F3",
        alpha = 0.6, color = "white") +
    ggplot2::geom_hline(yintercept = n_sbc / 10, linetype = "dashed",
        color = "#FF5722") +
    ggplot2::facet_wrap(~parameter, scales = "free_y") +
    ggplot2::labs(x = "Rank (quantile of true value in posterior)",
        y = "Count",
        title = "Simulation-based calibration",
        subtitle = paste0(n_sbc, " replications; dashed: uniform reference")) +
    theme_relab()
```

Uniform rank histograms indicate correct calibration. A U-shaped
histogram indicates the credible intervals are too narrow
(overconfident); an inverted-U indicates they are too wide
(underconfident); a skewed histogram indicates systematic bias.

## Advanced analysis

### Robust Bayesian sensitivity analysis

lucifer’s
[`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md)
function performs power-scaling sensitivity analysis [\[18\]](#ref18) to
assess how sensitive posterior inferences are to the prior and
likelihood. For SDE models, this is useful for diagnosing prior-data
conflict (e.g., an informative prior on \\\kappa\\ that contradicts the
data) or identifying parameters whose posteriors are dominated by the
prior rather than the data.

``` r

rb_ou <- RobustBayes(fit_m1)
summary(rb_ou)
plot(rb_ou)
```

The sensitivity plot shows how each parameter’s posterior shifts as the
prior or likelihood is power-scaled. Parameters with high prior
sensitivity and low likelihood sensitivity are poorly identified by the
data. Parameters with high likelihood sensitivity but low prior
sensitivity are well-informed by the data. The
[`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md)
function works with `sde_fit` objects exactly as with any `demonoid`
fit, since `sde_fit` inherits from `demonoid`.

lucifer’s
[`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md)
[\[18\]](#ref18) performs power-scaling sensitivity analysis on SDE
fits, assessing how posterior inferences depend on the prior and
likelihood:

``` r

rb <- RobustBayes(fit_ou)
summary(rb)
plot(rb)
```

The sensitivity plot reveals which parameters are prior-dominated (high
prior sensitivity, low likelihood sensitivity) versus data-informed (low
prior sensitivity, high likelihood sensitivity). For well-identified SDE
models, the observation noise \\\sigma\_{obs}\\ is typically
well-informed by the data, while the reversion speed \\\kappa\\ may be
more sensitive to the prior if the time series is short relative to the
mean-reversion timescale \\1/\kappa\\.

Power-scaling sensitivity analysis [\[18\]](#ref18) perturbs the prior
and likelihood by raising them to a power \\\alpha\\, then examines how
the posterior shifts. This diagnoses: (1) prior-data conflict (the prior
and likelihood pull in different directions), (2) prior sensitivity (the
posterior changes substantially when the prior is weakened), and (3)
likelihood sensitivity (the posterior changes when the likelihood is
weakened, indicating fragile inference).

``` r

# Full RobustBayes analysis on OU fit
rb <- RobustBayes(fit_ou)
summary(rb)
```

``` r

# Power-scaling sensitivity plot
plot(rb)
```

``` r

# Compare prior vs. likelihood sensitivity
rb_df <- data.frame(
    Parameter = rb$parm.names,
    Prior_sensitivity = round(rb$prior_sensitivity, 4),
    Likelihood_sensitivity = round(rb$likelihood_sensitivity, 4),
    Diagnosis = ifelse(rb$prior_sensitivity > 0.05, "Prior-sensitive",
        ifelse(rb$likelihood_sensitivity > 0.05, "Likelihood-sensitive",
            "Well-identified"))
)
knitr::kable(rb_df,
    caption = "RobustBayes parameter sensitivity classification")
```

Parameters classified as “prior-sensitive” have posteriors dominated by
the prior rather than the data. This typically occurs for parameters
that are weakly identified by the likelihood (e.g., the observation
noise \\\sigma\_{obs}\\ when it is small relative to the process noise
\\\sigma\\). Parameters classified as “well-identified” have posteriors
driven by the data regardless of prior perturbation.

### Frequentist-Bayesian convergence

For well-identified models with large samples, Bayesian posterior
summaries converge to frequentist maximum-likelihood estimates (the
Bernstein-von Mises theorem). lucifer’s frequentist bridge functions
allow direct comparison. This example fits an OU process with both
[`SDE.fit()`](https://robustecologies.github.io/lucifer/reference/SDE.fit.md)
(Bayesian) and
[`LaplaceApproximation()`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md)
(frequentist-like), then compares via
[`confint_compare()`](https://robustecologies.github.io/lucifer/reference/confint_compare.md).

``` r

# Large OU dataset for good asymptotics
set.seed(42)
true_k <- 3; true_m <- 5; true_s <- 0.8; true_o <- 0.15
times_fb <- seq(0, 20, by = 0.05)
n_fb <- length(times_fb)
x_fb <- numeric(n_fb); x_fb[1] <- true_m
for (k in 2:n_fb) {
    dt_fb <- times_fb[k] - times_fb[k - 1]
    ekt <- exp(-true_k * dt_fb)
    x_fb[k] <- rnorm(1, true_m + (x_fb[k - 1] - true_m) * ekt,
        true_s * sqrt((1 - exp(-2 * true_k * dt_fb)) / (2 * true_k)))
}
y_fb <- x_fb + rnorm(n_fb, 0, true_o)
```

``` r

# Bayesian fit via MCMC
sde_fb <- SDE(data = y_fb, times = times_fb, family = "ou")
fit_bayes <- SDE.fit(sde_fb, Algorithm = "HARM",
                      Iterations = 30000, Thinning = 10)

# Frequentist-like fit via Laplace approximation
# LaplaceApproximation requires a clean Data list (no internal SDE fields),
# so we build a minimal wrapper.
la_Model <- sde_fb$Model
la_Data <- list(
    N = length(y_fb),
    parm.names = sde_fb$parm.names,
    mon.names = "LP"
)
# Wrap Model to return scalar Monitor (LA expects single mon.names = "LP")
la_Model_wrap <- function(parm, Data) {
    res <- sde_fb$Model(parm, sde_fb$Data)
    res$Monitor <- res$LP
    res
}
fit_laplace <- LaplaceApproximation(
    la_Model_wrap, sde_fb$Initial.Values, la_Data,
    Method = "NM", Iterations = 5000)
```

``` r

# Compare confidence/credible intervals
# confint_compare expects demonoid class; strip sde_fit wrapper
fit_bayes_dem <- fit_bayes
class(fit_bayes_dem) <- "demonoid"
confint_compare(
    list("MCMC (Bayesian)" = fit_bayes_dem,
         "Laplace (frequentist)" = fit_laplace),
    true_values = c(kappa = true_k, mu = true_m,
                    sigma = true_s, obs_sd = true_o)
)
```

The
[`confint_compare()`](https://robustecologies.github.io/lucifer/reference/confint_compare.md)
forest plot shows point estimates with 95% intervals for each method.
For a well-identified model with 400 observations, the Bayesian credible
intervals and frequentist confidence intervals should nearly overlap,
confirming asymptotic equivalence. Discrepancies indicate either prior
influence (small sample) or multimodality (poor mixing).

``` r

# Profile likelihood: frequentist CI shape
prof <- profile_likelihood(fit_laplace, la_Model_wrap, la_Data,
                            parm = c("kappa", "mu"), n_grid = 40)
plot(prof)
```

The profile likelihood curves show the curvature of the log-likelihood
around each parameter’s MLE. Symmetric curves indicate well-behaved
quadratic approximation (Laplace approximation is accurate); asymmetric
curves suggest the posterior may be skewed and MCMC is preferable.

``` r

# Frequentist summary with p-values and CIs
fs <- freq_summary(fit_laplace)
knitr::kable(fs$coefficients, digits = 4,
    caption = "Frequentist summary of OU parameter estimates")
```

### LOO-PSIS for particle filter models

For SDE models fitted with the particle filter,
[`log_lik()`](https://robustecologies.github.io/lucifer/reference/log_lik.md)
extracts per-observation log-likelihood increments by re-running the PF
for each posterior sample, enabling LOO-PSIS [\[25\]](#ref25):

``` r

set.seed(42)
times_short <- seq(0, 5, by = 0.5)
sde_lv_s <- SDE(data = cbind(rpois(length(times_short), 50),
    rpois(length(times_short), 30)),
    times = times_short, family = "lotka_volterra")
fit_lv_s <- SDE.fit.pmmh(sde_lv_s, N.particles = 100, Iterations = 1000)

ll_mat <- log_lik(fit_lv_s)
cat("log_lik matrix:", nrow(ll_mat), "x", ncol(ll_mat), "\n")
loo_result <- LOO(fit_lv_s)
summary(loo_result)
```

### State smoothing

The forward-filtering backward-sampling (FFBSm) smoother produces
smoothed state estimates \\p(X\_{0:T} \mid y\_{1:T}, \theta)\\ that
condition on all observations (past and future), unlike the filter which
only conditions on past observations [\[11\]](#ref11). This is useful
for reconstructing latent state trajectories and for assessing the
uncertainty in the hidden state at each time point.

``` r

# Smooth the Lotka-Volterra fit
smoothed <- smooth.sde_fit(fit_lv, N.particles = 500,
                            n.trajectories = 50)

# Compare filtered vs smoothed means
par(mfrow = c(2, 1))
matplot(times_lv, smoothed$filtered_mean, type = "l", lty = 1,
        col = c("blue", "red"), xlab = "Time", ylab = "State",
        main = "Filtered state estimates")
legend("topright", c("Prey (filtered)", "Predator (filtered)"),
       col = c("blue", "red"), lty = 1)
matplot(times_lv, smoothed$smoothed_mean, type = "l", lty = 1,
        col = c("blue", "red"), xlab = "Time", ylab = "State",
        main = "Smoothed state estimates")
legend("topright", c("Prey (smoothed)", "Predator (smoothed)"),
       col = c("blue", "red"), lty = 1)
```

The smoothed estimates are generally smoother than the filtered
estimates and have narrower uncertainty bands, particularly at the
beginning of the time series where the filter has not yet accumulated
much information.

### Forecasting and prediction

The [`predict()`](https://rdrr.io/r/stats/predict.html) method for
`sde_fit` objects produces posterior predictive forecasts that propagate
both parameter uncertainty and process noise:

``` r

# Forecast the OU process
pred_ou <- predict(fit_ou, n.ahead = 50, n.paths = 500,
                    probs = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975))
print(pred_ou)
plot(pred_ou)
```

For the OU process, the prediction intervals converge to the stationary
distribution \\N(\mu, \sigma^2/(2\kappa))\\ as the horizon grows. For
non-stationary processes like GBM, the intervals grow without bound.

### Prior predictive analysis and simulation

Before fitting, simulate from the prior predictive to verify the prior
produces plausible dynamics:

``` r

sde_mod <- SDE(data = y_ou, times = times, family = "ou")
sim <- simulate(sde_mod, n.paths = 5, seed = 42)
cat("Trajectory array:", paste(dim(sim$trajectories), collapse = " x "), "\n")
plot(sde_mod, type = "prior_predictive", n.paths = 20)
```

## Practical workflow summary

A complete SDE analysis in lucifer follows this sequence:

1.  **Specify the model.** Use `SDE(family = "...")` for built-in
    families or custom drift/diffusion functions. For Stratonovich
    models, add `interpretation = "stratonovich"`.

2.  **Prior predictive check.** Run
    `plot(sde, type = "prior_predictive")` to verify the prior produces
    plausible dynamics.

3.  **Fit.** Use
    [`SDE.fit()`](https://robustecologies.github.io/lucifer/reference/SDE.fit.md)
    with HARM, AMWG, or twalk for exact/EKF models. Use
    [`SDE.fit.pmmh()`](https://robustecologies.github.io/lucifer/reference/SDE.fit.pmmh.md)
    for particle filter models.

4.  **Convergence diagnostics.** Check `plot(fit, type = "posterior")`
    for trace plots, density, and ACF. Inspect running means and Geweke
    diagnostics for each parameter.

5.  **Posterior predictive checks.** Run at minimum: trajectory plot
    (`type = "trajectory"`), predictive paths (`type = "predictive"`),
    residuals (`type = "residuals"`), and the ACF/marginal/spectral PPC
    suite.

6.  **Model comparison.** Compare DIC across candidates. For particle
    filter models, use
    [`LOO()`](https://robustecologies.github.io/lucifer/reference/LOO.md)
    with
    [`log_lik()`](https://robustecologies.github.io/lucifer/reference/log_lik.md).
    Visualize PPC faceted by model.

7.  **Sensitivity analysis.** Run
    [`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md)
    to check prior-data conflict and parameter identifiability.

8.  **Frequentist bridge.** For publication, verify that Bayesian
    credible intervals converge to frequentist confidence intervals via
    [`confint_compare()`](https://robustecologies.github.io/lucifer/reference/confint_compare.md).

9.  **Forecasting.** Use
    [`predict()`](https://rdrr.io/r/stats/predict.html) with appropriate
    `n.ahead` and visualize prediction intervals.

10. **Report.** The fitted `sde_fit` object stores everything needed:
    posterior samples, model specification, diagnostics, and the SDE
    closure for reproducibility.

## API reference

### Core functions

| Function | Purpose |
|----|----|
| [`SDE()`](https://robustecologies.github.io/lucifer/reference/SDE.md) | Construct an `sde_model` from drift, diffusion, observations, and priors |
| [`SDE.fit()`](https://robustecologies.github.io/lucifer/reference/SDE.fit.md) | Fit an `sde_model` via MCMC, returns `sde_fit` |
| [`SDE.fit.pmmh()`](https://robustecologies.github.io/lucifer/reference/SDE.fit.pmmh.md) | Particle MCMC wrapper with recommended defaults |
| [`simulate()`](https://rdrr.io/r/stats/simulate.html) | Forward simulation from an `sde_model` |
| [`predict()`](https://rdrr.io/r/stats/predict.html) | Posterior predictive forecasting from an `sde_fit` |
| [`smooth.sde_fit()`](https://robustecologies.github.io/lucifer/reference/smooth.sde_fit.md) | Forward-filtering backward-sampling smoother |

### S3 methods for `sde_model`

| Method | Description |
|----|----|
| [`print()`](https://rdrr.io/r/base/print.html) | One-screen summary of model specification |
| [`summary()`](https://rdrr.io/r/base/summary.html) | Extended output: dimensions, constraints, data statistics |
| `plot(type = "data")` | Observed time series |
| `plot(type = "phase")` | Phase portrait for 2D+ systems |
| `plot(type = "prior_predictive")` | Prior predictive paths overlaid on data |

### S3 methods for `sde_fit`

| Method | Description |
|----|----|
| [`print()`](https://rdrr.io/r/base/print.html) | Posterior parameter summary |
| [`summary()`](https://rdrr.io/r/base/summary.html) | Extended MCMC and posterior diagnostics |
| `plot(type = "trajectory")` | Posterior mean + credible bands on data |
| `plot(type = "predictive")` | Posterior predictive paths |
| `plot(type = "phase")` | Phase portrait with posterior trajectories |
| `plot(type = "residuals")` | One-step-ahead residuals |
| `plot(type = "posterior")` | MCMC trace, density, ACF |
| `plot(type = "filtered")` | Filtered state + ESS (particle filter) |

### Compatibility with lucifer ecosystem

| Function | Works with `sde_fit`? | Notes |
|----|----|----|
| [`LOO()`](https://robustecologies.github.io/lucifer/reference/LOO.md) | Yes | Via [`log_lik.sde_fit()`](https://robustecologies.github.io/lucifer/reference/log_lik.md) for particle filter models |
| [`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md) | Yes | Power-scaling sensitivity analysis |
| [`Prescribe()`](https://robustecologies.github.io/lucifer/reference/Prescribe.md) | Yes | Detects SDE models, warns about gradient algorithms |
| [`Crucible()`](https://robustecologies.github.io/lucifer/reference/Crucible.md) | Yes | Algorithm benchmarking |
| [`predict()`](https://rdrr.io/r/stats/predict.html) | Yes | Posterior predictive forecasting |
| [`PosteriorChecks()`](https://robustecologies.github.io/lucifer/reference/PosteriorChecks.md) | Yes | Inherited from `demonoid` |

## Appendix: Ito vs. Stratonovich interpretation

lucifer’s SDE engines use Ito calculus internally, but models can be
specified in either Ito or Stratonovich form. The `interpretation`
parameter triggers automatic symbolic conversion when the Stratonovich
form is used.

### The difference

For SDEs with state-dependent (multiplicative) noise, the Ito and
Stratonovich stochastic integrals produce different transition
densities. A Stratonovich SDE \\dX = f(X)dt + g(X) \circ dW\\ is
equivalent to the Ito SDE

\\dX = \left\[f(X) + \frac{1}{2}g(X)\frac{\partial g}{\partial
X}(X)\right\]dt + g(X)\\dW\\

The correction \\\frac{1}{2}g\\g'\\ is the noise-induced drift
[\[1\]](#ref1). For additive noise (\\g\\ constant, \\g' = 0\\), the two
interpretations are identical.

### When Stratonovich matters

Stratonovich calculus is the natural choice when the SDE arises as the
limit of an ODE driven by smooth colored noise (the Wong-Zakai theorem
[\[28\]](#ref28)). This includes mechanical systems (Langevin
equations), neural models with synaptic filtering, and ecological models
where environmental noise has finite correlation time. Ito SDEs arise
naturally in finance (Black-Scholes), population genetics
(Wright-Fisher), and any context where the noise is a genuine
white-noise innovation.

### Automatic conversion in lucifer

When drift and diffusion are specified as expressions (character strings
or [`quote()`](https://rdrr.io/r/base/substitute.html)), lucifer’s
symbolic differentiation engine computes \\g'(x)\\ automatically and
adds the correction to the drift. The user specifies
`interpretation = "stratonovich"` and the rest is handled internally.

``` r

# Stratonovich geometric growth: dX = mu*X dt + sigma*X o dW
# The Ito equivalent adds a noise-induced drift: (mu + 0.5*sigma^2)*X
# This means the Stratonovich drift mu produces HIGHER effective growth
# than the Ito drift mu (by 0.5*sigma^2).

# Simulate data from the Ito version (true model)
set.seed(42)
true_mu <- 0.1; true_sigma <- 0.3; true_obs_sd <- 0.2
times_s <- seq(0, 10, by = 0.1)
n_s <- length(times_s)
x_s <- numeric(n_s); x_s[1] <- 2
for (k in 2:n_s) {
    dt_s <- times_s[k] - times_s[k - 1]
    # Ito dynamics (drift = mu*X, no correction)
    x_s[k] <- x_s[k - 1] + true_mu * x_s[k - 1] * dt_s +
        true_sigma * x_s[k - 1] * sqrt(dt_s) * rnorm(1)
    x_s[k] <- max(x_s[k], 0.01)
}
y_s <- abs(x_s + rnorm(n_s, 0, true_obs_sd))
```

``` r

# Fit as Stratonovich model: lucifer auto-corrects the drift
sde_strat <- SDE(
    drift = "theta[1] * x[1]",       # Stratonovich drift (no manual correction)
    diffusion = "theta[2] * x[1]",   # multiplicative noise
    data = y_s, times = times_s, x0 = y_s[1],
    obs.model = "gaussian",
    obs.noise = function(theta) theta[3],
    prior = function(theta) {
        if (theta[2] <= 0 || theta[3] <= 0) return(-Inf)
        dnorm(theta[1], 0, 1, log = TRUE) +
            dnorm(theta[2], 0, 1, log = TRUE) +
            dexp(theta[3], 5, log = TRUE)
    },
    parm.names = c("mu", "sigma", "obs_sd"),
    method = "ekf",
    interpretation = "stratonovich"   # auto-conversion to Ito
)

# Fit as Ito model (manual correction for comparison)
sde_ito <- SDE(
    drift = "theta[1] * x[1] + 0.5 * theta[2]^2 * x[1]",  # manual Ito correction
    diffusion = "theta[2] * x[1]",
    data = y_s, times = times_s, x0 = y_s[1],
    obs.model = "gaussian",
    obs.noise = function(theta) theta[3],
    prior = function(theta) {
        if (theta[2] <= 0 || theta[3] <= 0) return(-Inf)
        dnorm(theta[1], 0, 1, log = TRUE) +
            dnorm(theta[2], 0, 1, log = TRUE) +
            dexp(theta[3], 5, log = TRUE)
    },
    parm.names = c("mu", "sigma", "obs_sd"),
    method = "ekf",
    interpretation = "ito"
)

# Both should produce identical log-posteriors
theta_test <- c(0.1, 0.3, 0.2)
lp_strat <- sde_strat$Model(theta_test, sde_strat$Data)$LP
lp_ito <- sde_ito$Model(theta_test, sde_ito$Data)$LP
knitr::kable(data.frame(
    Specification = c("Stratonovich (auto-corrected)", "Ito (manual correction)"),
    `Log-posterior` = round(c(lp_strat, lp_ito), 6),
    Identical = c("", ifelse(abs(lp_strat - lp_ito) < 1e-6, "Yes", "No")),
    check.names = FALSE
), caption = "Stratonovich vs. manual Ito correction: identical log-posteriors")
```

The two specifications produce identical log-posteriors because lucifer
symbolically differentiates \\g(x) = \sigma x\\ to obtain \\g'(x) =
\sigma\\, then adds \\\frac{1}{2}\sigma x \cdot \sigma =
\frac{1}{2}\sigma^2 x\\ to the drift, which matches the manual
correction exactly.

### Worked example: effect of interpretation on inference

The Ito-Stratonovich distinction has practical consequences for
parameter interpretation. Consider a population model with
multiplicative noise: under the Stratonovich interpretation, the drift
parameter \\\mu\\ represents the “physical” growth rate (what you would
measure in a deterministic experiment), while under Ito, the effective
growth rate is \\\mu - \frac{1}{2}\sigma^2\\ (the noise-induced drift
reduces net growth). Fitting the wrong interpretation biases the growth
rate estimate by \\\frac{1}{2}\sigma^2\\.

``` r

# Fit both interpretations to the same data
fit_strat <- SDE.fit(sde_strat, Algorithm = "HARM",
                      Iterations = 10000, Thinning = 3)
fit_ito <- SDE.fit(sde_ito, Algorithm = "HARM",
                    Iterations = 10000, Thinning = 3)

# Compare posterior means
post_strat <- colMeans(fit_strat$Posterior1)
post_ito <- colMeans(fit_ito$Posterior1)
knitr::kable(data.frame(
    Specification = c("Stratonovich (auto)", "Ito (manual)"),
    `Posterior mean mu` = round(c(post_strat[1], post_ito[1]), 4),
    `Posterior mean sigma` = round(c(post_strat[2], post_ito[2]), 4),
    `Posterior mean obs_sd` = round(c(post_strat[3], post_ito[3]), 4),
    check.names = FALSE
), caption = "Posterior parameter comparison: Stratonovich vs. manual Ito")
```

``` r

# Compare posterior distributions of mu
require_package("ggplot2")
post_df <- data.frame(
    mu = c(fit_strat$Posterior1[, 1], fit_ito$Posterior1[, 1]),
    model = rep(c("Stratonovich (auto-corrected)",
                  "Ito (manual correction)"),
                each = nrow(fit_strat$Posterior1))
)
ggplot2::ggplot(post_df, ggplot2::aes(x = .data$mu, fill = .data$model)) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::scale_fill_manual(values = c(
        "Stratonovich (auto-corrected)" = "#2196F3",
        "Ito (manual correction)" = "#FF9800")) +
    ggplot2::geom_vline(xintercept = true_mu, linetype = "dashed",
        color = "black") +
    ggplot2::annotate("text", x = true_mu + 0.02, y = 0,
        label = paste0("true mu = ", true_mu),
        hjust = 0, vjust = -0.5) +
    ggplot2::labs(x = expression(mu), y = "Posterior density",
        title = expression(paste(
            "Posterior of ", mu,
            ": Stratonovich vs. manual It\u00f4 correction")),
        fill = "Interpretation") +
    theme_relab()
```

The two posterior distributions should be nearly identical, confirming
that the automatic Stratonovich correction produces the same inference
as the manual approach. The practical advantage is that users working
with physical or ecological models in Stratonovich form can specify
their model naturally without worrying about the noise-induced drift
correction.

### Trajectory comparison

``` r

# Show both fits' trajectories on the data
plot(fit_strat, type = "trajectory")
```

### When no correction is needed

For additive noise (\\g\\ does not depend on \\x\\), the Ito and
Stratonovich interpretations are identical. lucifer’s symbolic
differentiator detects this case (\\g' = 0\\) and skips the correction:

``` r

# Additive noise: g(x) = sigma (constant, independent of x)
sde_add <- SDE(
    drift = "theta[1] * (theta[2] - x[1])",  # OU drift
    diffusion = "theta[3]",                    # additive noise
    data = y_ou, times = times, x0 = y_ou[1],
    obs.model = "gaussian",
    obs.noise = function(theta) theta[4],
    prior = function(theta) {
        if (theta[1] <= 0 || theta[3] <= 0 || theta[4] <= 0) return(-Inf)
        sum(dnorm(theta, 0, 10, log = TRUE))
    },
    parm.names = c("kappa", "mu", "sigma", "obs_sd"),
    method = "ekf",
    interpretation = "stratonovich"  # no correction needed: g' = 0
)

# Same LP as Ito (no correction applied)
theta_ou <- c(2, 3, 0.5, 0.1)
sde_add_ito <- SDE(
    drift = "theta[1] * (theta[2] - x[1])",
    diffusion = "theta[3]",
    data = y_ou, times = times, x0 = y_ou[1],
    obs.model = "gaussian",
    obs.noise = function(theta) theta[4],
    prior = function(theta) {
        if (theta[1] <= 0 || theta[3] <= 0 || theta[4] <= 0) return(-Inf)
        sum(dnorm(theta, 0, 10, log = TRUE))
    },
    parm.names = c("kappa", "mu", "sigma", "obs_sd"),
    method = "ekf", interpretation = "ito"
)
lp_add_s <- sde_add$Model(theta_ou, sde_add$Data)$LP
lp_add_i <- sde_add_ito$Model(theta_ou, sde_add_ito$Data)$LP
knitr::kable(data.frame(
    Interpretation = c("Stratonovich", "Ito"),
    `Log-posterior` = round(c(lp_add_s, lp_add_i), 6),
    `Correction applied` = c("None (g' = 0)", "None"),
    check.names = FALSE
), caption = "Additive noise: Stratonovich = Ito (no correction needed)")
```

## References

**\[1\]** Oksendal, B. (2003). *Stochastic Differential Equations: An
Introduction with Applications*, 6th ed. Springer.
[doi:10.1007/978-3-642-14394-6](https://doi.org/10.1007/978-3-642-14394-6)

**\[2\]** Kloeden, P.E. and Platen, E. (1992). *Numerical Solution of
Stochastic Differential Equations*. Springer.
[doi:10.1007/978-3-662-12616-5](https://doi.org/10.1007/978-3-662-12616-5)

**\[3\]** Iacus, S.M. (2008). *Simulation and Inference for Stochastic
Differential Equations: With R Examples*. Springer.
[doi:10.1007/978-0-387-75839-8](https://doi.org/10.1007/978-0-387-75839-8)

**\[4\]** Brouste, A., Fukasawa, M., Hino, H., Iacus, S.M., Kamatani,
K., Koike, Y., Masuda, H., Nomura, R., Ogihara, T., Shimuzu, Y., Uchida,
M., and Yoshida, N. (2014). “The YUIMA project: a computational
framework for simulation and inference of stochastic differential
equations.” *Journal of Statistical Software*, 57(4), p. 1–51.
[doi:10.18637/jss.v057.i04](https://doi.org/10.18637/jss.v057.i04)

**\[5\]** King, A.A., Nguyen, D., and Ionides, E.L. (2016). “Statistical
inference for partially observed Markov processes via the R package
pomp.” *Journal of Statistical Software*, 69(12), p. 1–43.
[doi:10.18637/jss.v069.i12](https://doi.org/10.18637/jss.v069.i12)

**\[6\]** Fuchs, C. (2013). *Inference for Diffusion Processes: With
Applications in Life Sciences*. Springer.
[doi:10.1007/978-3-642-25969-2](https://doi.org/10.1007/978-3-642-25969-2)

**\[7\]** Ito, K. (1944). “Stochastic integral.” *Proceedings of the
Imperial Academy*, 20(8), p. 519–524.
[doi:10.3792/pia/1195572786](https://doi.org/10.3792/pia/1195572786)

**\[8\]** Risken, H. (1989). *The Fokker-Planck Equation: Methods of
Solution and Applications*, 2nd ed. Springer.
[doi:10.1007/978-3-642-61544-3](https://doi.org/10.1007/978-3-642-61544-3)

**\[9\]** Uhlenbeck, G.E. and Ornstein, L.S. (1930). “On the theory of
the Brownian motion.” *Physical Review*, 36(5), p. 823–841.
[doi:10.1103/PhysRev.36.823](https://doi.org/10.1103/PhysRev.36.823)

**\[10\]** Cox, J.C., Ingersoll, J.E., and Ross, S.A. (1985). “A theory
of the term structure of interest rates.” *Econometrica*, 53(2),
p. 385–407. [doi:10.2307/1911242](https://doi.org/10.2307/1911242)

**\[11\]** Sarkka, S. (2013). *Bayesian Filtering and Smoothing*.
Cambridge University Press.
[doi:10.1017/CBO9781139344203](https://doi.org/10.1017/CBO9781139344203)

**\[12\]** Doucet, A., de Freitas, N., and Gordon, N.J., eds. (2001).
*Sequential Monte Carlo Methods in Practice*. Springer.
[doi:10.1007/978-1-4757-3437-9](https://doi.org/10.1007/978-1-4757-3437-9)

**\[13\]** Gordon, N.J., Salmond, D.J., and Smith, A.F.M. (1993). “Novel
approach to nonlinear/non-Gaussian Bayesian state estimation.” *IEE
Proceedings F*, 140(2), p. 107–113.
[doi:10.1049/ip-f-2.1993.0015](https://doi.org/10.1049/ip-f-2.1993.0015)

**\[14\]** Andrieu, C., Doucet, A., and Holenstein, R. (2010). “Particle
Markov chain Monte Carlo methods.” *Journal of the Royal Statistical
Society, Series B*, 72(3), p. 269–342.
[doi:10.1111/j.1467-9868.2009.00736.x](https://doi.org/10.1111/j.1467-9868.2009.00736.x)

**\[15\]** Pitt, M.K. and Shephard, N. (1999). “Filtering via
simulation: auxiliary particle filters.” *Journal of the American
Statistical Association*, 94(446), p. 590–599.
[doi:10.1080/01621459.1999.10474153](https://doi.org/10.1080/01621459.1999.10474153)

**\[16\]** Delyon, B. and Hu, Y. (2006). “Simulation of conditioned
diffusion and application to parameter estimation.” *Stochastic
Processes and their Applications*, 116(11), p. 1660–1675.
[doi:10.1016/j.spa.2006.04.004](https://doi.org/10.1016/j.spa.2006.04.004)

**\[17\]** Kramers, H.A. (1940). “Brownian motion in a field of force
and the diffusion model of chemical reactions.” *Physica*, 7(4),
p. 284–304.
[doi:10.1016/S0031-8914(40)90098-2](https://doi.org/10.1016/S0031-8914(40)90098-2)

**\[18\]** Kallioinen, N., Paananen, T., Buerkner, P.C., and Vehtari, A.
(2024). “Detecting and diagnosing prior and likelihood sensitivity with
power-scaling.” *Statistics and Computing*, 34, 57.
[doi:10.1007/s11222-023-10366-5](https://doi.org/10.1007/s11222-023-10366-5)

**\[19\]** Merton, R.C. (1976). “Option pricing when underlying stock
returns are discontinuous.” *Journal of Financial Economics*, 3(1–2),
p. 125–144.
[doi:10.1016/0304-405X(76)90022-2](https://doi.org/10.1016/0304-405X(76)90022-2)

**\[20\]** Kou, S.G. (2002). “A jump-diffusion model for option
pricing.” *Management Science*, 48(8), p. 1086–1101.
[doi:10.1287/mnsc.48.8.1086.166](https://doi.org/10.1287/mnsc.48.8.1086.166)

**\[21\]** Heston, S.L. (1993). “A closed-form solution for options with
stochastic volatility with applications to bond and currency options.”
*Review of Financial Studies*, 6(2), p. 327–343.
[doi:10.1093/rfs/6.2.327](https://doi.org/10.1093/rfs/6.2.327)

**\[22\]** FitzHugh, R. (1961). “Impulses and physiological states in
theoretical models of nerve membrane.” *Biophysical Journal*, 1(6),
p. 445–466.
[doi:10.1016/S0006-3495(61)86902-6](https://doi.org/10.1016/S0006-3495(61)86902-6)

**\[23\]** Lotka, A.J. (1925). *Elements of Physical Biology*. Williams
& Wilkins.

**\[24\]** Volterra, V. (1926). “Fluctuations in the abundance of a
species considered mathematically.” *Nature*, 118, p. 558–560.
[doi:10.1038/118558a0](https://doi.org/10.1038/118558a0)

**\[25\]** Vehtari, A., Gelman, A., and Gabry, J. (2017). “Practical
Bayesian model evaluation using leave-one-out cross-validation and
WAIC.” *Statistics and Computing*, 27(5), p. 1413–1432.
[doi:10.1007/s11222-016-9696-4](https://doi.org/10.1007/s11222-016-9696-4)

**\[26\]** Jazwinski, A.H. (1970). *Stochastic Processes and Filtering
Theory*. Academic Press. ISBN: 978-0-12-381550-7.

**\[27\]** Wilkinson, D.J. (2018). *Stochastic Modelling for Systems
Biology*, 3rd ed. Chapman and Hall/CRC.
[doi:10.1201/9781351000918](https://doi.org/10.1201/9781351000918)

**\[28\]** Wong, E. and Zakai, M. (1965). “On the convergence of
ordinary integrals to stochastic integrals.” *Annals of Mathematical
Statistics*, 36(5), p. 1560–1564.
[doi:10.1214/aoms/1177699916](https://doi.org/10.1214/aoms/1177699916)

**\[29\]** Stratonovich, R.L. (1966). “A new representation for
stochastic integrals and equations.” *SIAM Journal on Control*, 4(2),
p. 362–371. [doi:10.1137/0304028](https://doi.org/10.1137/0304028)

**\[30\]** May, R.M. (1974). *Stability and Complexity in Model
Ecosystems*, 2nd ed. Princeton University Press. ISBN:
978-0-691-08861-7.

**\[31\]** Chesson, P. (2000). “Mechanisms of maintenance of species
diversity.” *Annual Review of Ecology and Systematics*, 31, p. 343–366.
[doi:10.1146/annurev.ecolsys.31.1.343](https://doi.org/10.1146/annurev.ecolsys.31.1.343)

**\[32\]** Gelman, A., Meng, X.-L., and Stern, H. (1996). “Posterior
predictive assessment of model fitness via realized discrepancies.”
*Statistica Sinica*, 6(4), p. 733–760.

**\[33\]** Spiegelhalter, D.J., Best, N.G., Carlin, B.P., and van der
Linde, A. (2002). “Bayesian measures of model complexity and fit.”
*JRSS-B*, 64(4), p. 583–639.
[doi:10.1111/1467-9868.00353](https://doi.org/10.1111/1467-9868.00353)

**\[34\]** Roberts, G.O. and Rosenthal, J.S. (2001). “Optimal scaling
for various Metropolis-Hastings algorithms.” *Statistical Science*,
16(4), p. 351–367.
[doi:10.1214/ss/1015346320](https://doi.org/10.1214/ss/1015346320)

**\[35\]** Bernstein, S.N. (1917). “Attempt at an axiomatic foundation
of probability theory.” *Communications of the Kharkov Mathematical
Society*, 15, p. 209–274.

**\[36\]** Hodgkin, A.L. and Huxley, A.F. (1952). “A quantitative
description of membrane current and its application to conduction and
excitation in nerve.” *Journal of Physiology*, 117(4), p. 500–544.
[doi:10.1113/jphysiol.1952.sp004764](https://doi.org/10.1113/jphysiol.1952.sp004764)

**\[37\]** Black, F. and Scholes, M. (1973). “The pricing of options and
corporate liabilities.” *Journal of Political Economy*, 81(3),
p. 637–654. [doi:10.1086/260062](https://doi.org/10.1086/260062)
