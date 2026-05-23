# State-space model inference

Jointly samples static parameters and latent state trajectories for
state-space models. Accepts either an `ssm_model` object (from
[`SSM_model`](https://robustecologies.github.io/lucifer/reference/SSM_model.md))
or an `sde_model` object (from
[`SDE`](https://robustecologies.github.io/lucifer/reference/SDE.md)) and
dispatches to a specialized algorithm that samples both the parameter
vector theta and the full latent state trajectory \\x\_{0:T}\\ at each
iteration.

## Usage

``` r
SSM(
  model,
  Algorithm = "auto",
  Iterations = 5000L,
  Thinning = 1L,
  Specs = NULL,
  Status = 100L,
  verbose = FALSE
)
```

## Arguments

- model:

  An object of class `ssm_model` (from
  [`SSM_model`](https://robustecologies.github.io/lucifer/reference/SSM_model.md))
  or `sde_model` (from
  [`SDE`](https://robustecologies.github.io/lucifer/reference/SDE.md)).

- Algorithm:

  Character string specifying the inference algorithm: `"auto"`
  (default, selects based on model structure), `"FFBS"`
  (Forward-Filtering Backward-Sampling for linear Gaussian SSMs),
  `"PGAS"` (Particle Gibbs with Ancestor Sampling for general
  nonlinear/non-Gaussian SSMs), `"SMC2"` (SMC-squared for the most
  general case), `"KSC"` (Kim-Shephard-Chib for stochastic volatility
  models), `"UKF"` (Unscented Kalman Filter for nonlinear Gaussian
  SSMs), `"EnKF"` (Ensemble Kalman Filter for high-dimensional nonlinear
  Gaussian SSMs), or `"MS-FFBS"` (Markov-switching FFBS for
  regime-switching linear Gaussian SSMs with Kim 1994 approximate
  filter).

- Iterations:

  Integer number of MCMC iterations (default 5000).

- Thinning:

  Integer thinning interval (default 1).

- Specs:

  A list of algorithm-specific tuning parameters. See details for each
  algorithm's requirements.

- Status:

  Integer controlling how often progress is printed (default 100).

- verbose:

  Logical; if `FALSE` (default), suppress progress output.

## Value

An object of class `ssm_fit`, a list containing:

- Posterior:

  Matrix of parameter samples (n_samples x dim_theta).

- States:

  Array of sampled state trajectories (n_samples x T x dim_x).

- Log.Marginal.Likelihood:

  Estimate of log p(y) (if available).

- Algorithm:

  The algorithm used.

- Summary:

  Matrix of posterior summaries.

- Minutes:

  Wall-clock time in minutes.

- call:

  The matched call.

## Details

State-space model inference

The automatic algorithm selection proceeds as follows. For `ssm_model`
objects that carry a `build_ssm` function, FFBS is selected because the
model has a linear-Gaussian structure that permits exact state sampling.
For `sde_model` objects, the selection mirrors the continuous-time
logic: stochastic volatility (`"sv"`) family triggers KSC; exact
families with Gaussian observations trigger FFBS; all other models
default to PGAS.

All algorithms perform Gibbs-style inference that alternates between
sampling the latent state trajectory conditional on parameters and
sampling parameters conditional on the state trajectory. The user must
supply a `theta_update_fn` either in the model object (via
[`SSM_model`](https://robustecologies.github.io/lucifer/reference/SSM_model.md))
or in `Specs` that draws from the full conditional distribution of the
parameters given the current state trajectory:

`theta_update_fn(parm, states, Data)`

where `parm` is the current parameter vector, `states` is a T x dim_x
matrix of the current state trajectory, and `Data` is the data list from
the model. The function must return the updated parameter vector.

## References

Carter, C.K. and Kohn, R. (1994). "On Gibbs sampling for state space
models." *Biometrika*, 81(3), p. 541–553.
[doi:10.1093/biomet/81.3.541](https://doi.org/10.1093/biomet/81.3.541)

Lindsten, F., Jordan, M.I., and Schon, T.B. (2014). "Particle Gibbs with
ancestor sampling." *JMLR*, 15(1), p. 2145–2184.

Durbin, J. and Koopman, S.J. (2012). *Time Series Analysis by State
Space Methods*. Second edition, Oxford University Press.
[doi:10.1093/acprof:oso/9780199641178.001.0001](https://doi.org/10.1093/acprof%3Aoso/9780199641178.001.0001)

## See also

[`as.demonoid.ssm_fit`](https://robustecologies.github.io/lucifer/reference/as.demonoid.ssm_fit.md),
[`plot.ssm_fit`](https://robustecologies.github.io/lucifer/reference/plot.ssm_fit.md),
[`print.ssm_fit`](https://robustecologies.github.io/lucifer/reference/print.ssm_fit.md),
[`summary.ssm_fit`](https://robustecologies.github.io/lucifer/reference/summary.ssm_fit.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Discrete-time SSM via SSM_model + FFBS
y <- cumsum(rnorm(200, sd = 0.5)) + rnorm(200, sd = 1)
mod <- SSM_model(type = "local_level", data = y)
theta_fn <- function(parm, states, Data) parm
fit <- SSM(mod, Specs = list(theta_update_fn = theta_fn))
print(fit)
summary(fit)
plot(fit)

# Continuous-time SDE via SDE() + PGAS (backward compatible)
sde <- SDE(family = "ou", data = y, times = seq_along(y))
fit2 <- SSM(sde, Algorithm = "PGAS",
            Specs = list(N.particles = 100,
                         theta_update_fn = theta_fn))
} # }
```
