# State-space model specification

Constructs a state-space model for inference via
[`SSM`](https://robustecologies.github.io/lucifer/reference/SSM.md).
Unlike
[`SDE`](https://robustecologies.github.io/lucifer/reference/SDE.md),
which requires drift and diffusion functions for continuous-time models,
`SSM_model` accepts discrete-time specifications directly, either
through a named builder type or a user-supplied `build_ssm` function
that returns the standard Kalman system matrices (F, Q, H, R, m0, P0) as
a function of the parameter vector.

## Usage

``` r
SSM_model(
  type = NULL,
  data,
  times = NULL,
  build_ssm = NULL,
  theta_update_fn = NULL,
  parm.names = NULL,
  state.names = NULL,
  prior = NULL,
  constraints = NULL,
  Initial.Values = NULL,
  ...
)
```

## Arguments

- type:

  Character string naming a built-in model type from the builder
  registry: `"local_level"`, `"local_linear_trend"`, `"seasonal"`,
  `"bsm"`, or `"var_p"`. When `NULL`, the user must supply `build_ssm`.

- data:

  Numeric vector or matrix of observations. For multivariate
  observations, rows are time points and columns are variables.

- times:

  Numeric vector of observation times. Defaults to `1:T` when `NULL`.

- build_ssm:

  Optional function with signature `function(parm)` returning a list
  with components `F`, `Q`, `H`, `R`, `m0`, `P0`. Matrices may be 3D
  arrays (dim_x x dim_x x T) for time-varying models, or 2D matrices for
  time-invariant models. Overrides the builder supplied by `type`.

- theta_update_fn:

  Optional function with signature `function(parm, states, Data)` that
  draws from the full conditional of the parameters given the current
  state trajectory, returning the updated parameter vector. Required for
  Gibbs-based inference; if `NULL`, it must be supplied later via
  `Specs$theta_update_fn` when calling
  [`SSM`](https://robustecologies.github.io/lucifer/reference/SSM.md).

- parm.names:

  Character vector of parameter names. Required when `type = NULL`.

- state.names:

  Character vector of state variable names. Inferred from the builder
  when `type` is specified.

- prior:

  Optional function with signature `function(theta)` returning the
  log-prior density.

- constraints:

  Optional parameter constraints, either a named list where each element
  is `c(lower, upper)` keyed by parameter name, or a matrix with
  `dim_theta` rows and 2 columns. Same format as
  [`SDE`](https://robustecologies.github.io/lucifer/reference/SDE.md).

- Initial.Values:

  Optional numeric vector of starting parameter values. When `NULL` and
  a builder is used, defaults are drawn from the builder.

- ...:

  Additional arguments passed to the builder (e.g., `period = 12` for
  seasonal models, `lags = 2` for VAR(p) models).

## Value

An object of class `ssm_model`, a list containing:

- Data:

  A list with `.y`, `.times`, `N`, `parm.names`, `mon.names`, and
  internal fields.

- build_ssm:

  The system matrix builder function.

- theta_update_fn:

  The parameter update function (may be `NULL`).

- parm.names:

  Character vector of parameter names.

- state.names:

  Character vector of state variable names.

- dim_x:

  Integer state dimension.

- dim_y:

  Integer observation dimension.

- dim_theta:

  Integer parameter dimension.

- Initial.Values:

  Numeric starting parameter values.

- constraints:

  Constraint matrix (may be `NULL`).

- prior:

  Log-prior function (may be `NULL`).

- family:

  The model type string (may be `NULL`).

## Details

State-space model specification

The `SSM_model` constructor produces an S3 object of class `ssm_model`
that [`SSM`](https://robustecologies.github.io/lucifer/reference/SSM.md)
can consume directly. When a `build_ssm` function is present, the
natural algorithm is FFBS (Forward-Filtering Backward-Sampling), which
exploits the linear-Gaussian structure for exact state sampling. Models
without `build_ssm` but with drift/diffusion functions (i.e.,
continuous-time SDE models) are handled by PGAS instead.

Built-in model types provide default system matrix builders, parameter
names, state names, initial values, and optionally a conjugate
`theta_update_fn`. User-supplied arguments always override builder
defaults.

The five built-in types cover the most common structural time series
specifications. The local level model has a single random walk state
observed with noise; the local linear trend adds a stochastic slope; the
seasonal model adds `period - 1` seasonal states; the basic structural
model (BSM) combines trend, slope, and seasonal components; and `var_p`
implements a vector autoregressive model of order p cast as an SSM.

## References

Durbin, J. and Koopman, S.J. (2012). *Time Series Analysis by State
Space Methods*. Second edition, Oxford University Press.
[doi:10.1093/acprof:oso/9780199641178.001.0001](https://doi.org/10.1093/acprof%3Aoso/9780199641178.001.0001)

Harvey, A.C. (1989). *Forecasting, Structural Time Series Models and the
Kalman Filter*. Cambridge University Press. ISBN 978-0-521-40573-7.

Carter, C.K. and Kohn, R. (1994). "On Gibbs sampling for state space
models." *Biometrika*, 81(3), p. 541–553.
[doi:10.1093/biomet/81.3.541](https://doi.org/10.1093/biomet/81.3.541)

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Local level model with built-in builder
y <- cumsum(rnorm(200, sd = 0.5)) + rnorm(200, sd = 1)
mod <- SSM_model(type = "local_level", data = y)
print(mod)

# Custom SSM with user-supplied builder
build_fn <- function(parm) {
    sigma_w <- exp(parm[1])
    sigma_v <- exp(parm[2])
    list(
        F = matrix(1, 1, 1),
        Q = matrix(sigma_w^2, 1, 1),
        H = matrix(1, 1, 1),
        R = matrix(sigma_v^2, 1, 1),
        m0 = 0,
        P0 = matrix(10, 1, 1)
    )
}
mod2 <- SSM_model(
    data = y,
    build_ssm = build_fn,
    parm.names = c("log_sigma_w", "log_sigma_v"),
    state.names = "level",
    Initial.Values = c(0, 0)
)

# Fit via SSM
theta_fn <- function(parm, states, Data) parm
fit <- SSM(mod2, Algorithm = "FFBS",
           Specs = list(theta_update_fn = theta_fn))
print(fit)
summary(fit)
plot(fit)
} # }
```
