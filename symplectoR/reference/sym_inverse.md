# Build an objective for inverse modeling (parameter estimation)

Normalizes a model-plus-data parameter estimation problem into the
[`sym_objective()`](https://robustecologies.github.io/symplectoR/reference/sym_objective.md)
currency consumed by
[`sym_optim()`](https://robustecologies.github.io/symplectoR/reference/sym_optim.md).
Two model forms are supported: a plain prediction function of the
parameter vector, and a janos `system_spec` dynamical system fitted by
trajectory matching.

## Usage

``` r
sym_inverse(
  model,
  data,
  times = NULL,
  loss = c("lsq", "nll", "custom"),
  loss_fn = NULL,
  theta_bounds = NULL,
  obs_sigma = NULL,
  solver = NULL,
  weights = NULL,
  theta_init = NULL,
  constraint = c("reflect", "barrier"),
  obs_times = NULL,
  ...
)

# S3 method for class '`function`'
sym_inverse(
  model,
  data,
  times = NULL,
  loss = c("lsq", "nll", "custom"),
  loss_fn = NULL,
  theta_bounds = NULL,
  obs_sigma = NULL,
  solver = NULL,
  weights = NULL,
  theta_init = NULL,
  constraint = c("reflect", "barrier"),
  obs_times = NULL,
  ...
)

# S3 method for class 'system_spec'
sym_inverse(
  model,
  data,
  times = NULL,
  loss = c("lsq", "nll", "custom"),
  loss_fn = NULL,
  theta_bounds = NULL,
  obs_sigma = NULL,
  solver = NULL,
  weights = NULL,
  theta_init = NULL,
  constraint = c("reflect", "barrier"),
  obs_times = NULL,
  ...
)

# Default S3 method
sym_inverse(model, data, ...)
```

## Arguments

- model:

  Either a function `model(theta)` returning predictions conformable
  with `data`, or a
  [`janos::system_spec`](https://robustecologies.github.io/janos/reference/system_spec.html)
  whose simulated trajectory is matched to observed time series.

- data:

  Observations. For the function form, a numeric vector or matrix
  conformable with the model output. For the `system_spec` form, a data
  frame with a `time` column and one column per observed state variable
  (names matching the specification's state names).

- times:

  Observation times for the function form when the model needs them
  (passed through as the second argument of `model` when supplied).
  Ignored for the `system_spec` form, which reads times from `data`.

- loss:

  Loss functional. `"lsq"` is (weighted) least squares; `"nll"` the
  Gaussian negative log-likelihood with observation standard deviation
  `obs_sigma`; `"custom"` uses `loss_fn`.

- loss_fn:

  For `loss = "custom"`: a function `loss_fn(pred, data)` returning a
  scalar.

- theta_bounds:

  Named parameter box as `list(lo = c(...), hi = c(...))` (the
  convention of
  [`RElabverse::stiff_globopt()`](https://rdrr.io/pkg/RElabverse/man/stiff_globopt.html)).
  The names define which parameters are estimated and label the fit
  output. Required for the `system_spec` form; for the function form
  either `theta_bounds` or `theta_init` must be supplied to fix the
  dimension.

- obs_sigma:

  Observation standard deviation for `loss = "nll"`.

- solver:

  Optional janos solver specification (for example
  [`janos::solver_rk45()`](https://robustecologies.github.io/janos/reference/solver_rk45.html))
  overriding the default
  [`janos::solver_rk4()`](https://robustecologies.github.io/janos/reference/solver_rk4.html).

- weights:

  Optional non-negative weights conformable with `data` for
  `loss = "lsq"`.

- theta_init:

  Optional named parameter vector fixing the dimension and names when
  `theta_bounds` is absent (function form only).

- constraint:

  Constraint mode when `theta_bounds` is present: `"reflect"` (default;
  iterates fold back into the box and the objective is evaluated
  exactly) or `"barrier"` (smooth logarithmic barrier).

- obs_times:

  Optional numeric vector of observation times used only to label the
  horizontal axis of the `"observed"` plot view; unlike `times` it is
  never passed to `model`, so it is the argument to use when the
  prediction function already knows its own time grid. Defaults to
  `times` for the function form and to `data$time` for the `system_spec`
  form; when neither is available the observation index is used.

- ...:

  Reserved.

## Value

A `sym_objective` over the parameter vector, with names and box taken
from `theta_bounds`, ready for
[`sym_optim()`](https://robustecologies.github.io/symplectoR/reference/sym_optim.md);
gradients are supplied by central finite differences. The objective also
remembers the observations and the prediction map, so that any
[`sym_optim()`](https://robustecologies.github.io/symplectoR/reference/sym_optim.md)
fit of it supports the `"observed"` and `"residuals"` views of
[`plot.sym_fit()`](https://robustecologies.github.io/symplectoR/reference/plot.sym_fit.md)
and leads its dashboard with them.

## Details

Trajectory matching integrates the system at the candidate parameters
with the transient discard disabled, linearly interpolates the simulated
states at the observation times, and applies the loss to the observed
columns. Quantum Hamiltonian descent (`method = "qhd"` in
[`sym_optim()`](https://robustecologies.github.io/symplectoR/reference/sym_optim.md))
is the natural global engine for likelihood surfaces with at most three
free parameters; the trajectory methods handle any dimension and refine
local minima.

## References

Ramsay, J. O., Hooker, G., Campbell, D., & Cao, J. (2007). Parameter
estimation for differential equations: a generalized smoothing approach.
*Journal of the Royal Statistical Society: Series B*, 69(5), 741-796.
[doi:10.1111/j.1467-9868.2007.00610.x](https://doi.org/10.1111/j.1467-9868.2007.00610.x)

## See also

[`sym_objective()`](https://robustecologies.github.io/symplectoR/reference/sym_objective.md),
[`sym_optim()`](https://robustecologies.github.io/symplectoR/reference/sym_optim.md),
[`as_incumbent_solver()`](https://robustecologies.github.io/symplectoR/reference/as_incumbent_solver.md),
[`sym_control()`](https://robustecologies.github.io/symplectoR/reference/sym_control.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## Function form on a shipped record: the Nicholson blowfly delayed-recruitment
## model, estimated by conditional one-step-ahead prediction on the log scale.
main <- subset(nicholson_blowfly, population == "main")
N <- main$count
k <- 7L                                  # 7 sampling steps x 2 days = the development delay
idx <- seq.int(k + 1L, length(N) - 1L)
predict_step <- function(theta) {
  log(exp(theta[1]) * N[idx - k] * exp(-N[idx - k] / exp(theta[2])) +
      N[idx] * exp(-exp(theta[3])) + 1)
}
obj <- sym_inverse(predict_step, data = log(N[idx + 1L] + 1),
                   obs_times = main$day[idx + 1L],
                   theta_bounds = list(
                     lo = c(logP = log(0.05), logN0 = log(50),    logdelta = log(0.005)),
                     hi = c(logP = log(50),   logN0 = log(20000), logdelta = log(3))))

## Three free parameters is the dimension ceiling of the quantum solver, so
## the whole box can be scanned globally before any local method is started.
global  <- sym_optim(obj, method = "qhd", seed = 1,
                     control = sym_control("qhd", N_grid = 48, K = 1200))
refined <- sym_optim(obj, x0 = global$x_best, method = "slc_expo",
                     control = sym_control("slc_expo", C = 0.1, h = 1.5,
                                           max_iter = 3000, tol_grad = 1e-9))
summary(refined)
exp(refined$x_best)                      # recruitment P, scale N0, mortality delta

## Fits of a sym_inverse() objective remember their data, so they plot
## against it: observations and predictions, residuals, and a composite.
plot(refined, type = "observed")
plot(refined, type = "residuals")
plot(refined, type = "dashboard")

## Trajectory-matching form on the shipped predator-prey microcosm, through
## a janos system. Rate parameters are estimated on the logarithmic scale,
## because their natural magnitudes differ by orders of magnitude.
d <- paramecium_didinium
lv <- janos::system_spec(
  rhs = list(paramecium ~ exp(lr) * paramecium - exp(la) * paramecium * didinium,
             didinium   ~ exp(lb) * exp(la) * paramecium * didinium - exp(lm) * didinium),
  state_names = c("paramecium", "didinium"),
  parms = list(lr = 0, la = -4, lb = -0.7, lm = 0),
  init = c(paramecium = d$paramecium[1], didinium = d$didinium[1]))
obj2 <- sym_inverse(lv, data = d, loss = "custom",
                    loss_fn = function(pred, data) sum((log(pred + 1) - log(data + 1))^2),
                    theta_bounds = list(lo = c(lr = -3, la = -9, lb = -5, lm = -3),
                                        hi = c(lr =  3, la =  0, lb =  2, lm =  3)))
fit2 <- sym_optim(obj2, x0 = c(lr = 0, la = -4, lb = -0.7, lm = 0), method = "slc_expo")
plot(fit2, type = "observed")            # both species against the fitted trajectory
plot(fit2, type = "dashboard")
} # }
```
