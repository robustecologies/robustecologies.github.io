# Forecast from a fitted neural ODE model

Forward simulation via fourth-order Runge-Kutta (RK4) integration using
the fitted process model ensemble. The integration operates in log space
to enforce positivity and returns predictions on the original scale. The
result is a `node_forecast` object with `print` and `plot` methods.

S3 method: apply [`predict()`](https://rdrr.io/r/stats/predict.html) to
objects of class `node_fit`.

## Usage

``` r
NODE_predict(object, horizon, n.steps = 100L, ...)

# S3 method for class 'node_fit'
predict(object, ...)
```

## Arguments

- object:

  See Details.

- horizon:

  Numeric. Time units to forecast beyond the last observation.

- n.steps:

  Integer. Number of output time steps in the forecast trajectory.
  Default 100.

- ...:

  See Details.

## Value

An object of class `node_forecast`, a list containing:

- times:

  Numeric vector of forecast time points.

- states:

  Matrix (n.steps x n_vars) of ensemble mean predicted values on the
  original (untransformed) scale.

- states_sd:

  Matrix (n.steps x n_vars) of ensemble standard deviation at each
  forecast step.

- states_lower:

  Matrix (n.steps x n_vars) of 5th percentile (90% CI lower bound) from
  the ensemble.

- states_upper:

  Matrix (n.steps x n_vars) of 95th percentile (90% CI upper bound) from
  the ensemble.

- x0:

  Named numeric vector of the initial state (last fitted value).

- obs_times:

  Numeric vector of observed time points (from the parent fit).

- obs_data:

  Matrix of original observations (from the parent fit).

- obs_fitted:

  Matrix of fitted values (from the parent fit).

- var.names:

  Character vector of variable names.

- n_vars:

  Number of variables.

- horizon:

  Forecast horizon used.

See Details.

## Details

The forecast starts from the last interpolated state of the fitted model
(not the last raw observation), ensuring continuity between the fitted
trajectory and the forecast. The process model ensemble provides
multiple trajectory samples; the returned states are the ensemble mean
prediction. Trajectories are clamped to prevent divergence when the
fitted dynamics have unrealistic eigenvalues far from the training data.

Implementation of `predict.node_fit`. Refer to the package vignettes and
the cited references for a complete algorithmic and mathematical
description.

## See also

[`NODE`](https://robustecologies.github.io/lucifer/reference/NODE.md)
for model fitting,
[`plot.node_forecast`](https://robustecologies.github.io/lucifer/reference/plot.node_forecast.md)
for visualization.

[`NODE`](https://robustecologies.github.io/lucifer/reference/NODE.md),
[`plot.node_fit`](https://robustecologies.github.io/lucifer/reference/plot.node_fit.md),
[`print.node_fit`](https://robustecologies.github.io/lucifer/reference/print.node_fit.md),
[`summary.node_fit`](https://robustecologies.github.io/lucifer/reference/summary.node_fit.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate and fit
set.seed(42)
n <- 100; dt <- 0.1
prey <- predator <- numeric(n)
prey[1] <- 1.0; predator[1] <- 0.5
for (i in 2:n) {
  prey[i] <- prey[i-1] + dt * prey[i-1] * (1 - 0.5 * predator[i-1])
  predator[i] <- predator[i-1] + dt * predator[i-1] * (-0.5 + 0.3 * prey[i-1])
}
fit <- NODE(cbind(prey = prey, predator = predator),
            times = seq(dt, n * dt, by = dt),
            obs.ensemble = 50, proc.ensemble = 50)

# Forecast 5 time units ahead
pred <- NODE_predict(fit, horizon = 5, n.steps = 80)

# Inspect
print(pred)

# Trajectory plot with forecast zone shading
plot(pred)

# Phase space with forecast extension
plot(pred, type = "phase")
} # }

if (FALSE) { # \dontrun{
## see package vignettes for full examples involving predict.node_fit
} # }
```
