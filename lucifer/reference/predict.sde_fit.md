# Predict future trajectories from a fitted SDE model

Samples parameter vectors from the posterior distribution of a fitted
SDE model and simulates forward trajectories beyond the last
observation, producing prediction intervals that account for both
parameter uncertainty and process noise.

## Usage

``` r
# S3 method for class 'sde_fit'
predict(
  object,
  n.ahead = 20L,
  dt.ahead = NULL,
  n.paths = 100L,
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  ...
)
```

## Arguments

- object:

  An object of class `sde_fit`, as returned by
  [`SDE.fit`](https://robustecologies.github.io/lucifer/reference/SDE.fit.md).

- n.ahead:

  Integer number of time steps to forecast. Default 20.

- dt.ahead:

  Numeric time step for the forecast grid. Defaults to the median
  observation interval from the fitted model.

- n.paths:

  Integer number of posterior predictive trajectories to simulate.
  Default 100.

- probs:

  Numeric vector of quantiles for prediction intervals. Default
  `c(0.025, 0.25, 0.5, 0.75, 0.975)`.

- ...:

  Additional arguments (currently unused).

## Value

A list of class `sde_prediction` containing:

- times:

  Numeric vector of forecast times.

- trajectories:

  Array of dimension `(n.ahead, dim_x, n.paths)` with simulated paths.

- quantiles:

  Matrix of dimension `(n.ahead, length(probs))` with prediction
  quantiles for the first state variable.

- mean:

  Numeric vector of posterior predictive means.

- state.names:

  Character vector of state variable names.

## Details

Forecast from a fitted SDE model

Posterior predictive sampling from a stochastic differential equation
fit produced by
[`SDE.fit`](https://robustecologies.github.io/lucifer/reference/SDE.fit.md).
Summary of the content is given below. For each posterior draw, the user
Model is re-evaluated and `yhat` is collected. The resulting posterior
predictive matrix enables posterior predictive checks via
[`PosteriorChecks`](https://robustecologies.github.io/lucifer/reference/PosteriorChecks.md)
and diagnostic plotting via the `ppc` family of plot methods.

## References

Iacus, S. M. (2008). *Simulation and Inference for Stochastic
Differential Equations*. Springer.
[doi:10.1007/978-0-387-75839-8](https://doi.org/10.1007/978-0-387-75839-8)

## See also

[`SDE.fit`](https://robustecologies.github.io/lucifer/reference/SDE.fit.md),
[`LOO.sde_fit`](https://robustecologies.github.io/lucifer/reference/LOO.md),
[`WAIC.sde_fit`](https://robustecologies.github.io/lucifer/reference/WAIC.md),
[`log_lik.sde_fit`](https://robustecologies.github.io/lucifer/reference/log_lik.md),
[`plot.sde_fit`](https://robustecologies.github.io/lucifer/reference/plot.sde_fit.md),
[`print.sde_fit`](https://robustecologies.github.io/lucifer/reference/print.sde_fit.md),
[`summary.sde_fit`](https://robustecologies.github.io/lucifer/reference/summary.sde_fit.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate and fit an OU process
set.seed(42)
n <- 100; dt <- 0.1
x <- numeric(n); x[1] <- 5
for (i in 2:n) x[i] <- x[i-1] + 2*(5-x[i-1])*dt + 0.8*sqrt(dt)*rnorm(1)
y <- x + rnorm(n, 0, 0.3)
times <- seq(0, by = dt, length.out = n)

sde <- SDE(data = y, times = times, family = "ou")
fit <- SDE.fit(sde, Algorithm = "AMWG", Iterations = 30000,
    Thinning = 10)

# Forecast 30 steps ahead
pred <- predict(fit, n.ahead = 30, n.paths = 200)
print(pred)
plot(pred)
} # }
```
