# Simulate from an SDE model

Forward simulation from a stochastic differential equation model
specification using the Euler-Maruyama scheme or exact transition
densities for pre-built families.

## Usage

``` r
# S3 method for class 'sde_model'
simulate(
  object,
  nsim = 1,
  seed = NULL,
  theta = NULL,
  n.paths = nsim,
  times_out = NULL,
  add.obs.noise = TRUE,
  ...
)
```

## Arguments

- object:

  An object of class `sde_model`.

- nsim:

  Number of simulation paths (default 1). Passed as `n.paths` for
  consistency with base `simulate`.

- seed:

  Optional random seed for reproducibility.

- theta:

  Numeric vector of parameter values. If `NULL`, uses
  `object$Initial.Values`.

- n.paths:

  Integer number of independent paths to simulate.

- times_out:

  Optional numeric vector of output times. If `NULL`, uses the
  observation times from the model.

- add.obs.noise:

  Logical. If `TRUE` (default), observation noise is added to the
  simulated latent trajectories according to the model's observation
  model. Set to `FALSE` to return latent state trajectories only (useful
  for phase-space plots).

- ...:

  Additional arguments (unused).

## Value

A list with components:

- trajectories:

  Array of dimension (n_times, dim_x, n.paths). Contains observed-scale
  values when `add.obs.noise = TRUE`, latent states otherwise.

- times:

  Numeric vector of simulation times.

- theta:

  Parameter values used for simulation.

- state.names:

  State variable names.

## Details

Implementation of `simulate.sde_model`. Refer to the package vignettes
and the cited references for a complete algorithmic and mathematical
description.

## See also

[`SDE`](https://robustecologies.github.io/lucifer/reference/SDE.md),
[`compile.sde_model`](https://robustecologies.github.io/lucifer/reference/compile.sde_model.md),
[`plot.sde_model`](https://robustecologies.github.io/lucifer/reference/plot.sde_model.md),
[`print.sde_model`](https://robustecologies.github.io/lucifer/reference/print.sde_model.md),
[`summary.sde_model`](https://robustecologies.github.io/lucifer/reference/summary.sde_model.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Build an OU model and simulate forward paths
set.seed(42)
y <- rnorm(100, 5, 0.5)
times <- seq(0, 10, length.out = 100)
sde <- SDE(data = y, times = times, family = "ou")

# Simulate 5 trajectories with observation noise
sim <- simulate(sde, theta = c(2, 5, 0.8, 0.3), n.paths = 5)
matplot(sim$times, sim$trajectories[, 1, ], type = "l",
    xlab = "Time", ylab = "X")
points(times, y, pch = 16, cex = 0.5)

# Latent trajectories only (no obs noise)
sim_lat <- simulate(sde, theta = c(2, 5, 0.8, 0.3), n.paths = 5,
    add.obs.noise = FALSE)
} # }
```
