# Batch simulation over forcing values

Runs parallel simulations across a range of F (thermal forcing) values
to explore the bifurcation structure and blocking frequency of the
Charney-DeVore model.

## Usage

``` r
charney_devore_batch(
  F_values,
  n_steps = 50000L,
  dt = 0.01,
  x0 = 1,
  y0 = 0.1,
  z0 = 0.1,
  k = 0.1,
  alpha = 1,
  beta = 0.5,
  delta = 1,
  sigma = 0.15,
  ou_phi = 0,
  ou_mu = 0,
  stochastic = TRUE,
  seed = NULL,
  thin = 1L,
  n_threads = 1L
)
```

## Arguments

- F_values:

  Numeric vector of thermal forcing parameter values. Values 0.5-3.0
  span the full range from stable zonal flow through the bistable regime
  to strong zonal dominance.

- n_steps:

  Number of integration steps per simulation.

- dt:

  Time step.

- x0, y0, z0:

  Initial conditions.

- k, alpha, beta, delta:

  Model parameters.

- sigma:

  Noise amplitude.

- ou_phi:

  OU relaxation rate.

- ou_mu:

  OU process mean.

- stochastic:

  Enable stochastic forcing.

- seed:

  Base random seed (each simulation uses seed + index).

- thin:

  Thinning factor.

- n_threads:

  Number of OpenMP threads for parallel execution.

## Value

Object of class "charney_devore_batch" containing simulation results for
all parameter values.

## See also

[`charney_devore_sim()`](https://robustecologies.github.io/tuRbulence/reference/charney_devore_sim.md),
[`print.charney_devore_batch()`](https://robustecologies.github.io/tuRbulence/reference/print.charney_devore_batch.md),
[`summary.charney_devore_batch()`](https://robustecologies.github.io/tuRbulence/reference/summary.charney_devore_batch.md),
[`plot.charney_devore_batch()`](https://robustecologies.github.io/tuRbulence/reference/plot.charney_devore_batch.md),
[`plot_bifurcation_panel()`](https://robustecologies.github.io/tuRbulence/reference/plot_bifurcation_panel.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Bifurcation sweep across thermal forcing values
F_seq <- seq(0.5, 3.0, length.out = 26)
batch <- charney_devore_batch(F_seq, n_steps = 100000, n_threads = 4, seed = 42)

# Visualize results
plot(batch)                                  # Bifurcation diagram
plot(batch, type = "density")                # Density ridges
plot(batch, var = "wave", type = "summary")  # Wave amplitude statistics

# Blocking frequency vs thermal forcing
# Low F: frequent blocking (weak thermal gradient)
# High F: rare blocking (strong zonal jet dominates)
plot(batch, type = "blocking")
} # }
```
