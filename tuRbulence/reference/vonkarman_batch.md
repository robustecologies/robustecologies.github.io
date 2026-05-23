# Batch simulation over multiple control parameter values

Runs parallel simulations across a range of μ values to explore the
bifurcation structure of the stochastic Duffing attractor representing
von Kármán turbulent flow dynamics.

## Usage

``` r
vonkarman_batch(
  mu_values,
  n_steps = 50000L,
  dt = 0.01,
  x0 = 0.1,
  y0 = 0,
  z0 = 0,
  a = 0.2,
  phi = 0.9,
  sigma = 0.2,
  omega = 1,
  seed = NULL,
  thin = 1L,
  n_threads = 1L
)
```

## Arguments

- mu_values:

  Numeric vector of control parameter values. Values near 0 produce
  random point attractors; ±0.02-0.04 gives noisy periodic motion; \|μ\|
  \> 0.06 produces chaotic attractors.

- n_steps:

  Number of integration steps per simulation.

- dt:

  Time step.

- x0, y0, z0:

  Initial conditions.

- a, phi, sigma, omega:

  Model parameters (see
  [`vonkarman_sim`](https://robustecologies.github.io/tuRbulence/reference/vonkarman_sim.md)).

- seed:

  Base random seed (each simulation uses seed + index).

- thin:

  Thinning factor.

- n_threads:

  Number of OpenMP threads for parallel execution.

## Value

Object of class "vonkarman_batch" containing simulation results for all
parameter values.

## See also

[`vonkarman_sim()`](https://robustecologies.github.io/tuRbulence/reference/vonkarman_sim.md),
[`print.vonkarman_batch()`](https://robustecologies.github.io/tuRbulence/reference/print.vonkarman_batch.md),
[`summary.vonkarman_batch()`](https://robustecologies.github.io/tuRbulence/reference/summary.vonkarman_batch.md),
[`plot.vonkarman_batch()`](https://robustecologies.github.io/tuRbulence/reference/plot.vonkarman_batch.md),
[`plot_bifurcation_panel()`](https://robustecologies.github.io/tuRbulence/reference/plot_bifurcation_panel.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Bifurcation sweep across control parameter
mu_seq <- seq(-0.3, 0.3, length.out = 31)
batch <- vonkarman_batch(mu_seq, n_steps = 100000, n_threads = 4, seed = 42)

# Visualize bifurcation diagram
plot(batch)

# Analyze attractor characteristics vs control parameter
# Compute variance (proxy for attractor size) at each μ
var_x <- sapply(mu_seq, function(m) {
    var(batch$x[batch$mu == m])
})

plot(mu_seq, var_x, type = "l",
     xlab = expression(mu), ylab = "Variance of x")

# Near μ = 0: random switching between point attractors (high variance)
# Intermediate μ: periodic orbits (moderate variance)
# Large |μ|: chaotic attractor (structured variance)

# Use plot_bifurcation_panel for publication-quality visualization
plot_bifurcation_panel(batch)
} # }
```
