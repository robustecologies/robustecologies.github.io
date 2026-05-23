# Batch simulation over freshwater flux values

Runs parallel simulations across a range of η₂ (freshwater flux) values
to explore the bifurcation structure and bistability of the Stommel
model.

## Usage

``` r
stommel_batch(
  eta2_values,
  n_steps = 50000L,
  dt = 0.01,
  T0 = 2,
  S0 = 1,
  eta1 = 3,
  eta3 = 0.3,
  sigma_T = 0.2,
  sigma_S = 0.2,
  ou_phi = 0,
  ou_mu_T = 0,
  ou_mu_S = 0,
  stochastic = TRUE,
  seed = NULL,
  thin = 1L,
  n_threads = 1L
)
```

## Arguments

- eta2_values:

  Numeric vector of freshwater flux parameter values. Values 0.5-2.0
  span the full range from thermal-dominated to salinity-dominated
  circulation.

- n_steps:

  Number of integration steps per simulation.

- dt:

  Time step.

- T0, S0:

  Initial conditions.

- eta1, eta3:

  Model parameters.

- sigma_T, sigma_S:

  Noise amplitudes.

- ou_phi:

  OU relaxation rate.

- ou_mu_T, ou_mu_S:

  OU process means for temperature and salinity noise.

- stochastic:

  Enable stochastic forcing.

- seed:

  Base random seed (each simulation uses seed + index).

- thin:

  Thinning factor.

- n_threads:

  Number of OpenMP threads for parallel execution.

## Value

Object of class "stommel_batch" containing simulation results for all
parameter values.

## See also

[`stommel_sim()`](https://robustecologies.github.io/tuRbulence/reference/stommel_sim.md),
[`print.stommel_batch()`](https://robustecologies.github.io/tuRbulence/reference/print.stommel_batch.md),
[`summary.stommel_batch()`](https://robustecologies.github.io/tuRbulence/reference/summary.stommel_batch.md),
[`plot.stommel_batch()`](https://robustecologies.github.io/tuRbulence/reference/plot.stommel_batch.md),
[`plot_bifurcation_panel()`](https://robustecologies.github.io/tuRbulence/reference/plot_bifurcation_panel.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Bifurcation sweep across freshwater flux values
eta2_seq <- seq(0.5, 2.0, length.out = 31)
batch <- stommel_batch(eta2_seq, n_steps = 100000, n_threads = 4, seed = 42)

# Visualize bifurcation diagram
plot(batch)

# Compare thermal vs haline circulation regimes
# Low eta2: thermal circulation dominates (q > 0)
# High eta2: haline circulation possible (q < 0)

# Extract mean flow for each parameter value
mean_flows <- sapply(seq_along(eta2_seq), function(i) {
    mean(batch$q[batch$eta2 == eta2_seq[i]])
})
plot(eta2_seq, mean_flows, type = "l",
     xlab = expression(eta[2]), ylab = "Mean flow q")
abline(h = 0, lty = 2)
} # }
```
