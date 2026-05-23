# Batch simulation over forcing values

Runs parallel simulations across a range of F (radiative forcing) values
to explore climate sensitivity and ice-albedo feedback behavior in the
Hasselmann model.

## Usage

``` r
hasselmann_batch(
  F_values,
  n_steps = 50000L,
  dt = 0.1,
  T0 = 0,
  Td0 = 0,
  I0 = 0.5,
  lambda = 0.1,
  gamma = 0.05,
  kappa = 0.01,
  alpha_ice = 0.1,
  beta_ice = 0.02,
  mu_ice = 0.5,
  sigma_T = 0.3,
  sigma_I = 0.1,
  ou_phi_T = 0,
  ou_phi_I = 0,
  ou_mu_T = 0,
  ou_mu_I = 0,
  stochastic = TRUE,
  seed = NULL,
  thin = 1L,
  n_threads = 1L
)
```

## Arguments

- F_values:

  Numeric vector of radiative forcing parameter values. Negative values
  represent cooling (e.g., volcanic aerosols), zero is the reference
  climate, positive values represent warming (e.g., CO2). F = 1
  corresponds to roughly 2×CO2 forcing (~3.7 W/m²).

- n_steps:

  Number of integration steps per simulation.

- dt:

  Time step.

- T0, Td0, I0:

  Initial conditions.

- lambda, gamma, kappa:

  Climate parameters.

- alpha_ice, beta_ice, mu_ice:

  Ice parameters.

- sigma_T, sigma_I:

  Noise amplitudes.

- ou_phi_T, ou_phi_I:

  OU relaxation rates.

- ou_mu_T, ou_mu_I:

  OU process means.

- stochastic:

  Enable stochastic forcing.

- seed:

  Base random seed (each simulation uses seed + index).

- thin:

  Thinning factor.

- n_threads:

  Number of OpenMP threads for parallel execution.

## Value

Object of class "hasselmann_batch" containing simulation results for all
forcing values.

## See also

[`hasselmann_sim()`](https://robustecologies.github.io/tuRbulence/reference/hasselmann_sim.md),
[`print.hasselmann_batch()`](https://robustecologies.github.io/tuRbulence/reference/print.hasselmann_batch.md),
[`summary.hasselmann_batch()`](https://robustecologies.github.io/tuRbulence/reference/summary.hasselmann_batch.md),
[`plot.hasselmann_batch()`](https://robustecologies.github.io/tuRbulence/reference/plot.hasselmann_batch.md),
[`plot_bifurcation_panel()`](https://robustecologies.github.io/tuRbulence/reference/plot_bifurcation_panel.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Climate sensitivity experiment: sweep forcing from cooling to warming
F_seq <- seq(-1.0, 2.0, length.out = 31)
batch <- hasselmann_batch(F_seq, n_steps = 100000, n_threads = 4, seed = 42)

# Visualize forcing-response relationship
plot(batch)

# Extract equilibrium temperature for each forcing level
# (using last half of simulation as pseudo-equilibrium)
eq_temp <- sapply(F_seq, function(f) {
    idx <- batch$F == f
    T_vals <- batch$T[idx]
    mean(T_vals[(length(T_vals)/2):length(T_vals)])
})

plot(F_seq, eq_temp, type = "l",
     xlab = "Radiative forcing F", ylab = "Equilibrium temperature anomaly")
abline(h = 0, v = 0, lty = 2)

# Estimate climate sensitivity (slope dT/dF)
sensitivity <- diff(eq_temp) / diff(F_seq)
cat("Climate sensitivity:", mean(sensitivity), "K per forcing unit\n")

# Note: nonlinearity from ice-albedo feedback visible at negative F
# where ice growth amplifies cooling
} # }
```
