# Simulate Huisman-Weissing competition dynamics

Integrates the Huisman-Weissing resource competition model and returns
time series of species abundances and resource concentrations.

## Usage

``` r
simulate_competition(
  system,
  N0 = NULL,
  R0 = NULL,
  t_max = 4000,
  dt = 0.01,
  t_out = 1
)
```

## Arguments

- system:

  A multispecies competition system object.

- N0:

  Optional numeric vector of initial species abundances. If NULL, uses
  default from system parameters.

- R0:

  Optional numeric vector of initial resource concentrations. If NULL,
  uses default from system parameters.

- t_max:

  Maximum simulation time (days).

- dt:

  Integration time step (default 0.01).

- t_out:

  Output time interval (default 1.0).

## Value

A list containing:

- time:

  Vector of time points

- N:

  Matrix of species abundances (time x species)

- R:

  Matrix of resource concentrations (time x resources)

## Examples

``` r
if (FALSE) { # \dontrun{
hw <- multispecies_competition(scenario = "5species")

# Simulate with default initial conditions
result <- simulate_competition(hw, t_max = 4000)

# Plot time series
matplot(result$time, result$N, type = "l", lty = 1,
        xlab = "Time (days)", ylab = "Species abundance",
        main = "5-species competition dynamics")
legend("topright", paste0("Sp", 1:5), col = 1:5, lty = 1)

# Different initial conditions
N0_custom <- c(0.5, 1.0, 0.1, 0.8, 0.1)
result2 <- simulate_competition(hw, N0 = N0_custom, t_max = 4000)
} # }
```
