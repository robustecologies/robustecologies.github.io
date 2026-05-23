# Compute Lyapunov exponents across parameter range

Performs bifurcation analysis by computing Lyapunov exponents for
multiple values of the control parameter.

## Usage

``` r
turbulence_bifurcation(
  system = c("vonkarman", "stommel", "charney_devore", "hasselmann", "lorenz", "rossler",
    "lorenz84"),
  param_values,
  n_steps = 100000L,
  embed_dim = 10L,
  tau = 10L,
  stochastic = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- system:

  Character string specifying the system.

- param_values:

  Numeric vector of control parameter values.

- n_steps:

  Simulation length for each parameter value.

- embed_dim:

  Embedding dimension for Lyapunov estimation.

- tau:

  Time delay.

- stochastic:

  Enable stochastic forcing (default TRUE).

- verbose:

  Print progress messages.

- ...:

  Additional arguments passed to simulation functions.

## Value

Object of class "turbulence_bifurcation" containing parameter values,
Lyapunov exponents, and metadata.

## Details

This function performs parameter sweeps across the control parameter,
computing Lyapunov exponents at each value to identify bifurcation
points and regime transitions. Positive \\\lambda_1\\ indicates chaos;
negative indicates stable fixed points or limit cycles.

## See also

[`print.turbulence_bifurcation()`](https://robustecologies.github.io/tuRbulence/reference/print.turbulence_bifurcation.md),
[`summary.turbulence_bifurcation()`](https://robustecologies.github.io/tuRbulence/reference/summary.turbulence_bifurcation.md),
[`plot.turbulence_bifurcation()`](https://robustecologies.github.io/tuRbulence/reference/plot.turbulence_bifurcation.md),
[`turbulence_lyapunov()`](https://robustecologies.github.io/tuRbulence/reference/turbulence_lyapunov.md),
[`plot_bifurcation_panel()`](https://robustecologies.github.io/tuRbulence/reference/plot_bifurcation_panel.md),
[`turbulence_sim()`](https://robustecologies.github.io/tuRbulence/reference/turbulence_sim.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Von Kármán bifurcation analysis
bif <- turbulence_bifurcation(
    system = "vonkarman",
    param_values = seq(0.1, 0.5, by = 0.05),
    n_steps = 50000,
    verbose = TRUE
)
print(bif)
summary(bif)

# Visualize bifurcation diagram
plot(bif)
} # }
```
