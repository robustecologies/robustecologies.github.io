# Create bifurcation panel showing attractor evolution

Generates a multi-panel figure displaying phase space projections at
different control parameter values, illustrating the system's
bifurcation structure across all supported dynamical systems.

## Usage

``` r
plot_bifurcation_panel(
  system = c("vonkarman", "stommel", "charney_devore", "hasselmann", "lorenz", "rossler",
    "lorenz84"),
  param_values = NULL,
  n_steps = 50000L,
  ncol = 5,
  n_plot = 5000,
  stochastic = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- system:

  Character string specifying the dynamical system:

  "vonkarman"

  :   Von Kármán turbulent flow (stochastic Duffing)

  "stommel"

  :   Stommel thermohaline circulation

  "charney_devore"

  :   Charney-DeVore atmospheric blocking

  "hasselmann"

  :   Hasselmann stochastic climate

  "lorenz"

  :   Lorenz (1963) atmospheric convection

  "rossler"

  :   Rössler chaotic attractor

  "lorenz84"

  :   Lorenz-84 atmospheric model

- param_values:

  Numeric vector of control parameter values. If NULL, uses
  system-specific defaults spanning the interesting regime.

- n_steps:

  Number of integration steps per simulation.

- ncol:

  Number of columns in the panel layout.

- n_plot:

  Maximum points per panel (for performance).

- stochastic:

  Logical; enable stochastic forcing where applicable.

- verbose:

  Logical; print progress messages.

- ...:

  Additional arguments passed to the underlying simulation function.

## Value

A ggplot2 object with faceted phase space plots.

## Details

This function runs simulations at each specified parameter value and
displays 2D phase space projections in a faceted layout. The phase
variables and parameter labels are automatically selected based on the
system type.

Default parameter ranges by system:

- vonkarman:

  μ: -0.3 to 0.3 (symmetric/asymmetric transition)

- stommel:

  η₂: 0.6 to 1.4 (bistable regime)

- charney_devore:

  F: 1.0 to 2.0 (blocking transitions)

- hasselmann:

  F: -0.5 to 0.5 (forcing scenarios)

- lorenz:

  ρ: 20 to 30 (onset of chaos)

- rossler:

  c: 4 to 8 (period-doubling to chaos)

- lorenz84:

  F: 6 to 10 (seasonal forcing range)

## See also

[`turbulence_sim`](https://robustecologies.github.io/tuRbulence/reference/turbulence_sim.md)
for the unified simulation interface,
[`turbulence_bifurcation`](https://robustecologies.github.io/tuRbulence/reference/turbulence_bifurcation.md)
for Lyapunov exponent analysis.

## Examples

``` r
if (FALSE) { # \dontrun{
# Von Kármán bifurcation sequence
plot_bifurcation_panel("vonkarman",
    param_values = c(-0.2, 0, 0.1, 0.2, 0.3),
    n_steps = 30000)

# Lorenz system across Rayleigh number values
plot_bifurcation_panel("lorenz",
    param_values = seq(20, 28, by = 2),
    n_steps = 50000)

# Stommel bistability
plot_bifurcation_panel("stommel",
    param_values = seq(0.8, 1.2, by = 0.1),
    n_steps = 100000)
} # }
```
