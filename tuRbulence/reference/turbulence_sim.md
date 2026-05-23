# Simulate a turbulent dynamical system

Unified interface for simulating any dynamical system in the tuRbulence
package. All systems support both deterministic and stochastic variants.

## Usage

``` r
turbulence_sim(
  system = c("vonkarman", "stommel", "charney_devore", "hasselmann", "lorenz", "rossler",
    "lorenz84"),
  control_param,
  n_steps = 100000L,
  dt = 0.01,
  stochastic = TRUE,
  seed = NULL,
  ...
)
```

## Arguments

- system:

  Character string specifying the system to simulate:

  "vonkarman"

  :   Von Kármán turbulent flow (stochastic Duffing oscillator)

  "stommel"

  :   Stommel thermohaline circulation model

  "charney_devore"

  :   Charney-DeVore atmospheric blocking model

  "hasselmann"

  :   Hasselmann stochastic climate model

  "lorenz"

  :   Lorenz (1963) atmospheric convection model

  "rossler"

  :   Rössler (1976) chaotic attractor

  "lorenz84"

  :   Lorenz-84 low-order atmospheric model

- control_param:

  Control parameter value. Meaning depends on system:

  vonkarman

  :   \\\mu\\, asymmetry parameter (default 0.3)

  stommel

  :   \\\eta_2\\, freshwater flux (default 1.0)

  charney_devore

  :   F, thermal forcing (default 1.5)

  hasselmann

  :   F, external radiative forcing (default 0.0)

  lorenz

  :   \\\rho\\, Rayleigh number (default 28)

  rossler

  :   c, third parameter (default 5.7)

  lorenz84

  :   F, symmetric thermal forcing (default 8.0)

- n_steps:

  Number of integration steps.

- dt:

  Time step.

- stochastic:

  Logical; if TRUE, enable stochastic forcing. Note: vonkarman is
  inherently stochastic; lorenz, rossler, lorenz84 use optional
  Ornstein-Uhlenbeck noise.

- seed:

  Random seed for reproducibility.

- ...:

  Additional system-specific parameters passed to the underlying
  simulation function.

## Value

Object of the appropriate simulation class.

## Details

This function provides a unified entry point for all dynamical systems
in the package. Each system has its own control parameter with physical
meaning.

The geophysical models (vonkarman, stommel, charney_devore, hasselmann)
default to stochastic forcing. The classical chaotic systems (lorenz,
rossler, lorenz84) support optional stochastic forcing via
Ornstein-Uhlenbeck noise but default to deterministic dynamics.

Reference Lyapunov exponents for standard parameters:

- Lorenz (\\\sigma\\=10, \\\rho\\=28, \\\beta\\=8/3): \\\lambda_1
  \approx 0.906\\

- Rössler (a=0.2, b=0.2, c=5.7): \\\lambda_1 \approx 0.07\\

- Lorenz-84: \\\lambda_1 \approx 0.07\\

## See also

[`simulate_system()`](https://robustecologies.github.io/tuRbulence/reference/simulate_system.md),
[`get_primary_series()`](https://robustecologies.github.io/tuRbulence/reference/get_primary_series.md),
[`lorenz_sim()`](https://robustecologies.github.io/tuRbulence/reference/lorenz_sim.md),
[`vonkarman_sim()`](https://robustecologies.github.io/tuRbulence/reference/vonkarman_sim.md),
[`hasselmann_sim()`](https://robustecologies.github.io/tuRbulence/reference/hasselmann_sim.md),
[`turbulence_lyapunov()`](https://robustecologies.github.io/tuRbulence/reference/turbulence_lyapunov.md),
[`turbulence_cao()`](https://robustecologies.github.io/tuRbulence/reference/turbulence_cao.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Von Kármán chaotic attractor
sim1 <- turbulence_sim("vonkarman", control_param = 0.3)
print(sim1)
plot(sim1)

# Lorenz butterfly attractor
sim2 <- turbulence_sim("lorenz", control_param = 28)
plot(sim2, type = "attractor")

# Stommel bistable circulation
sim3 <- turbulence_sim("stommel", control_param = 1.0)
summary(sim3)

# Rössler attractor
sim4 <- turbulence_sim("rossler", control_param = 5.7)
plot(sim4, type = "phase")
} # }
```
