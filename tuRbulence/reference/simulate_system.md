# Simulate a dynamical system (legacy interface)

This function provides a legacy interface for simulating dynamical
systems. New code should use
[`turbulence_sim`](https://robustecologies.github.io/tuRbulence/reference/turbulence_sim.md)
instead, which provides a unified interface for all systems in the
package.

## Usage

``` r
simulate_system(
  system = c("vonkarman", "lorenz", "rossler", "lorenz84", "stommel", "charney_devore",
    "hasselmann"),
  params = list(),
  n_steps = 100000L,
  dt = 0.01,
  seed = NULL,
  thin = 1L,
  transient = 0L,
  stochastic = FALSE,
  ...
)
```

## Arguments

- system:

  Character string specifying the system to simulate.

- params:

  Named list of system-specific parameters.

- n_steps:

  Number of integration steps.

- dt:

  Time step.

- seed:

  Random seed.

- thin:

  Thinning factor.

- transient:

  Number of transient steps to discard.

- stochastic:

  Logical; add stochastic forcing.

- ...:

  Additional arguments passed to the specific simulator.

## Value

Object of class "turbulence_system" with system-specific subclass.

## See also

[`turbulence_sim`](https://robustecologies.github.io/tuRbulence/reference/turbulence_sim.md)
for the recommended unified interface.
