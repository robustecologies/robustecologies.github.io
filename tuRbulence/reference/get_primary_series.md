# Extract primary time series from simulation

Extracts the primary observable time series from any turbulence
simulation object, suitable for embedding and analysis.

## Usage

``` r
get_primary_series(sim, var = NULL)
```

## Arguments

- sim:

  Simulation object from any system simulator.

- var:

  Optional variable name to extract; if NULL, uses the primary
  observable for each system.

## Value

Numeric vector of the time series.

## Details

Primary observables by system:

- vonkarman:

  theta (transformed symmetry parameter)

- stommel:

  q (flow strength = T - S)

- charney_devore:

  x (zonal flow amplitude)

- hasselmann:

  T (surface temperature anomaly)

- lorenz:

  x (convective intensity)

- rossler:

  x (first state variable)

- lorenz84:

  x (westerly wind strength)

## See also

[`turbulence_sim()`](https://robustecologies.github.io/tuRbulence/reference/turbulence_sim.md),
[`simulate_system()`](https://robustecologies.github.io/tuRbulence/reference/simulate_system.md),
[`turbulence_embed()`](https://robustecologies.github.io/tuRbulence/reference/turbulence_embed.md),
[`turbulence_lyapunov()`](https://robustecologies.github.io/tuRbulence/reference/turbulence_lyapunov.md)
