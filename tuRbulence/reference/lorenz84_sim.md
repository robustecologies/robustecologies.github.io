# Simulate the Lorenz-84 atmospheric model

Integrates the Lorenz-84 low-order atmospheric circulation model using
4th-order Runge-Kutta, with optional stochastic forcing.

## Usage

``` r
lorenz84_sim(
  a = 0.25,
  b = 4,
  F = 8,
  G = 1,
  n_steps = 500000L,
  dt = 0.001,
  x0 = 1,
  y0 = 0,
  z0 = 0,
  seed = NULL,
  thin = 5L,
  transient = 0L,
  stochastic = FALSE,
  noise_phi = 1,
  noise_sigma = 0.1,
  noise_mu = 0
)
```

## Arguments

- a:

  Damping parameter representing mechanical and thermal damping of the
  large-scale flow (default 0.25). Controls the rate at which the
  westerly wind current relaxes toward equilibrium.

- b:

  Advection coefficient for the wave-flow interaction (default 4.0).
  Represents the strength of nonlinear coupling between the westerly
  wind and the planetary waves. Larger values increase wave-mean flow
  interaction.

- F:

  Symmetric thermal forcing representing the equator-to-pole temperature
  gradient (default 8.0). This is the primary control parameter: F \< 4
  gives stable fixed points; F ≈ 6-8 produces periodic or chaotic
  behavior; F \> 8 enhances chaos. Physically represents radiative
  heating.

- G:

  Asymmetric thermal forcing representing land-sea contrast or
  topographic effects (default 1.0). G = 0 gives a symmetric system; G ≠
  0 breaks the y ↔ -y symmetry. Typical values 0-2; G \> 0 favors
  positive y.

- n_steps:

  Number of integration steps (default 500000). With dt = 0.001 and thin
  = 5, this gives 500 time units with 100000 output points.

- dt:

  Time step for RK4 integration (default 0.001). Smaller values give
  smoother trajectories for publication-quality plots.

- x0, y0, z0:

  Initial conditions. Default (1, 0, 0) starts near the westerly
  equilibrium.

- seed:

  Random seed for reproducibility.

- thin:

  Thinning factor for output (default 5). Only every thin-th point is
  stored, reducing memory while maintaining smooth visualization.

- transient:

  Number of initial steps to discard as transient.

- stochastic:

  Logical; if TRUE, add Ornstein-Uhlenbeck noise to the x-equation,
  representing unresolved small-scale atmospheric variability.

- noise_phi:

  OU relaxation rate (default 1.0). Higher values give
  faster-decorrelating noise.

- noise_sigma:

  OU noise amplitude (default 0.1).

- noise_mu:

  OU process mean (default 0.0).

## Value

Object of class "lorenz84_sim" containing:

- t:

  Time vector

- x:

  Westerly wind strength (zonally averaged westerly current)

- y:

  Cosine phase amplitude of planetary wave (north-south displacement)

- z:

  Sine phase amplitude of planetary wave (east-west displacement)

- parameters:

  List of all simulation parameters

## Details

The Lorenz-84 model is a severely truncated spectral model of
mid-latitude atmospheric circulation, representing the interaction
between the zonally averaged westerly jet stream and large-scale
planetary (Rossby) waves:

\$\$\frac{dx}{dt} = -ax - y^2 - z^2 + aF\$\$ \$\$\frac{dy}{dt} = -y +
xy - bxz + G\$\$ \$\$\frac{dz}{dt} = -z + bxy + xz\$\$

**Physical interpretation of state variables:**

- x:

  Intensity of the symmetric, zonally averaged westerly wind current.
  Positive x indicates strong westerlies; x ≈ F at equilibrium without
  waves.

- y, z:

  Amplitudes of a large-scale planetary wave with cosine (y) and
  sine (z) longitudinal phases. The wave amplitude is \\\sqrt{y^2 +
  z^2}\\.

**Parameter regimes:**

- F \< 4:

  Stable fixed point (permanent westerly regime)

- F ≈ 4-6:

  Periodic oscillations (regular seasonal cycle)

- F ≈ 6-8, G ≈ 1:

  Chaotic behavior with intermittent "blocking" events

- F \> 9:

  Strong chaos with rapid fluctuations

The model captures qualitative features of atmospheric variability
including the irregular alternation between zonal (westerly-dominated)
and blocked (wave-amplified) flow patterns. The largest Lyapunov
exponent for standard parameters (F=8, G=1) is approximately \\\lambda_1
\approx 0.07\\.

## References

Lorenz, E.N. (1984). Irregularity: a fundamental property of the
atmosphere. *Tellus A*, 36(2), 98-110.
[doi:10.1111/j.1600-0870.1984.tb00230.x](https://doi.org/10.1111/j.1600-0870.1984.tb00230.x)

Lorenz, E.N. (1990). Can chaos and intransitivity lead to interannual
variability? *Tellus A*, 42(3), 378-389.
[doi:10.1034/j.1600-0870.1990.t01-2-00005.x](https://doi.org/10.1034/j.1600-0870.1990.t01-2-00005.x)

## See also

[`print.lorenz84_sim()`](https://robustecologies.github.io/tuRbulence/reference/print.lorenz84_sim.md),
[`summary.lorenz84_sim()`](https://robustecologies.github.io/tuRbulence/reference/summary.lorenz84_sim.md),
[`plot.lorenz84_sim()`](https://robustecologies.github.io/tuRbulence/reference/plot.lorenz84_sim.md),
[`lorenz_sim()`](https://robustecologies.github.io/tuRbulence/reference/lorenz_sim.md),
[`rossler_sim()`](https://robustecologies.github.io/tuRbulence/reference/rossler_sim.md),
[`simulate_system()`](https://robustecologies.github.io/tuRbulence/reference/simulate_system.md),
[`turbulence_sim()`](https://robustecologies.github.io/tuRbulence/reference/turbulence_sim.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Chaotic Lorenz-84 attractor (standard parameters)
sim <- lorenz84_sim(a = 0.25, b = 4.0, F = 8.0, G = 1.0, n_steps = 100000)

# Inspect results
print(sim)
summary(sim)

# Visualize time series
plot(sim)

# Phase space and 3D attractor
plot(sim, type = "phase")
plot(sim, type = "attractor")

# Explore different forcing regimes
sim_weak <- lorenz84_sim(F = 4.0, G = 0.5, n_steps = 50000)   # Periodic
sim_strong <- lorenz84_sim(F = 10.0, G = 1.0, n_steps = 50000) # Strong chaos
} # }
```
