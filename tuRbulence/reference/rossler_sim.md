# Simulate the Rössler system

Integrates the Rössler attractor equations using 4th-order Runge-Kutta,
with optional Ornstein-Uhlenbeck stochastic forcing.

## Usage

``` r
rossler_sim(
  a = 0.2,
  b = 0.2,
  c = 5.7,
  n_steps = 500000L,
  dt = 0.001,
  x0 = 1,
  y0 = 1,
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

  Linear feedback parameter controlling the spiraling behavior (default
  0.2). This parameter determines the rate at which trajectories spiral
  outward in the x-y plane. Increasing a above ~0.3 begins a
  period-doubling cascade; a ≈ 0.1 gives a simple limit cycle.

- b:

  Offset parameter for the z-dynamics (default 0.2). Controls the
  baseline excitation of the z variable. In chemical oscillator
  analogies, b represents a constant feed rate of reactant. Small
  changes in b have limited effect on dynamics.

- c:

  Threshold parameter controlling the reinjection mechanism (default
  5.7). When z exceeds approximately c, the z(x - c) term becomes
  strongly positive, producing a rapid excursion that reinjects the
  trajectory back toward the attractor. This is the primary bifurcation
  parameter: c \< 4 gives periodic orbits; c ≈ 4-6 shows
  period-doubling; c \> 5.7 produces fully developed chaos; c \> 12 can
  produce hyperchaos or divergence.

- n_steps:

  Number of integration steps (default 500000). With dt = 0.001 and thin
  = 5, this gives 500 time units (~80 orbits) with 100000 output points.

- dt:

  Time step for RK4 integration (default 0.001). Smaller values give
  smoother trajectories; the Rössler system is less stiff than Lorenz.

- x0, y0, z0:

  Initial conditions. Defaults (1, 1, 0) place the trajectory near the
  inner loop of the attractor.

- seed:

  Random seed for reproducibility.

- thin:

  Thinning factor for output (default 5). Only every thin-th point is
  stored, reducing memory while maintaining smooth visualization.

- transient:

  Number of initial steps to discard as transient.

- stochastic:

  Logical; if TRUE, add Ornstein-Uhlenbeck noise to the x-equation,
  representing external perturbations or model uncertainty.

- noise_phi:

  OU relaxation rate (default 1.0). Higher values give
  faster-decorrelating noise.

- noise_sigma:

  OU noise amplitude (default 0.1).

- noise_mu:

  OU process mean (default 0.0).

## Value

Object of class "rossler_sim" containing:

- t:

  Time vector

- x:

  First state variable (oscillatory component)

- y:

  Second state variable (oscillatory component, phase-shifted from x)

- z:

  Third state variable (slow pulsing component, the "funnel")

- parameters:

  List of all simulation parameters

## Details

The Rössler system was designed as a minimal chaotic flow, simpler than
the Lorenz system but exhibiting similar strange attractor dynamics:

\$\$\frac{dx}{dt} = -y - z\$\$ \$\$\frac{dy}{dt} = x + ay\$\$
\$\$\frac{dz}{dt} = b + z(x - c)\$\$

Unlike Lorenz which has two symmetric lobes, the Rössler attractor has a
single-band structure: trajectories spiral outward in the x-y plane
until z triggers a reinjection pulse, creating the characteristic
"funnel" shape.

**Physical interpretation:** While purely mathematical in origin, the
equations can model certain chemical reaction kinetics (e.g., Belousov-
Zhabotinsky reaction variants) where x and y represent oscillating
concentrations and z represents a slowly accumulating inhibitor that
periodically "fires" to reset the oscillation.

**Parameter regimes:**

- c \< 3:

  Simple limit cycle (periodic)

- c ≈ 3-4:

  Period-2 orbit

- c ≈ 4-5:

  Period-4 and higher period orbits (period-doubling cascade)

- c ≈ 5.7:

  Onset of chaos, standard chaotic parameters

- c \> 12:

  Risk of unbounded growth or hyperchaos

For standard parameters (a = 0.2, b = 0.2, c = 5.7), the system exhibits
a strange attractor with Lyapunov spectrum approximately \\\lambda_1
\approx 0.07\\, \\\lambda_2 = 0\\, \\\lambda_3 \approx -5.4\\, giving a
Kaplan-Yorke dimension of approximately 2.01, slightly lower than the
Lorenz attractor.

## References

Rössler, O.E. (1976). An equation for continuous chaos. *Physics Letters
A*, 57(5), 397-398.
[doi:10.1016/0375-9601(76)90101-8](https://doi.org/10.1016/0375-9601%2876%2990101-8)

Rössler, O.E. (1979). An equation for hyperchaos. *Physics Letters A*,
71(2-3), 155-157.
[doi:10.1016/0375-9601(79)90150-6](https://doi.org/10.1016/0375-9601%2879%2990150-6)

## See also

[`print.rossler_sim()`](https://robustecologies.github.io/tuRbulence/reference/print.rossler_sim.md),
[`summary.rossler_sim()`](https://robustecologies.github.io/tuRbulence/reference/summary.rossler_sim.md),
[`plot.rossler_sim()`](https://robustecologies.github.io/tuRbulence/reference/plot.rossler_sim.md),
[`lorenz_sim()`](https://robustecologies.github.io/tuRbulence/reference/lorenz_sim.md),
[`lorenz84_sim()`](https://robustecologies.github.io/tuRbulence/reference/lorenz84_sim.md),
[`simulate_system()`](https://robustecologies.github.io/tuRbulence/reference/simulate_system.md),
[`turbulence_sim()`](https://robustecologies.github.io/tuRbulence/reference/turbulence_sim.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Classic Rössler attractor (chaotic regime)
sim <- rossler_sim(a = 0.2, b = 0.2, c = 5.7, n_steps = 100000)

# Inspect results
print(sim)
summary(sim)

# Visualize time series and attractor
plot(sim)
plot(sim, type = "phase")
plot(sim, type = "attractor")

# Explore period-doubling route to chaos
sim_periodic <- rossler_sim(c = 3.0, n_steps = 50000)   # Simple limit cycle
sim_p2 <- rossler_sim(c = 4.0, n_steps = 50000)         # Period-2
sim_chaos <- rossler_sim(c = 5.7, n_steps = 50000)      # Chaotic

plot(sim_periodic, type = "attractor")
plot(sim_p2, type = "attractor")
plot(sim_chaos, type = "attractor")
} # }
```
