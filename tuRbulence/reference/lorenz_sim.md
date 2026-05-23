# Simulate the Lorenz system

Integrates the classic Lorenz equations for atmospheric convection using
4th-order Runge-Kutta, with optional Ornstein-Uhlenbeck stochastic
forcing.

## Usage

``` r
lorenz_sim(
  sigma = 10,
  rho = 28,
  beta = 8/3,
  n_steps = 500000L,
  dt = 0.001,
  x0 = 1,
  y0 = 1,
  z0 = 1,
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

- sigma:

  Prandtl number (default 10), the ratio of kinematic viscosity to
  thermal diffusivity: \\\sigma = \nu / \kappa\\. Controls how quickly
  momentum diffuses relative to heat. Higher σ increases the coupling
  between velocity and temperature fields.

- rho:

  Rayleigh number (default 28), the ratio of buoyancy to viscous forces:
  \\\rho = g \alpha \Delta T H^3 / (\nu \kappa)\\. This is the primary
  control parameter. The system becomes chaotic for ρ \> 24.74 (the
  critical value for σ = 10, β = 8/3). Higher ρ means stronger
  convective driving.

- beta:

  Geometric factor (default 8/3), related to the aspect ratio of
  convection cells: \\\beta = 4 / (1 + a^2)\\ where a is the horizontal
  wave number. Controls the rate of energy dissipation in the vertical
  mode.

- n_steps:

  Number of integration steps (default 500000). With dt = 0.001 and thin
  = 5, this gives 500 time units with 100000 output points.

- dt:

  Time step for RK4 integration (default 0.001). Smaller values give
  smoother trajectories; 0.001 is recommended for publication-quality
  plots.

- x0, y0, z0:

  Initial conditions. Defaults (1, 1, 1) are near the unstable fixed
  point to allow transient dynamics toward the attractor.

- seed:

  Random seed for reproducibility.

- thin:

  Thinning factor for output (default 5). Only every thin-th point is
  stored, reducing memory while maintaining smooth visualization.

- transient:

  Number of initial steps to discard as transient.

- stochastic:

  Logical; if TRUE, add Ornstein-Uhlenbeck noise to the x-equation,
  representing unresolved turbulent fluctuations.

- noise_phi:

  OU relaxation rate (default 1.0). Higher values give
  faster-decorrelating noise.

- noise_sigma:

  OU noise amplitude (default 0.1).

- noise_mu:

  OU process mean (default 0.0).

## Value

Object of class "lorenz_sim" containing:

- t:

  Time vector

- x:

  Convective intensity (proportional to circulatory fluid velocity)

- y:

  Horizontal temperature gradient (temperature difference between
  ascending and descending currents)

- z:

  Vertical temperature stratification (deviation from linear vertical
  temperature profile)

- parameters:

  List of all simulation parameters

## Details

The Lorenz system is a simplified model of Rayleigh-Bénard convection,
describing the motion of a fluid layer heated from below:

\$\$\frac{dx}{dt} = \sigma(y - x)\$\$ \$\$\frac{dy}{dt} = x(\rho - z) -
y\$\$ \$\$\frac{dz}{dt} = xy - \beta z\$\$

**Physical interpretation of state variables:**

- x:

  Rate of convective overturning (positive = clockwise circulation)

- y:

  Horizontal temperature variation (temperature contrast across cell)

- z:

  Departure of vertical temperature profile from linearity (positive =
  more stable stratification)

**Parameter regimes:**

- ρ \< 1:

  No convection (heat conduction only)

- 1 \< ρ \< 24.74:

  Stable convection (steady rolls)

- ρ \> 24.74:

  Chaotic convection (butterfly attractor)

For standard parameters (σ = 10, ρ = 28, β = 8/3), the system exhibits
chaotic behavior on a strange attractor. The trajectory alternates
unpredictably between two lobes representing clockwise and
counterclockwise convection. The Lyapunov spectrum is approximately
\\\lambda_1 \approx 0.906\\, \\\lambda_2 = 0\\, \\\lambda_3 \approx
-14.57\\, giving a Kaplan-Yorke dimension of approximately 2.06.

## References

Lorenz, E.N. (1963). Deterministic nonperiodic flow. *Journal of the
Atmospheric Sciences*, 20(2), 130-141.
[doi:10.1175/1520-0469(1963)020\<0130:DNF\>2.0.CO;2](https://doi.org/10.1175/1520-0469%281963%29020%3C0130%3ADNF%3E2.0.CO%3B2)

Sparrow, C. (1982). *The Lorenz Equations: Bifurcations, Chaos, and
Strange Attractors*. Springer-Verlag.
[doi:10.1007/978-1-4612-5767-7](https://doi.org/10.1007/978-1-4612-5767-7)

## See also

[`print.lorenz_sim()`](https://robustecologies.github.io/tuRbulence/reference/print.lorenz_sim.md),
[`summary.lorenz_sim()`](https://robustecologies.github.io/tuRbulence/reference/summary.lorenz_sim.md),
[`plot.lorenz_sim()`](https://robustecologies.github.io/tuRbulence/reference/plot.lorenz_sim.md),
[`rossler_sim()`](https://robustecologies.github.io/tuRbulence/reference/rossler_sim.md),
[`lorenz84_sim()`](https://robustecologies.github.io/tuRbulence/reference/lorenz84_sim.md),
[`simulate_system()`](https://robustecologies.github.io/tuRbulence/reference/simulate_system.md),
[`turbulence_sim()`](https://robustecologies.github.io/tuRbulence/reference/turbulence_sim.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Classic chaotic Lorenz attractor
sim <- lorenz_sim(sigma = 10, rho = 28, beta = 8/3, n_steps = 100000)

# Inspect results
print(sim)
summary(sim)

# Visualize time series and attractor
plot(sim)
plot(sim, type = "phase")
plot(sim, type = "attractor")

# Explore subcritical regime (stable convection)
sim_stable <- lorenz_sim(rho = 20, n_steps = 50000)
plot(sim_stable, type = "attractor")  # Converges to fixed point
} # }
```
