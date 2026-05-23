# Simulate Hasselmann stochastic climate model

Integrates the extended Hasselmann model describing the climate system
as a slow variable (ocean/ice) forced by fast atmospheric fluctuations,
which appear as stochastic forcing at climate timescales.

## Usage

``` r
hasselmann_sim(
  F_forcing = 0,
  n_steps = 100000L,
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
  transient = 0L
)
```

## Arguments

- F_forcing:

  External radiative forcing anomaly (default 0.0). Represents changes
  in the planetary energy budget from greenhouse gases, solar
  variability, or volcanic aerosols. Positive F represents net warming
  (e.g., increased CO2); negative F represents cooling. Units are
  non-dimensional, but F = 1 corresponds roughly to a forcing of ~3-4
  W/m² (approximately 2×CO2).

- n_steps:

  Number of integration steps.

- dt:

  Time step in non-dimensional units (default 0.1). One time unit
  corresponds to approximately the thermal relaxation time of the
  surface ocean mixed layer (~years to decades). For climate
  applications, dt = 0.1 represents sub-decadal resolution.

- T0:

  Initial surface temperature anomaly (default 0.0). Positive values
  indicate warming relative to a reference climate.

- Td0:

  Initial deep ocean temperature anomaly (default 0.0). The deep ocean
  responds much slower than the surface due to its large heat capacity.

- I0:

  Initial ice extent (default 0.5). Normalized ice coverage where I = 0
  is ice-free and I = 1 is maximum glaciation.

- lambda:

  Climate feedback parameter (default 0.1). The aggregate effect of all
  fast feedbacks (Planck, water vapor, lapse rate, clouds) on surface
  temperature. Higher λ means stronger negative feedback (faster
  relaxation to equilibrium). Climate sensitivity ∝ 1/λ.

- gamma:

  Ocean heat uptake coefficient (default 0.05). Controls the rate of
  heat exchange between surface and deep ocean: γ = C_s/(C_d × τ_m)
  where C_s, C_d are heat capacities and τ_m is the mixing timescale.
  Larger γ means faster equilibration between surface and deep ocean.

- kappa:

  Deep ocean relaxation rate (default 0.01). Controls how quickly the
  deep ocean temperature relaxes toward the surface temperature: κ ∝
  1/τ_deep where τ_deep is the deep ocean turnover time (~centuries).
  Small κ reflects the long memory of the deep ocean.

- alpha_ice:

  Ice-albedo feedback strength (default 0.1). How strongly ice extent
  affects surface temperature through changes in planetary albedo. The
  ice-albedo feedback is positive: more ice → higher albedo → cooling →
  more ice. This is captured as −α×I in the temperature equation.

- beta_ice:

  Ice dynamics timescale (default 0.02). The rate at which ice extent
  adjusts toward its equilibrium value: β ∝ 1/τ_ice where τ_ice is the
  ice sheet response time (millennia for continental ice sheets, decades
  for sea ice).

- mu_ice:

  Ice-temperature sensitivity (default 0.5). Slope of the ice
  equilibrium curve: I_eq = max(0, −μ×T). Larger μ means ice responds
  more sensitively to temperature changes. This parameter controls the
  strength of the ice-albedo feedback loop.

- sigma_T:

  Temperature noise amplitude (default 0.3). Represents the integrated
  effect of high-frequency atmospheric variability (weather) on the slow
  climate state. Hasselmann's key insight: weather appears as stochastic
  forcing to the climate system.

- sigma_I:

  Ice noise amplitude (default 0.1). Represents stochastic variability
  in ice dynamics from factors like precipitation variability, ocean
  heat transport fluctuations, or ice calving events.

- ou_phi_T, ou_phi_I:

  Ornstein-Uhlenbeck relaxation rates for temperature and ice noise
  (defaults 0.0). When φ = 0, white noise is used; φ \> 0 produces
  temporally correlated (red) noise with decorrelation time 1/φ. Red
  noise better represents persistent atmospheric patterns.

- ou_mu_T, ou_mu_I:

  Ornstein-Uhlenbeck process means for temperature and ice noise
  (defaults 0.0). The OU process reverts to these values with rate φ.
  Non-zero means introduce systematic bias representing long-term
  forcing.

- stochastic:

  Logical; if TRUE (default), enable stochastic forcing. The
  deterministic version shows the underlying attractor structure.

- seed:

  Random seed for reproducibility.

- thin:

  Thinning factor for output (store every thin-th point).

- transient:

  Number of initial steps to discard as transient.

## Value

Object of class "hasselmann_sim" containing:

- t:

  Time vector (non-dimensional, ~decades per unit)

- T:

  Surface temperature anomaly relative to reference climate

- Td:

  Deep ocean temperature anomaly

- I:

  Ice extent (0 = ice-free, 1 = maximum glaciation)

- parameters:

  List of all simulation parameters

## Details

The extended Hasselmann model equations are: \$\$\frac{dT}{dt} =
-\lambda T - \gamma(T - T_d) - \alpha I + F + \sigma_T \xi_T\$\$
\$\$\frac{dT_d}{dt} = -\kappa(T_d - T)\$\$ \$\$\frac{dI}{dt} =
-\beta(I - I\_{eq}(T)) + \sigma_I \xi_I\$\$

where T is the surface temperature anomaly, Td is the deep ocean
temperature, I is ice extent, and F is external forcing. The ice
equilibrium function \\I\_{eq}(T) = \max(0, -\mu T)\\ models the
ice-albedo feedback.

The model produces red noise spectra characteristic of climate
variability, with the deep ocean acting as a low-pass filter. The ice
component can generate glacial-interglacial oscillations through the
albedo feedback.

This is inherently a stochastic model; the deterministic version
(stochastic=FALSE) shows the underlying attractor structure without
atmospheric noise.

## References

Hasselmann, K. (1976). Stochastic climate models. Part I: Theory.
*Tellus*, 28, 473-485.
[doi:10.3402/tellusa.v28i6.11316](https://doi.org/10.3402/tellusa.v28i6.11316)

Saltzman, B. (2002). *Dynamical paleoclimatology: Generalized theory of
global climate change*. Academic Press. ISBN: 978-0126173314.

## See also

[`print.hasselmann_sim()`](https://robustecologies.github.io/tuRbulence/reference/print.hasselmann_sim.md),
[`summary.hasselmann_sim()`](https://robustecologies.github.io/tuRbulence/reference/summary.hasselmann_sim.md),
[`plot.hasselmann_sim()`](https://robustecologies.github.io/tuRbulence/reference/plot.hasselmann_sim.md),
[`hasselmann_batch()`](https://robustecologies.github.io/tuRbulence/reference/hasselmann_batch.md),
[`simulate_system()`](https://robustecologies.github.io/tuRbulence/reference/simulate_system.md),
[`turbulence_sim()`](https://robustecologies.github.io/tuRbulence/reference/turbulence_sim.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate stochastic climate variability
sim <- hasselmann_sim(F_forcing = 0, n_steps = 200000, seed = 42)

# Inspect results
print(sim)
summary(sim)

# Visualize climate dynamics
plot(sim)
plot(sim, type = "phase")
plot(sim, type = "spectrum")
plot(sim, type = "ice_temp")

# Forced warming scenario
sim_warm <- hasselmann_sim(F_forcing = 0.5, n_steps = 200000)
plot(sim_warm)
} # }
```
