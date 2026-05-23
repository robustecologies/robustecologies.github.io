# Simulate Charney-DeVore atmospheric blocking model

Integrates the Charney-DeVore truncated spectral model describing the
interaction between zonal flow and planetary waves, with optional
stochastic forcing representing synoptic-scale variability.

## Usage

``` r
charney_devore_sim(
  F = 1.5,
  n_steps = 100000L,
  dt = 0.01,
  x0 = 1,
  y0 = 0.1,
  z0 = 0.1,
  k = 0.1,
  alpha = 1,
  beta = 0.5,
  delta = 1,
  sigma = 0.15,
  ou_phi = 0,
  ou_mu = 0,
  stochastic = TRUE,
  seed = NULL,
  thin = 1L,
  transient = 0L
)
```

## Arguments

- F:

  Non-dimensional external thermal forcing (default 1.5). Represents the
  meridional (equator-to-pole) temperature gradient that drives the
  zonal wind through thermal wind balance: \\F \propto \partial T /
  \partial y\\. This is the primary control parameter. F \< 1 gives
  stable zonal flow; F ≈ 1-2 produces multiple equilibria (bistability);
  F \> 2.5 tends toward the high-index state only.

- n_steps:

  Number of integration steps.

- dt:

  Time step for Euler-Maruyama integration.

- x0, y0, z0:

  Initial conditions for state variables. x0 = 1 starts near the
  high-index (zonal) state; smaller x0 with larger y0, z0 starts closer
  to blocking.

- k:

  Ekman (frictional) damping coefficient (default 0.1). Represents the
  rate of momentum dissipation through boundary layer friction: \\k \sim
  f / (2H)\\ where f is Coriolis parameter and H is the Ekman depth.
  Higher k stabilizes the system and favors the zonal state.

- alpha:

  Nonlinear wave-mean flow interaction coefficient (default 1.0).
  Controls the strength of the triad coupling between zonal flow and
  waves. Derived from the vorticity equation as \\\alpha \sim k_w^2\\
  where k_w is the wavenumber. Larger α enhances nonlinear interactions
  and can destabilize the zonal flow.

- beta:

  Topographic/thermal wave forcing coefficient (default 0.5). Represents
  how the external thermal gradient forces planetary waves through
  baroclinic conversion. In the original model, β is related to
  topographic height and wavenumber structure. Larger β means stronger
  wave excitation by the mean flow gradient.

- delta:

  Wave-wave (self) interaction coefficient (default 1.0). Controls the
  internal nonlinear coupling between wave modes y and z. Derived from
  the spectral projection of the advection term. δ affects the wave
  amplitude saturation mechanism.

- sigma:

  Noise amplitude (default 0.15). Represents unresolved synoptic-scale
  variability (weather) that perturbs the large-scale flow. Noise can
  trigger transitions between metastable states.

- ou_phi:

  Ornstein-Uhlenbeck relaxation rate (default 0.0). When φ = 0, white
  noise is used; φ \> 0 produces temporally correlated noise with
  decorrelation time 1/φ, representing persistent weather patterns.

- ou_mu:

  Ornstein-Uhlenbeck process mean (default 0.0). The OU process reverts
  to this value with rate φ. Non-zero μ introduces systematic bias in
  the noise forcing, which can shift the preferred state occupancy.

- stochastic:

  Logical; if TRUE (default), enable stochastic forcing.

- seed:

  Random seed for reproducibility.

- thin:

  Thinning factor for output (store every thin-th point).

- transient:

  Number of initial steps to discard as transient.

## Value

Object of class "charney_devore_sim" containing:

- t:

  Time vector (non-dimensional, scaled by 1/k)

- x:

  Zonal flow amplitude (strength of the westerly jet). High x
  corresponds to strong, straight jet stream (high-index pattern).

- y:

  Amplitude of planetary wave (cosine phase). Represents the north-south
  displacement of the jet stream.

- z:

  Amplitude of planetary wave (sine phase). Represents the east-west
  displacement of wave phase.

- wave_energy:

  Total wave energy \\y^2 + z^2\\, measuring the strength of the
  blocking pattern.

- parameters:

  List of all simulation parameters

## Details

The simplified Charney-DeVore equations are: \$\$\frac{dx}{dt} = k(F -
x) - \alpha yz + \beta y\$\$ \$\$\frac{dy}{dt} = -ky + \alpha xz - \beta
x - \delta z\$\$ \$\$\frac{dz}{dt} = -kz + \delta y\$\$

where x represents the zonal flow amplitude that relaxes toward the
external forcing F, and (y, z) represent the amplitudes of planetary
wave modes. The parameter F controls the strength of the meridional
temperature gradient (thermal forcing) and determines the equilibrium
zonal flow in the absence of wave activity.

This simplified 3-mode version captures the essential wave-mean flow
interaction. The full 6-mode Charney-DeVore model exhibits bistability
(multiple equilibria) under specific parameter regimes:

- High-index state: Strong zonal flow (large x), weak waves

- Low-index state (blocking): Weak zonal flow, amplified waves

Stochastic forcing induces variability around the equilibrium and can
create blocking-like events when wave amplitudes exceed threshold
values, characteristic of mid-latitude weather patterns.

## References

Charney, J.G. & DeVore, J.G. (1979). Multiple flow equilibria in the
atmosphere and blocking. *Journal of the Atmospheric Sciences*, 36,
1205-1216.
[doi:10.1175/1520-0469(1979)036\<1205:MFEITA\>2.0.CO;2](https://doi.org/10.1175/1520-0469%281979%29036%3C1205%3AMFEITA%3E2.0.CO%3B2)

## See also

[`print.charney_devore_sim()`](https://robustecologies.github.io/tuRbulence/reference/print.charney_devore_sim.md),
[`summary.charney_devore_sim()`](https://robustecologies.github.io/tuRbulence/reference/summary.charney_devore_sim.md),
[`plot.charney_devore_sim()`](https://robustecologies.github.io/tuRbulence/reference/plot.charney_devore_sim.md),
[`charney_devore_6mode()`](https://robustecologies.github.io/tuRbulence/reference/charney_devore_6mode.md),
[`charney_devore_batch()`](https://robustecologies.github.io/tuRbulence/reference/charney_devore_batch.md),
[`simulate_system()`](https://robustecologies.github.io/tuRbulence/reference/simulate_system.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate atmospheric blocking dynamics
sim <- charney_devore_sim(F = 1.5, n_steps = 200000, seed = 42)

# Inspect results
print(sim)
summary(sim)

# Visualize dynamics
plot(sim)
plot(sim, type = "phase")
plot(sim, type = "attractor")
plot(sim, type = "blocking")

# Deterministic comparison
sim_det <- charney_devore_sim(F = 1.5, stochastic = FALSE, n_steps = 100000)
plot(sim_det, type = "phase")
} # }
```
