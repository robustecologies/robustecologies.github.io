# Simulate the full 6-mode Charney-DeVore model with bistability

Integrates the complete 6-mode Charney-DeVore spectral model that
exhibits genuine bistability between high-index (zonal) and low-index
(blocked) atmospheric states. This is the full model from the original
1979 paper, preserving all six interacting spectral modes.

## Usage

``` r
charney_devore_6mode(
  x1_star = 0.95,
  x4_star = -0.76095,
  n_steps = 100000L,
  dt = 0.01,
  x1_0 = 0.95,
  x2_0 = 0,
  x3_0 = 0,
  x4_0 = -0.76,
  x5_0 = 0,
  x6_0 = 0,
  C = 0.1,
  b = 0.5,
  beta = 1.25,
  gamma = 0.2,
  sigma = 0,
  ou_phi = 0,
  ou_mu = 0,
  stochastic = FALSE,
  seed = NULL,
  thin = 1L,
  transient = 0L
)
```

## Arguments

- x1_star, x4_star:

  Forced equilibrium states for zonal modes (defaults 0.95 and
  -0.76095). These represent the externally forced zonal flow that the
  system relaxes toward in the absence of wave activity. Higher x1_star
  represents stronger thermal forcing.

- n_steps:

  Number of integration steps.

- dt:

  Time step for Euler-Maruyama integration.

- x1_0, x2_0, x3_0, x4_0, x5_0, x6_0:

  Initial conditions for the six spectral modes. For high-index initial
  condition use x1_0 near x1_star with small wave amplitudes; for
  low-index use smaller x1_0 with larger wave modes.

- C:

  Thermal relaxation rate / Ekman damping (default 0.1). Corresponds to
  a damping time of approximately 10 days. Higher C means faster
  damping.

- b:

  Channel half-width parameter (default 0.5). Controls the aspect ratio
  of the beta-plane channel geometry. Standard value gives a channel of
  approximately 6300 km x 1600 km.

- beta:

  Planetary vorticity gradient (default 1.25). Related to the variation
  of the Coriolis parameter with latitude. Standard value corresponds to
  a central latitude of 45 degrees.

- gamma:

  Topographic amplitude (default 0.2). Represents the height of bottom
  topography in non-dimensional units. Standard value corresponds to
  approximately 200 m amplitude.

- sigma:

  Noise amplitude (default 0.0). Amplitude of stochastic forcing
  representing unresolved atmospheric variability.

- ou_phi:

  Ornstein-Uhlenbeck relaxation rate (default 0.0). When phi = 0, white
  noise is used; phi \> 0 produces temporally correlated (red) noise.

- ou_mu:

  Ornstein-Uhlenbeck process mean (default 0.0).

- stochastic:

  Logical; if TRUE, enable stochastic forcing.

- seed:

  Random seed for reproducibility.

- thin:

  Thinning factor for output (store every thin-th point).

- transient:

  Number of initial steps to discard as transient.

## Value

Object of class "charney_devore_6mode" containing:

- t:

  Time vector (non-dimensional)

- x1:

  First zonal flow mode (primary indicator of regime)

- x2,x3:

  First wave modes (wavenumber m=1)

- x4:

  Second zonal flow mode

- x5,x6:

  Second wave modes (wavenumber m=2)

- zonal_index:

  Combined zonal index (x1 + x4)

- wave_energy:

  Total wave energy (x2^2 + x3^2 + x5^2 + x6^2)

- parameters:

  List of all model parameters including derived coefficients

## Details

The full 6-mode Charney-DeVore equations are: \$\$\dot{x}\_1 = -C(x_1 -
x_1^\*) + \tilde{\gamma}\_1 x_3\$\$ \$\$\dot{x}\_2 = -Cx_2 + \beta_1
x_3 - \alpha_1 x_1 x_3 - \delta_1 x_4 x_6\$\$ \$\$\dot{x}\_3 = -Cx_3 -
\beta_1 x_2 + \alpha_1 x_1 x_2 + \delta_1 x_4 x_5 - \gamma_1 x_1\$\$
\$\$\dot{x}\_4 = -C(x_4 - x_4^\*) + \varepsilon(x_2 x_6 - x_3 x_5) +
\tilde{\gamma}\_2 x_6\$\$ \$\$\dot{x}\_5 = -Cx_5 + \beta_2 x_6 -
\alpha_2 x_1 x_6 - \delta_2 x_4 x_3\$\$ \$\$\dot{x}\_6 = -Cx_6 - \beta_2
x_5 + \alpha_2 x_1 x_5 + \delta_2 x_4 x_2 - \gamma_2 x_4\$\$

The derived coefficients (alpha, beta, gamma, delta, epsilon) depend on
the channel geometry (b), planetary vorticity gradient (beta), and
topographic amplitude (gamma) according to the spectral truncation
formulas from Charney and DeVore (1979).

Unlike the simplified 3-mode version, this full model can exhibit
genuine bistability under appropriate parameter regimes. The original
Charney-DeVore (1979) paper demonstrated that for specific topographic
amplitudes and forcing strengths, two stable equilibria coexist: the
high-index (strong zonal flow, weak waves) and low-index (weak zonal
flow, amplified waves, blocking) states. Finding these bistable regimes
requires careful parameter tuning or bifurcation analysis; the standard
parameters produce a single stable equilibrium. Stochastic forcing can
induce variability and intermittent transitions even without
deterministic bistability.

## References

Charney, J.G. & DeVore, J.G. (1979). Multiple flow equilibria in the
atmosphere and blocking. *Journal of the Atmospheric Sciences*, 36,
1205-1216.
[doi:10.1175/1520-0469(1979)036\<1205:MFEITA\>2.0.CO;2](https://doi.org/10.1175/1520-0469%281979%29036%3C1205%3AMFEITA%3E2.0.CO%3B2)

Kwasniok, F. (2023). On the interaction of stochastic forcing and regime
dynamics. *Nonlinear Processes in Geophysics*, 30, 49-62.
[doi:10.5194/npg-30-49-2023](https://doi.org/10.5194/npg-30-49-2023)

## See also

[`print.charney_devore_6mode()`](https://robustecologies.github.io/tuRbulence/reference/print.charney_devore_6mode.md),
[`summary.charney_devore_6mode()`](https://robustecologies.github.io/tuRbulence/reference/summary.charney_devore_6mode.md),
[`plot.charney_devore_6mode()`](https://robustecologies.github.io/tuRbulence/reference/plot.charney_devore_6mode.md),
[`charney_devore_sim()`](https://robustecologies.github.io/tuRbulence/reference/charney_devore_sim.md),
[`charney_devore_batch()`](https://robustecologies.github.io/tuRbulence/reference/charney_devore_batch.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# High-index initial condition (near zonal equilibrium)
sim_hi <- charney_devore_6mode(
    x1_0 = 0.95, x2_0 = 0.0, x3_0 = 0.0,
    x4_0 = -0.76, x5_0 = 0.0, x6_0 = 0.0,
    n_steps = 100000, stochastic = FALSE
)
print(sim_hi)

# Low-index initial condition (near blocked equilibrium)
sim_lo <- charney_devore_6mode(
    x1_0 = 0.2, x2_0 = 0.3, x3_0 = 0.3,
    x4_0 = -0.5, x5_0 = 0.1, x6_0 = 0.1,
    n_steps = 100000, stochastic = FALSE
)
print(sim_lo)

# Compare final states - should show bistability
cat("High-index final x1:", tail(sim_hi$x1, 1), "\n")
cat("Low-index final x1:", tail(sim_lo$x1, 1), "\n")

# Stochastic simulation with regime transitions
sim_stoch <- charney_devore_6mode(
    n_steps = 500000, sigma = 0.05, stochastic = TRUE, seed = 42
)
plot(sim_stoch)
} # }
```
