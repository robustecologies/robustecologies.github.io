# Simulate Stommel box model for thermohaline circulation

Integrates the Stommel two-box model describing the competition between
thermally-driven and salinity-driven ocean circulation, with optional
stochastic forcing representing atmospheric variability.

## Usage

``` r
stommel_sim(
  eta2 = 1,
  n_steps = 100000L,
  dt = 0.01,
  T0 = 2,
  S0 = 1,
  eta1 = 3,
  eta3 = 0.3,
  sigma_T = 0.2,
  sigma_S = 0.2,
  ou_phi = 0,
  ou_mu_T = 0,
  ou_mu_S = 0,
  stochastic = TRUE,
  seed = NULL,
  thin = 1L,
  transient = 0L
)
```

## Arguments

- eta2:

  Non-dimensional freshwater flux at high latitudes (default 1.0). This
  is the primary bifurcation control parameter, representing the rate of
  freshwater input from precipitation minus evaporation in the polar
  box. \\\eta_2 = E / (\alpha \Delta T^\*)\\ where E is freshwater flux
  and \\\Delta T^\*\\ is the reference temperature difference. Values
  0.8-1.5 span the bistable regime where both thermal and haline
  circulation modes are stable.

- n_steps:

  Number of integration steps.

- dt:

  Time step for Euler-Maruyama integration.

- T0, S0:

  Initial conditions for non-dimensional temperature and salinity
  differences between equatorial and polar boxes. T0 = 2, S0 = 1 starts
  in the thermal circulation regime.

- eta1:

  Non-dimensional surface thermal forcing (default 3.0). Represents the
  restoring temperature difference imposed by atmosphere-ocean heat
  flux: \\\eta_1 = T^\*\_{eq} / \Delta T^\*\\ where \\T^\*\_{eq}\\ is
  the equilibrium equator-pole temperature difference. Larger η₁ favors
  the thermal mode.

- eta3:

  Ratio of thermal to haline relaxation timescales (default 0.3).
  Physically, \\\eta_3 = \tau_T / \tau_S\\ where τ_T and τ_S are the
  e-folding times for surface temperature and salinity anomalies. Since
  heat exchanges faster than salt (τ_T \< τ_S), η₃ \< 1 in realistic
  oceans. Smaller η₃ makes salinity more persistent and favors haline
  circulation.

- sigma_T, sigma_S:

  Noise amplitudes for temperature and salinity (defaults 0.2).
  Represent unresolved atmospheric variability affecting air-sea fluxes.
  Typical ratio σ_S/σ_T depends on whether precipitation or temperature
  variability dominates.

- ou_phi:

  Ornstein-Uhlenbeck relaxation rate (default 0.0). When φ = 0, white
  noise is used; φ \> 0 produces temporally correlated (red) noise with
  decorrelation time 1/φ.

- ou_mu_T, ou_mu_S:

  Ornstein-Uhlenbeck process means for temperature and salinity noise
  (defaults 0.0). The OU process reverts to these values with rate φ.
  Non-zero means introduce systematic bias in the noise forcing.

- stochastic:

  Logical; if TRUE (default), enable stochastic forcing.

- seed:

  Random seed for reproducibility.

- thin:

  Thinning factor for output (store every thin-th point).

- transient:

  Number of initial steps to discard as transient.

## Value

Object of class "stommel_sim" containing:

- t:

  Time vector (non-dimensional, scaled by τ_T)

- T:

  Temperature difference (equator minus pole), positive = normal
  gradient

- S:

  Salinity difference (equator minus pole)

- q:

  Flow strength q = T - S, proportional to density difference. Positive
  q indicates thermally-driven circulation (like modern Atlantic);
  negative q indicates salinity-driven (reversed) circulation.

- parameters:

  List of all simulation parameters

## Details

The non-dimensional Stommel equations are: \$\$\frac{dT}{dt} = \eta_1 -
T(1 + \|T - S\|)\$\$ \$\$\frac{dS}{dt} = \eta_2 - S(\eta_3 + \|T -
S\|)\$\$

With stochastic forcing: \$\$dT = \[\eta_1 - T(1 + \|q\|)\] dt +
\sigma_T dW_T\$\$ \$\$dS = \[\eta_2 - S(\eta_3 + \|q\|)\] dt + \sigma_S
dW_S\$\$

where \\q = T - S\\ is the flow strength (proportional to density
difference).

The model exhibits bistability: for certain parameter values, both a
thermally-dominated circulation (\\q \> 0\\, like the modern Atlantic)
and a salinity-dominated circulation (\\q \< 0\\) are stable. Stochastic
forcing can induce transitions between these states, creating
intermittent dynamics reminiscent of paleoclimate records.

## References

Stommel, H. (1961). Thermohaline convection with two stable regimes of
flow. *Tellus*, 13(2), 224-230.
[doi:10.3402/tellusa.v13i2.9491](https://doi.org/10.3402/tellusa.v13i2.9491)

Cessi, P. (1994). A simple box model of stochastically forced
thermohaline flow. *Journal of Physical Oceanography*, 24(9), 1911-1920.
[doi:10.1175/1520-0485(1994)024\<1911:ASBMOS\>2.0.CO;2](https://doi.org/10.1175/1520-0485%281994%29024%3C1911%3AASBMOS%3E2.0.CO%3B2)

## See also

[`print.stommel_sim()`](https://robustecologies.github.io/tuRbulence/reference/print.stommel_sim.md),
[`summary.stommel_sim()`](https://robustecologies.github.io/tuRbulence/reference/summary.stommel_sim.md),
[`plot.stommel_sim()`](https://robustecologies.github.io/tuRbulence/reference/plot.stommel_sim.md),
[`stommel_batch()`](https://robustecologies.github.io/tuRbulence/reference/stommel_batch.md),
[`simulate_system()`](https://robustecologies.github.io/tuRbulence/reference/simulate_system.md),
[`turbulence_sim()`](https://robustecologies.github.io/tuRbulence/reference/turbulence_sim.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate bistable thermohaline circulation
sim <- stommel_sim(eta2 = 1.0, n_steps = 200000, seed = 42)

# Inspect results
print(sim)
summary(sim)

# Visualize flow dynamics
plot(sim)
plot(sim, type = "phase")
plot(sim, type = "flow")
plot(sim, type = "density")

# Deterministic comparison
sim_det <- stommel_sim(eta2 = 1.0, stochastic = FALSE, n_steps = 100000)
plot(sim_det)
} # }
```
