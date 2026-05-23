# Simulate stochastic Duffing oscillator for von Kármán flow

Integrates the stochastic Duffing equations with Ornstein-Uhlenbeck
forcing that model the turbulent attractor dynamics in von Kármán
swirling flow.

## Usage

``` r
vonkarman_sim(
  mu,
  n_steps = 100000L,
  dt = 0.01,
  x0 = 0.1,
  y0 = 0,
  z0 = NULL,
  a = 0.2,
  phi = 0.9,
  sigma = 0.2,
  omega = 1,
  ou_mu = NULL,
  seed = NULL,
  thin = 1L,
  transient = 0L
)
```

## Arguments

- mu:

  Control parameter (asymmetry in angular momentum flux). Corresponds to
  experimental γ. Values near 0 produce random point attractors; values
  around 0.02-0.04 produce noisy periodic motion; values \> 0.06 produce
  chaotic attractors.

- n_steps:

  Number of integration steps.

- dt:

  Time step for Euler-Maruyama integration.

- x0, y0, z0:

  Initial conditions for state variables.

- a:

  Damping coefficient (default 0.2 from paper).

- phi:

  Ornstein-Uhlenbeck relaxation rate (default 0.9 from paper).

- sigma:

  Noise amplitude (default 0.2 from paper).

- omega:

  Forcing frequency (default 1.0 from paper).

- ou_mu:

  Ornstein-Uhlenbeck process mean (default NULL uses mu). The OU process
  reverts to this value with rate phi. By default, the noise mean equals
  the control parameter mu, which is physically motivated. Set
  explicitly for independent control of noise bias.

- seed:

  Random seed for reproducibility.

- thin:

  Thinning factor for output (store every thin-th point).

- transient:

  Number of initial steps to discard as transient.

## Value

Object of class "vonkarman_sim" containing time series and parameters.

## Details

**What this function actually integrates.** The `vonkarman_*` family
does NOT solve the Navier-Stokes equations, nor does it perform a
Galerkin truncation of the von-Karman vortex street. The underlying
model is the *stochastic Duffing oscillator* that Faranda, Sato,
Saint-Michel et al. (2017) introduced as a phenomenological reduced
model for the experimental observables of the turbulent von-Karman
swirling-flow cell (counter-rotating impellers). The function name
reflects the experimental setup the model was built to describe, not the
fluid-mechanical phenomenon of vortex shedding.

The model equations are: \$\$dx = y \\ dt\$\$ \$\$dy = (-ay + x - x^3 +
z\sin(\omega t)) \\ dt\$\$ \$\$dz = -\phi(z - \mu\_{ou}) \\ dt + \sigma
\\ dW_t\$\$

where \\W_t\\ is a Wiener process integrated by Euler-Maruyama. The
variable \\x\\ corresponds to the experimental observable \\\theta\\
(symmetry order parameter measuring the asymmetry in impeller rotation
frequencies), and \\\mu\\ corresponds to the experimental control
parameter \\\gamma\\ (asymmetry in applied torques).

The quasistationary states of the system are located at \\x_s = \pm 1\\.
For visualization matching the experimental attractors, the
transformation \\x_m = \text{sign}(\mu)(x - 1)\\ shifts the attractor
appropriately.

## References

Faranda, D., Sato, Y., Saint-Michel, B., Wiertel, C., Padilla, V.,
Dubrulle, B., & Daviaud, F. (2017). Stochastic chaos in a turbulent
swirling flow. *Physical Review Letters*, 119(1), 014502.
[doi:10.1103/PhysRevLett.119.014502](https://doi.org/10.1103/PhysRevLett.119.014502)

## See also

[`print.vonkarman_sim()`](https://robustecologies.github.io/tuRbulence/reference/print.vonkarman_sim.md),
[`summary.vonkarman_sim()`](https://robustecologies.github.io/tuRbulence/reference/summary.vonkarman_sim.md),
[`plot.vonkarman_sim()`](https://robustecologies.github.io/tuRbulence/reference/plot.vonkarman_sim.md),
[`vonkarman_batch()`](https://robustecologies.github.io/tuRbulence/reference/vonkarman_batch.md),
[`vonkarman_peaks()`](https://robustecologies.github.io/tuRbulence/reference/vonkarman_peaks.md),
[`vonkarman_attractor()`](https://robustecologies.github.io/tuRbulence/reference/vonkarman_attractor.md),
[`simulate_system()`](https://robustecologies.github.io/tuRbulence/reference/simulate_system.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate chaotic attractor (μ = 0.3)
sim <- vonkarman_sim(mu = 0.3, n_steps = 100000, seed = 42)

# Inspect results
print(sim)
summary(sim)

# Visualize time series
plot(sim)
plot(sim, type = "timeseries", var = "theta")

# Phase space and attractor
plot(sim, type = "phase")
plot(sim, type = "attractor")

# Power spectral density
plot(sim, type = "psd")
} # }
```
