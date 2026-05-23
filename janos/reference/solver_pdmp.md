# PDMP solver (piecewise deterministic Markov process)

Creates a solver specification for piecewise deterministic Markov
processes. The solver uses RK4 for intra-regime ODE integration and the
Lewis-Shedler thinning algorithm for generating inter-event (regime
switching) times. The thinning bound is updated at each putative event
time using the actual transition rates.

## Usage

``` r
solver_pdmp(
  dt_ode = 0.01,
  output_dt = 1,
  seed = 42L,
  max_transitions = 1000000L
)
```

## Arguments

- dt_ode:

  Time step for RK4 integration within each regime (default: 0.01).

- output_dt:

  Interval between stored output points (default: 1.0).

- seed:

  Random seed for reproducible simulations (default: 42).

- max_transitions:

  Maximum number of regime transitions before stopping (default: 1e6).

## Value

A list of class `solver_spec` with `method = "pdmp"`.

## Details

Between switching events the state evolves deterministically according
to the current regime's ODE. At each step, the solver draws a putative
event time from an exponential distribution with rate equal to the sum
of all outgoing transition rates from the current regime (thinning
bound). The ODE is integrated to that time, rates are recomputed, and
the event is accepted with probability lambda_actual / lambda_bar. If
accepted, the target regime is selected proportional to individual
rates.

## References

Lewis, P. A. W. & Shedler, G. S. (1979). Simulation of nonhomogeneous
Poisson processes by thinning. *Naval Research Logistics Quarterly*,
26(3), 403-413.
[doi:10.1002/nav.3800260304](https://doi.org/10.1002/nav.3800260304)

Davis, M. H. A. (1984). Piecewise-deterministic Markov processes: A
general class of non-diffusion stochastic models. *Journal of the Royal
Statistical Society, Series B*, 46(3), 353-388.

## See also

[`system_spec`](https://robustecologies.github.io/janos/reference/system_spec.md),
[`dyn_sim`](https://robustecologies.github.io/janos/reference/dyn_sim.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Two-state bacterial growth model
model <- system_spec(
    regimes = list(
        active  = list(x ~ r * x * (1 - x / K)),
        dormant = list(x ~ -d * x)
    ),
    transitions = list(
        list(from = "active", to = "dormant", rate = ~ lambda_ad),
        list(from = "dormant", to = "active", rate = ~ lambda_da)
    ),
    state_names = "x",
    parms = list(r = 0.5, K = 100, d = 0.01,
                 lambda_ad = 0.1, lambda_da = 0.05),
    init = c(x = 50),
    initial_regime = "active"
)
result <- dyn_sim(model, t_max = 200, solver = solver_pdmp())
print(result)
summary(result)
plot(result)
} # }
```
