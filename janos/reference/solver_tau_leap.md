# Adaptive tau-leap solver

Creates a solver specification for the adaptive tau-leaping method of
Cao, Gillespie, and Petzold (2006). At each step, the algorithm selects
the largest time increment tau that keeps the expected relative change
in each propensity below `epsilon`. When propensities are too small or
the selected tau is shorter than a few exact-SSA steps, the algorithm
falls back to `n_critical` exact Gillespie steps for accuracy.

## Usage

``` r
solver_tau_leap(
  seed = 42,
  output_dt = 1,
  max_steps = 1e+07,
  epsilon = 0.03,
  n_critical = 10L,
  postleap_check = TRUE,
  max_rejections = 10L
)
```

## Arguments

- seed:

  Random seed for reproducible simulations (default: 42).

- output_dt:

  Time interval between stored output points (default: 1.0).

- max_steps:

  Maximum number of tau-leap steps (default: 1e7).

- epsilon:

  Relative change bound for propensities (default: 0.03). Smaller values
  increase accuracy but reduce step sizes.

- n_critical:

  Number of exact SSA steps to take when tau is too small (default: 10).

- postleap_check:

  Logical; clamp negative populations to zero after each leap when step
  rejection is exhausted (default: TRUE).

- max_rejections:

  Maximum number of tau-halving rejections per step when a leap produces
  negative populations (default: 10). Setting to 0 disables step
  rejection and reverts to the clamping-only behavior.

## Value

A list of class `solver_spec` with `method = "tau_leap"`.

## Details

The adaptive tau is computed as:

\$\$\tau = \min_i \left\\ \frac{\max(\epsilon x_i / g_i,
1)}{\|\mu_i\|},\\ \frac{\max(\epsilon x_i / g_i, 1)^2}{\sigma_i^2}
\right\\\$\$

where \\\mu_i\\ and \\\sigma_i^2\\ are the expected mean and variance of
the state change per unit time, \\x_i\\ is the current population of
species i, and \\g_i\\ is the highest order of reactant i across all
reactions.

## References

Cao, Y., Gillespie, D. T., & Petzold, L. R. (2006). Efficient step size
selection for the tau-leaping simulation method. *Journal of Chemical
Physics*, 124, 044109.
[doi:10.1063/1.2187339](https://doi.org/10.1063/1.2187339)

## See also

[`solver_tau_leap_midpoint`](https://robustecologies.github.io/janos/reference/solver_tau_leap_midpoint.md),
[`solver_hybrid`](https://robustecologies.github.io/janos/reference/solver_hybrid.md),
[`solver_ssa_direct`](https://robustecologies.github.io/janos/reference/solver_ssa_direct.md),
[`dyn_sim`](https://robustecologies.github.io/janos/reference/dyn_sim.md)

## Examples

``` r
if (FALSE) { # \dontrun{
sir <- model_spec(
    stoichiometry = list(
        infection = c(S = -1L, I = 1L, R = 0L),
        recovery  = c(S = 0L,  I = -1L, R = 1L)
    ),
    propensities = list(
        infection ~ beta * S * I,
        recovery  ~ gamma * I
    ),
    state_names = c("S", "I", "R"),
    parms = list(beta = 0.001, gamma = 0.1),
    init = c(S = 999, I = 1, R = 0)
)
result <- dyn_sim(sir, t_max = 200, solver = solver_tau_leap(epsilon = 0.03),
                  discard_transient = 0)
plot(result)
} # }
```
