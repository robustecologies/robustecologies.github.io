# Midpoint tau-leap solver

Creates a solver specification for the midpoint tau-leaping method. This
variant evaluates propensities at the midpoint of the tau interval for
improved accuracy compared to the explicit tau-leap, at the cost of an
additional propensity evaluation per step.

## Usage

``` r
solver_tau_leap_midpoint(
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

A list of class `solver_spec` with `method = "tau_leap_midpoint"`.

## Details

The midpoint method first computes the deterministic drift to obtain a
midpoint state \\x\_{mid} = x + (\tau/2) \cdot \sum_j a_j(x) \nu_j\\,
then draws Poisson-distributed reaction counts using propensities
evaluated at \\x\_{mid}\\. This reduces the bias inherent in explicit
tau-leaping by one order.

## References

Gillespie, D. T. (2001). Approximate accelerated stochastic simulation
of chemically reacting systems. *Journal of Chemical Physics*, 115(4),
1716-1733. [doi:10.1063/1.1378322](https://doi.org/10.1063/1.1378322)

## See also

[`solver_tau_leap`](https://robustecologies.github.io/janos/reference/solver_tau_leap.md),
[`solver_hybrid`](https://robustecologies.github.io/janos/reference/solver_hybrid.md),
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
result <- dyn_sim(sir, t_max = 200,
                  solver = solver_tau_leap_midpoint(epsilon = 0.03),
                  discard_transient = 0)
plot(result)
} # }
```
