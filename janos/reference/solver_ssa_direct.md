# Gillespie direct method SSA solver

Creates a solver specification for the Gillespie direct method
(Gillespie, 1977). At each step, all propensities are computed, the time
to the next event is drawn from an exponential distribution with rate
equal to the sum of propensities, and the reaction is selected
proportional to its propensity. Computational cost is O(M) per step
where M is the number of reactions.

## Usage

``` r
solver_ssa_direct(seed = 42, output_dt = 1, max_reactions = 1e+08)
```

## Arguments

- seed:

  Random seed for reproducible simulations (default: 42).

- output_dt:

  Time interval between stored output points (default: 1.0).

- max_reactions:

  Maximum number of reaction events before stopping (default: 1e8).
  Prevents infinite loops in systems that do not reach an absorbing
  state.

## Value

A list of class `solver_spec` with `method = "ssa_direct"`.

## References

Gillespie, D. T. (1977). Exact stochastic simulation of coupled chemical
reactions. *Journal of Physical Chemistry*, 81(25), 2340-2361.
[doi:10.1021/j100540a008](https://doi.org/10.1021/j100540a008)

## See also

[`solver_ssa_nrm`](https://robustecologies.github.io/janos/reference/solver_ssa_nrm.md),
[`solver_ssa_mnrm`](https://robustecologies.github.io/janos/reference/solver_ssa_mnrm.md),
[`dyn_sim`](https://robustecologies.github.io/janos/reference/dyn_sim.md)

## Examples

``` r
if (FALSE) { # \dontrun{
sir <- system_spec(
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
result <- dyn_sim(sir, t_max = 200, solver = solver_ssa_direct(),
                  discard_transient = 0)
plot(result)
} # }
```
