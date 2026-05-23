# Next-reaction method SSA solver (Gibson-Bruck)

Creates a solver specification for the Gibson-Bruck next-reaction method
(Gibson & Bruck, 2000). Uses a priority queue of putative firing times,
achieving O(log M) cost per step where M is the number of reactions.
After firing a reaction, only affected propensities are recomputed and
their putative times rescaled, avoiding full recomputation.

## Usage

``` r
solver_ssa_nrm(seed = 42, output_dt = 1, max_reactions = 1e+08)
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

A list of class `solver_spec` with `method = "ssa_nrm"`.

## References

Gibson, M. A., & Bruck, J. (2000). Efficient exact stochastic simulation
of chemical systems with many species and many channels. *Journal of
Physical Chemistry A*, 104(9), 1876-1889.
[doi:10.1021/jp993732q](https://doi.org/10.1021/jp993732q)

## See also

[`solver_ssa_direct`](https://robustecologies.github.io/janos/reference/solver_ssa_direct.md),
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
result <- dyn_sim(sir, t_max = 200, solver = solver_ssa_nrm(),
                  discard_transient = 0)
plot(result)
} # }
```
