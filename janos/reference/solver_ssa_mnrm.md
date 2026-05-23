# Modified next-reaction method SSA solver (Anderson)

Creates a solver specification for Anderson's modified next-reaction
method (Anderson, 2007). Uses internal times via unit-rate Poisson
processes and a dependency graph for selective propensity updates,
achieving O(log M) cost per step with improved numerical stability
compared to the standard NRM.

## Usage

``` r
solver_ssa_mnrm(seed = 42, output_dt = 1, max_reactions = 1e+08)
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

A list of class `solver_spec` with `method = "ssa_mnrm"`.

## References

Anderson, D. F. (2007). A modified next reaction method for simulating
chemical systems with time dependent propensities and delays. *Journal
of Chemical Physics*, 127, 214107.
[doi:10.1063/1.2799998](https://doi.org/10.1063/1.2799998)

## See also

[`solver_ssa_direct`](https://robustecologies.github.io/janos/reference/solver_ssa_direct.md),
[`solver_ssa_nrm`](https://robustecologies.github.io/janos/reference/solver_ssa_nrm.md),
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
result <- dyn_sim(sir, t_max = 200, solver = solver_ssa_mnrm(),
                  discard_transient = 0)
plot(result)
} # }
```
