# Hybrid SSA/CLE solver

Creates a solver specification for a hybrid stochastic simulation that
partitions reactions into fast and slow channels. Fast reactions are
approximated using the chemical Langevin equation (CLE, Euler-Maruyama
discretization), while slow reactions are simulated exactly via the
Gillespie direct method. The partition is updated every
`repartition_interval` steps.

## Usage

``` r
solver_hybrid(
  seed = 42,
  output_dt = 1,
  max_steps = 1e+07,
  dt_cle = 0.01,
  threshold = 10,
  repartition_interval = 100L
)
```

## Arguments

- seed:

  Random seed for reproducible simulations (default: 42).

- output_dt:

  Time interval between stored output points (default: 1.0).

- max_steps:

  Maximum number of hybrid steps (default: 1e7).

- dt_cle:

  CLE timestep for fast reactions (default: 0.01).

- threshold:

  Propensity times `dt_cle` must exceed this value for a reaction to be
  classified as fast (default: 10.0).

- repartition_interval:

  Recompute the fast/slow partition every this many steps (default:
  100).

## Value

A list of class `solver_spec` with `method = "hybrid"`.

## Details

At each time step of duration `dt_cle`, the algorithm first draws the
time to the next slow reaction. If that time exceeds `dt_cle`, the
entire step uses CLE. Otherwise, CLE runs until the slow event fires,
one slow reaction is selected and applied, and CLE continues for the
remainder. The CLE update for each fast reaction j is:

\$\$x \leftarrow x + \nu_j \left\[ a_j(x) \Delta t + \sqrt{a_j(x) \Delta
t} \\ Z_j \right\]\$\$

where \\Z_j \sim N(0,1)\\. Because the CLE produces continuous-valued
state, the output trajectory contains double values (not integers).

## References

Haseltine, E. L., & Rawlings, J. B. (2002). Approximate simulation of
coupled fast and slow reactions for stochastic chemical kinetics.
*Journal of Chemical Physics*, 117(15), 6959-6969.
[doi:10.1063/1.1505860](https://doi.org/10.1063/1.1505860)

## See also

[`solver_tau_leap`](https://robustecologies.github.io/janos/reference/solver_tau_leap.md),
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
result <- dyn_sim(sir, t_max = 200,
                  solver = solver_hybrid(dt_cle = 0.01),
                  discard_transient = 0)
plot(result)
} # }
```
