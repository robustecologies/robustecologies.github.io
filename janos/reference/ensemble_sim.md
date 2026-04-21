# Simulate an ensemble of independent replicates

Runs multiple independent simulations of a dynamical system model and
collects trajectory-level and terminal-state statistics. Two
parallelization backends are available: compiled C++ OpenMP batch
templates provide maximum throughput for SSA Direct, SDE Euler-Maruyama,
and adaptive tau-leap solvers, while an R-level parallel fallback
supports every solver type in the package.

## Usage

``` r
ensemble_sim(
  model,
  n_replicates = 100,
  t_max = 50,
  solver = solver_ssa_direct(),
  init = NULL,
  parms = NULL,
  vary = NULL,
  seed = 42,
  parallel = TRUE,
  n_cores = NULL,
  backend = c("auto", "openmp", "mclapply", "future"),
  store_trajectories = TRUE,
  summary_quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975),
  discard_transient = 0,
  verbose = TRUE
)
```

## Arguments

- model:

  A
  [`model_spec`](https://robustecologies.github.io/janos/reference/model_spec.md)
  object.

- n_replicates:

  Number of independent replicates (default: 100).

- t_max:

  Simulation duration (default: 50).

- solver:

  A `solver_spec` object (default:
  [`solver_ssa_direct()`](https://robustecologies.github.io/janos/reference/solver_ssa_direct.md)).

- init:

  Optional named numeric vector of initial conditions overriding the
  model defaults.

- parms:

  Optional named list of parameters overriding the model defaults.

- vary:

  Optional list with elements `init` and/or `parms`, each a function of
  replicate index. See Details.

- seed:

  Master random seed (default: 42). Replicate i uses seed + i - 1.

- parallel:

  Logical; use parallel execution (default: TRUE).

- n_cores:

  Number of parallel workers. Defaults to `parallel::detectCores() - 1`,
  capped at `n_replicates`.

- backend:

  Backend selection: `"auto"` (default), `"openmp"`, `"mclapply"`, or
  `"future"`.

- store_trajectories:

  Logical; store all trajectory matrices (default: TRUE). Set to FALSE
  for memory-efficient large ensembles.

- summary_quantiles:

  Numeric vector of quantile levels for summary statistics (default:
  `c(0.025, 0.25, 0.5, 0.75, 0.975)`).

- discard_transient:

  Time to discard as transient dynamics (default: 0).

- verbose:

  Logical; print progress messages (default: TRUE).

## Value

An S3 object of class `ensemble_sim` containing:

- trajectories:

  List of trajectory matrices, or NULL if `store_trajectories = FALSE`

- summary:

  List with `time`, `mean`, `sd`, and `quantiles` matrices, or NULL if
  trajectories not stored

- terminal_states:

  Matrix (n_replicates x n_states) of final values

- extinction:

  Named list of logical vectors per state (did state reach zero at any
  point), only for CTMC models

- n_replicates:

  Number of replicates

- t_max:

  Simulation duration

- model:

  The model_spec used

- parms:

  Parameters used

- solver:

  The solver_spec used

- metadata:

  List with backend, n_cores, elapsed_seconds, n_extinctions per state

## Details

When `backend = "auto"` (the default), the function selects the fastest
available backend. The OpenMP batch path is chosen when the solver is
one of
[`solver_ssa_direct()`](https://robustecologies.github.io/janos/reference/solver_ssa_direct.md),
[`solver_euler_maruyama()`](https://robustecologies.github.io/janos/reference/solver_euler_maruyama.md),
or
[`solver_tau_leap()`](https://robustecologies.github.io/janos/reference/solver_tau_leap.md),
and `vary` is NULL (all replicates share identical initial conditions
and parameters). In all other cases the function falls back to R-level
parallelism via
[`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html) or,
when requested,
[`future.apply::future_lapply()`](https://future.apply.futureverse.org/reference/future_lapply.html).

Reproducibility is guaranteed regardless of backend or thread count:
each replicate i uses seed = `seed + i - 1`, so results depend only on
the master seed and replicate index.

For large ensembles where memory is a concern, setting
`store_trajectories = FALSE` avoids allocating the full 3D trajectory
array. In this mode only terminal states and running mean and SD are
kept.

The `vary` argument allows per-replicate initial conditions or
parameters. When non-NULL it must be a list with optional elements
`init` (a function of replicate index returning a named numeric vector)
and `parms` (a function of replicate index returning a named list).
Using `vary` forces the R-level parallel backend because each replicate
may differ structurally.

## See also

[`model_spec`](https://robustecologies.github.io/janos/reference/model_spec.md),
[`dyn_sim`](https://robustecologies.github.io/janos/reference/dyn_sim.md),
[`solver_ssa_direct`](https://robustecologies.github.io/janos/reference/solver_ssa_direct.md),
[`solver_euler_maruyama`](https://robustecologies.github.io/janos/reference/solver_euler_maruyama.md),
[`solver_tau_leap`](https://robustecologies.github.io/janos/reference/solver_tau_leap.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Birth-death SSA ensemble
bd <- model_spec(
    stoichiometry = list(
        birth = c(N = 1L),
        death = c(N = -1L)
    ),
    propensities = list(
        birth ~ lambda * N,
        death ~ mu * N
    ),
    state_names = "N",
    parms = list(lambda = 1.0, mu = 1.1),
    init = c(N = 100)
)

ens <- ensemble_sim(bd, n_replicates = 200, t_max = 20,
                    solver = solver_ssa_direct())
print(ens)
summary(ens)
plot(ens)
plot(ens, type = "terminal")
plot(ens, type = "extinction")
} # }
```
