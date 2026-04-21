# Compute a bifurcation diagram by parameter sweep

Traces equilibrium branches as a function of a single parameter by
solving for all equilibria at each parameter value using Newton's method
with multiple initial guesses. This is a robust alternative to
pseudo-arclength continuation
([`continuation`](https://robustecologies.github.io/janos/reference/continuation.md))
that handles stiff systems, fold bifurcations, and disconnected branches
without difficulty.

## Usage

``` r
bifurcation_sweep(
  model,
  par_name,
  par_range,
  n_points = 200L,
  init_guesses = NULL,
  state_range = NULL,
  eq_tol = 1e-08,
  merge_tol = NULL,
  max_equilibria = 20L,
  parallel = TRUE,
  n_cores = NULL,
  verbose = TRUE
)
```

## Arguments

- model:

  A
  [`model_spec`](https://robustecologies.github.io/janos/reference/model_spec.md)
  object with formula-based RHS.

- par_name:

  Character string naming the continuation parameter.

- par_range:

  Numeric vector of length 2 giving the parameter range.

- n_points:

  Integer; number of parameter values to sample. Default: 200.

- init_guesses:

  Optional list of named numeric vectors providing initial guesses for
  equilibrium search at each parameter value. These are used in addition
  to warm-start from previous solutions.

- state_range:

  Optional named list with elements for each state variable giving
  `c(min, max)` bounds for the automatic guess grid. If `NULL`,
  estimated from the initial model.

- eq_tol:

  Numeric; convergence tolerance for Newton iteration. Default: 1e-8.

- merge_tol:

  Numeric; distance below which equilibria are merged as duplicates. If
  `NULL`, set to one percent of the state space diagonal.

- max_equilibria:

  Integer; maximum equilibria per parameter value. Default: 20.

- parallel:

  Logical; dispatch the per-parameter Newton solves across workers. On
  Unix this uses
  [`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html); on
  Windows a PSOCK cluster is created. When `TRUE`, warm-start from the
  previous parameter value is disabled (each worker is independent), so
  convergence may occasionally miss branches that are only reachable by
  continuation from a nearby solution. Default: TRUE.

- n_cores:

  Integer or NULL; number of workers. Defaults to
  `parallel::detectCores() - 1`. Forced to 1 inside R CMD check.

- verbose:

  Logical; print progress messages. Default: TRUE.

## Value

An S3 object of class `bifurcation_sweep` with elements `branches` (data
frame with columns `par`, `branch_id`, state variables, `stable`,
`classification`), `bifurcations` (data frame of detected bifurcation
points), `par_name`, `par_range`, `state_names`, `n_points`, and
`model`.

## Details

The algorithm proceeds as follows. For each value of the continuation
parameter, the model is instantiated with the updated parameter and
equilibria are found via Newton iteration from a pool of initial
guesses. The guess pool consists of user-supplied seeds, solutions found
at the previous parameter value (warm-start for branch continuity), and
an automatic grid over the state space. Duplicate solutions are merged
within `merge_tol` Euclidean distance.

After the sweep, equilibria across parameter values are connected into
branches by nearest-neighbour matching in state space. Bifurcation
points are detected where branches appear, disappear, or change
stability.

## See also

[`continuation`](https://robustecologies.github.io/janos/reference/continuation.md),
[`equilibrium`](https://robustecologies.github.io/janos/reference/equilibrium.md),
[`bifurcation_diagram`](https://robustecologies.github.io/janos/reference/bifurcation_diagram.md),
[`model_spec`](https://robustecologies.github.io/janos/reference/model_spec.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Pitchfork bifurcation in double-well system
dw <- model_spec(
    rhs = list(x ~ a * x - x^3),
    state_names = "x",
    parms = list(a = 1),
    init = c(x = 0)
)

bs <- bifurcation_sweep(dw, "a", c(-1, 2), n_points = 100)
print(bs)
summary(bs)
plot(bs)
} # }
```
