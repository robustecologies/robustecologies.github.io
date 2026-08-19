# Specify an objective function for symplectic optimization

Constructs the objective object consumed by every solver in symplectoR.
The objective wraps either a plain R function or a compiled function
pointer produced by
[`sym_compile()`](https://robustecologies.github.io/symplectoR/reference/sym_compile.md),
together with the problem dimension, optional box constraints, and the
gradient. When no gradient is supplied, a central finite-difference
gradient is generated automatically.

## Usage

``` r
sym_objective(
  f,
  grad = NULL,
  dim,
  lower = NULL,
  upper = NULL,
  params = NULL,
  vectorized = FALSE,
  constraint = c("none", "barrier", "reflect"),
  barrier_weight = 1e-06,
  fd_step = NULL,
  name = NULL
)
```

## Arguments

- f:

  Objective function. Either an R function taking a numeric vector of
  length `dim` and returning a scalar, a `sym_compiled_objective`
  produced by
  [`sym_compile()`](https://robustecologies.github.io/symplectoR/reference/sym_compile.md),
  or a `sym_benchmark` object produced by
  [`sym_benchmark()`](https://robustecologies.github.io/symplectoR/reference/sym_benchmark.md)
  (in which case the benchmark's gradient, dimension and box are
  inherited unless overridden).

- grad:

  Optional gradient function taking and returning a numeric vector of
  length `dim`. When `NULL`, central finite differences are used with
  step `fd_step`.

- dim:

  Problem dimension (positive integer). Inferred from `f` when `f` is a
  benchmark or compiled objective.

- lower, upper:

  Optional box constraints; scalars are recycled to length `dim`. `NULL`
  means unbounded.

- params:

  Optional numeric vector of fixed parameters forwarded to compiled
  objectives (ignored for R closures, which should capture their
  parameters lexically).

- vectorized:

  Logical; set to `TRUE` when `f` accepts a matrix with one point per
  row and returns a vector of values. Used to accelerate grid
  evaluations in the quantum Hamiltonian descent solver.

- constraint:

  Constraint handling mode. `"none"` evaluates `f` as given; `"barrier"`
  adds a logarithmic barrier `-log(x - lower) - log(upper - x)` scaled
  by `barrier_weight` to the objective and its gradient; `"reflect"`
  folds iterates back into the box by reflection after every position
  update.

- barrier_weight:

  Positive weight of the logarithmic barrier. The barrier displaces the
  constrained minimizer by an amount vanishing with the weight, so
  interior solutions should be checked against the unconstrained
  problem.

- fd_step:

  Finite-difference step. `NULL` selects the automatic step
  `eps^(1/3) * (1 + abs(x))` per coordinate, with `eps` the machine
  epsilon.

- name:

  Optional display name used by print and plot methods.

## Value

An object of class `sym_objective`: a list carrying the callable,
gradient, dimension, box, constraint specification and provenance, ready
to be passed to
[`sym_optim()`](https://robustecologies.github.io/symplectoR/reference/sym_optim.md).

## Details

The objective is the single currency of the package: solvers never see
raw functions. Two evaluation routes exist. The closure route calls back
into R once per evaluation and is never parallelized, because the R
interpreter is not thread safe. The compiled route
([`sym_compile()`](https://robustecologies.github.io/symplectoR/reference/sym_compile.md))
evaluates a C++ function pointer with no R API access, which enables
OpenMP parallelism over multi-start ensembles, parameter sweeps and
pseudospectral grids.

The automatic finite-difference gradient uses central differences,
\$\$g_i \approx \frac{f(x + h_i e_i) - f(x - h_i e_i)}{2 h_i},\$\$ whose
truncation error is second order in the step. Each gradient costs
`2 * dim` objective evaluations; supply an analytic gradient whenever
one is available.

## References

Nocedal, J., & Wright, S. J. (2006). *Numerical optimization* (2nd ed.).
Springer.
[doi:10.1007/978-0-387-40065-5](https://doi.org/10.1007/978-0-387-40065-5)

## See also

[`sym_compile()`](https://robustecologies.github.io/symplectoR/reference/sym_compile.md),
[`sym_benchmark()`](https://robustecologies.github.io/symplectoR/reference/sym_benchmark.md),
[`sym_inverse()`](https://robustecologies.github.io/symplectoR/reference/sym_inverse.md),
[`sym_optim()`](https://robustecologies.github.io/symplectoR/reference/sym_optim.md),
[`print.sym_objective()`](https://robustecologies.github.io/symplectoR/reference/print.sym_objective.md),
[`summary.sym_objective()`](https://robustecologies.github.io/symplectoR/reference/summary.sym_objective.md),
[`plot.sym_objective()`](https://robustecologies.github.io/symplectoR/reference/plot.sym_objective.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## Rosenbrock objective with analytic gradient
obj <- sym_objective(
  f = function(x) 100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2,
  grad = function(x) c(
    -400 * x[1] * (x[2] - x[1]^2) - 2 * (1 - x[1]),
    200 * (x[2] - x[1]^2)
  ),
  dim = 2, name = "rosenbrock"
)
print(obj)
fit <- sym_optim(obj, x0 = c(-1.2, 1), method = "rgd")
summary(fit)
plot(fit, type = "trace")
plot(fit, type = "dashboard")
} # }
```
