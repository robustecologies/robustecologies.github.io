# Sweep solver parameters over a grid

Runs one solver once per point of a parameter grid and tabulates the
outcome, reproducing in one call the (C, h) stability phase diagrams
that guide the tuning of every method in the package. Grid rows are
independent and run in parallel under OpenMP for compiled objectives.

## Usage

``` r
sym_sweep(
  objective,
  method,
  grid,
  x0 = NULL,
  control = sym_control(),
  n_threads = 1L,
  metric = c("f_best", "iters_to_tol", "success"),
  tol = 1e-06,
  seed = NULL
)
```

## Arguments

- objective:

  A `sym_objective`; compile it with
  [`sym_compile()`](https://robustecologies.github.io/symplectoR/reference/sym_compile.md)
  to enable the parallel path.

- method:

  A single method name accepted by
  [`sym_optim()`](https://robustecologies.github.io/symplectoR/reference/sym_optim.md)
  (the trajectory methods; `"qhd"` is not sweepable through this
  driver).

- grid:

  Named list of numeric vectors, one per swept parameter (for example
  `list(C = 10^seq(-3, 1), h = 10^seq(-3, 0))`); the full factorial
  expansion is evaluated. Valid names are the tunables of the chosen
  method as documented in
  [`sym_control()`](https://robustecologies.github.io/symplectoR/reference/sym_control.md).

- x0:

  Start vector shared by every run; resolved as in
  [`sym_optim()`](https://robustecologies.github.io/symplectoR/reference/sym_optim.md).

- control:

  Base `sym_control`; swept parameters override it per row.

- n_threads:

  OpenMP threads over grid rows (compiled objectives only).

- metric:

  Primary metric for the default plot: final incumbent value, iterations
  to the tolerance, or a success indicator.

- tol:

  Objective threshold defining `"success"` when the benchmark minimum is
  known to be zero; a run succeeds when it converged and its incumbent
  value is at most `tol`.

- seed:

  Reserved for future stochastic methods; the trajectory sweeps are
  deterministic.

## Value

An object of class `sym_sweep`: the expanded grid data frame with
columns `f_best`, `f_final`, `n_iter`, `converged`, `diverged`,
`success`, plus the sweep metadata and the incumbent row.

## Details

The convergent region of the `(C, h)` plane for the Bregman leapfrog
methods is nearly independent of the problem dimension (Duruisseaux &
Leok 2023), so a sweep on a low-dimensional analogue transfers to the
production problem. Sweeps are also the honest way to compare stability
margins across methods: structure-preserving integrators hold larger
stable step-size regions than the Nesterov discretization at equal cost.

## References

Duruisseaux, V., & Leok, M. (2023). Practical perspectives on symplectic
accelerated optimization. *Optimization Methods and Software*, 38(6),
1230-1268.
[doi:10.1080/10556788.2023.2214837](https://doi.org/10.1080/10556788.2023.2214837)

## See also

[`sym_optim()`](https://robustecologies.github.io/symplectoR/reference/sym_optim.md),
[`sym_control()`](https://robustecologies.github.io/symplectoR/reference/sym_control.md),
[`sym_compile()`](https://robustecologies.github.io/symplectoR/reference/sym_compile.md),
[`print.sym_sweep()`](https://robustecologies.github.io/symplectoR/reference/print.sym_sweep.md),
[`summary.sym_sweep()`](https://robustecologies.github.io/symplectoR/reference/summary.sym_sweep.md),
[`plot.sym_sweep()`](https://robustecologies.github.io/symplectoR/reference/plot.sym_sweep.md)

## Examples

``` r
if (FALSE) { # \dontrun{
co <- sym_compile("return 0.5 * arma::dot(x, x);", "return x;")
obj <- sym_objective(co, dim = 20)
sw <- sym_sweep(obj, "slc_expo",
                grid = list(C = 10^seq(-2, 2, length.out = 20),
                            h = 10^seq(-2, 1, length.out = 20)),
                x0 = rep(3, 20), n_threads = 8)
print(sw)
plot(sw, type = "heatmap")
plot(sw, type = "dashboard")
} # }
```
