# Plot a parameter sweep

Views of the sweep outcome over the parameter grid.

## Usage

``` r
# S3 method for class 'sym_sweep'
plot(x, type = c("heatmap", "slice", "pareto", "marginal", "dashboard"), ...)
```

## Arguments

- x:

  A `sym_sweep` object.

- type:

  View: `"heatmap"` (metric over a two-parameter grid; the stability
  phase diagram), `"slice"` (metric against a single swept parameter),
  `"pareto"` (incumbent value against iterations, colored by
  convergence), `"marginal"` (success fraction against each swept
  parameter value), or `"dashboard"` (a composite of the views available
  for this grid).

- ...:

  Ignored.

## Value

The ggplot object, invisibly; for `type = "dashboard"` the assembled
grid graphical object.

## See also

[`sym_sweep()`](https://robustecologies.github.io/symplectoR/reference/sym_sweep.md),
[`print.sym_sweep()`](https://robustecologies.github.io/symplectoR/reference/print.sym_sweep.md),
[`summary.sym_sweep()`](https://robustecologies.github.io/symplectoR/reference/summary.sym_sweep.md)

## Examples

``` r
if (FALSE) { # \dontrun{
sw <- sym_sweep(sym_objective(sym_benchmark("quadratic", d = 4)), "slc_expo",
                grid = list(C = 10^seq(-2, 2, 0.25), h = 10^seq(-2, 1, 0.25)),
                x0 = rep(2, 4))
plot(sw, type = "heatmap")
plot(sw, type = "dashboard")
} # }
```
