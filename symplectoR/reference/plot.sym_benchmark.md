# Plot a benchmark problem

Landscape views of a benchmark
([`plot.sym_objective()`](https://robustecologies.github.io/symplectoR/reference/plot.sym_objective.md)
applied to the wrapped objective) plus the `"analytic"` view showing the
exact continuous-time trace where one exists.

## Usage

``` r
# S3 method for class 'sym_benchmark'
plot(
  x,
  type = c("contour", "surface", "profiles", "analytic", "dashboard"),
  t_max = 20,
  ...
)
```

## Arguments

- x:

  A `sym_benchmark` object.

- type:

  View: `"contour"`, `"surface"`, `"profiles"` (via the objective plot),
  `"analytic"` (the closed-form trajectory of the oscillator
  benchmarks), or `"dashboard"` (a composite adding the analytic trace
  to the landscape views).

- t_max:

  Time horizon for the analytic view.

- ...:

  Passed on to
  [`plot.sym_objective()`](https://robustecologies.github.io/symplectoR/reference/plot.sym_objective.md).

## Value

The ggplot object, invisibly; for `type = "dashboard"` the assembled
grid graphical object.

## See also

[`sym_benchmark()`](https://robustecologies.github.io/symplectoR/reference/sym_benchmark.md),
[`print.sym_benchmark()`](https://robustecologies.github.io/symplectoR/reference/print.sym_benchmark.md),
[`summary.sym_benchmark()`](https://robustecologies.github.io/symplectoR/reference/summary.sym_benchmark.md)

## Examples

``` r
if (FALSE) { # \dontrun{
plot(sym_benchmark("damped_oscillator"), type = "analytic")
plot(sym_benchmark("damped_oscillator"), type = "dashboard")
} # }
```
