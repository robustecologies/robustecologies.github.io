# Plot method for an observed_dyn_sim object

Overlays true and observed time series for each observed variable in a
faceted ggplot layout. Other plot types are delegated to
[`plot.dyn_sim()`](https://robustecologies.github.io/janos/reference/plot.dyn_sim.md).

## Usage

``` r
# S3 method for class 'observed_dyn_sim'
plot(x, type = c("comparison", "timeseries", "attractor", "phase"), ...)
```

## Arguments

- x:

  An `observed_dyn_sim` object.

- type:

  Plot type. `"comparison"` (default) overlays true and observed; other
  values are passed through to
  [`plot.dyn_sim()`](https://robustecologies.github.io/janos/reference/plot.dyn_sim.md).

- ...:

  Additional arguments passed to the delegated plotter.

## Value

A ggplot object.

## See also

[`observe()`](https://robustecologies.github.io/janos/reference/observe.md),
[`print.observed_dyn_sim()`](https://robustecologies.github.io/janos/reference/print.observed_dyn_sim.md),
[`summary.observed_dyn_sim()`](https://robustecologies.github.io/janos/reference/summary.observed_dyn_sim.md).

## Examples

``` r
if (FALSE) { # \dontrun{
run <- dyn_sim(lorenz, t_max = 10, solver = solver_rk45())
obs <- observe(run, obs_model = list(type = "gaussian", sd = 0.1))
plot(obs)
} # }
```
