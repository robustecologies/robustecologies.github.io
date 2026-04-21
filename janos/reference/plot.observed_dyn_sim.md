# Plot method for observed_dyn_sim objects

Overlays true and observed time series for the observed variables.

## Usage

``` r
# S3 method for class 'observed_dyn_sim'
plot(x, type = c("comparison", "timeseries", "attractor", "phase"), ...)
```

## Arguments

- x:

  An `observed_dyn_sim` object

- type:

  Plot type: `"comparison"` (default) overlays true and observed; other
  types are passed to `plot.dyn_sim`.

- ...:

  Additional arguments passed to the plot function

## Value

A ggplot object.
