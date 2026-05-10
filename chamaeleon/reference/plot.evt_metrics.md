# Plot EVT metrics

Visualize the instantaneous dimension d and inverse persistence theta
computed by extreme value theory analysis. This S3 method works with
objects of class `evt_metrics` returned by
[`evt_metrics`](https://robustecologies.github.io/chamaeleon/reference/evt_metrics.md).

## Usage

``` r
# S3 method for class 'evt_metrics'
plot(x, type = c("timeseries", "histogram", "scatter"), ...)
```

## Arguments

- x:

  Object of class `evt_metrics`.

- type:

  Character. Type of plot: "timeseries" (default) shows both metrics
  over time, "histogram" shows distributions, "scatter" shows d vs theta
  scatter plot.

- ...:

  Additional arguments passed to the plotting functions.

## Value

A patchwork ggplot2 composition, returned invisibly.

## Details

The plot shows two panels: the upper panel displays the time series of
instantaneous dimension d(t), while the lower panel shows the inverse
persistence theta(t). Horizontal lines indicate the mean values.

The instantaneous dimension d measures the local complexity of the
attractor, the number of active degrees of freedom around each state.
The inverse persistence theta measures trajectory stability: theta near
1 indicates stochastic behavior, while theta \< 1 indicates persistent
dynamics.

## See also

[`evt_metrics`](https://robustecologies.github.io/chamaeleon/reference/evt_metrics.md)
for computing these metrics,
[`plot.takens_embedding`](https://robustecologies.github.io/chamaeleon/reference/plot.takens_embedding.md)
for attractor visualization.

## Examples

``` r
if (FALSE) { # \dontrun{
embedded <- takens_embed(x, m = 3, tau = 20)
metrics <- evt_metrics(embedded)
plot(metrics)
plot(metrics, type = "histogram")
plot(metrics, type = "scatter")
} # }
```
