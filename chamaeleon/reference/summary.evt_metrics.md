# Summary method for evt_metrics objects

Compute detailed diagnostics for EVT-based dynamical metrics, including
distribution quantiles, NA proportions, and correlation between
dimension and persistence. These diagnostics help interpret the
attractor properties and identify potential estimation issues.

## Usage

``` r
# S3 method for class 'evt_metrics'
summary(object, ...)
```

## Arguments

- object:

  Object of class `evt_metrics`.

- ...:

  Additional arguments (ignored).

## Value

An object of class `summary.evt_metrics` containing:

- n:

  Number of trajectory points.

- p:

  Embedding dimension.

- q:

  Quantile used for threshold selection.

- d_quantiles:

  Named vector of quantiles (5%, 25%, 50%, 75%, 95%) for d.

- theta_quantiles:

  Named vector of quantiles for theta.

- d_na_prop:

  Proportion of NA values in d.

- theta_na_prop:

  Proportion of NA values in theta.

- d_theta_cor:

  Correlation between d and theta (excluding NAs).

- interpretation:

  List with diagnostic flags and guidance.

## Details

The summary provides interpretation guidance based on the computed
metrics:

**NA proportions:** High NA rates (\>10%) indicate estimation failures,
often due to insufficient data, inappropriate quantile choice, or highly
irregular dynamics. Consider increasing series length or adjusting q.

**Dimension interpretation:** If median d is close to the embedding
dimension p, the dynamics may be noise-dominated or the embedding
dimension is too low. Values much smaller than p suggest low-dimensional
structure.

**Persistence interpretation:** theta near 1 indicates stochastic,
memoryless dynamics. theta significantly below 1 indicates persistent
behavior where trajectories tend to remain near their current state.

**d-theta correlation:** Positive correlation suggests regions of high
dimension also show less persistence (more stochastic). Negative
correlation indicates the opposite pattern.

## See also

[`evt_metrics`](https://robustecologies.github.io/chamaeleon/reference/evt_metrics.md)
for computing the metrics;
[`print.evt_metrics`](https://robustecologies.github.io/chamaeleon/reference/print.evt_metrics.md),
[`plot.evt_metrics`](https://robustecologies.github.io/chamaeleon/reference/plot.evt_metrics.md).
