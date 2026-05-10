# Plot scale-dependent metrics

Visualize how instantaneous dimension D(t,f) and inverse persistence
theta(t,f) vary across frequency scales. This S3 method works with
objects of class `scale_metrics` returned by
[`scale_dependent_metrics`](https://robustecologies.github.io/chamaeleon/reference/scale_dependent_metrics.md).

## Usage

``` r
# S3 method for class 'scale_metrics'
plot(
  x,
  show_psd = TRUE,
  show_ci = TRUE,
  ci_level = 0.95,
  reference_lines = TRUE,
  log_freq = TRUE,
  interactive = FALSE,
  ...
)
```

## Arguments

- x:

  Object of class `scale_metrics` or `chameleon_analysis`.

- show_psd:

  Logical. Whether to overlay power spectral density (not yet
  implemented).

- show_ci:

  Logical. Whether to show confidence intervals around the mean values.
  Default is TRUE.

- ci_level:

  Numeric. Confidence level for intervals, between 0 and 1. Default is
  0.95 (95 percent CI).

- reference_lines:

  Logical. Add reference lines at d=3 and theta=1 to help interpret the
  results. Default is TRUE.

- log_freq:

  Logical. Use logarithmic scale for frequency axis. Default is TRUE,
  which is appropriate for most natural signals.

- interactive:

  Logical. Use plotly for interactive plots. Default is FALSE (ggplot2 +
  patchwork).

- ...:

  Additional arguments passed to the plotting functions.

## Value

If `interactive=TRUE`, returns a plotly htmlwidget. Otherwise returns a
patchwork ggplot2 composition invisibly.

## Details

The plot shows two panels: the upper panel displays the mean
instantaneous dimension across frequency scales, while the lower panel
shows the mean inverse persistence. Confidence intervals can be
displayed around the means, and reference lines indicate typical values
(d=3 for 3D embedding, theta=1 for stochastic behavior).

The scale-dependent metrics reveal whether the attractor exhibits
chameleon behavior. For a scale-invariant attractor, both D(f) and
theta(f) remain relatively constant across frequencies. For a chameleon
attractor, these metrics vary significantly with frequency scale.

Interpretation guidelines:

- **Dimension D(f)**: Higher values indicate more complex dynamics at
  that frequency scale. D approaching the embedding dimension suggests
  noise-dominated behavior.

- **Persistence theta(f)**: theta near 1 indicates stochastic
  (unpredictable) behavior. theta \< 1 indicates persistent
  (predictable) dynamics.

## See also

[`scale_dependent_metrics`](https://robustecologies.github.io/chamaeleon/reference/scale_dependent_metrics.md)
for computing these metrics,
[`plot.chameleon_analysis`](https://robustecologies.github.io/chamaeleon/reference/plot.chameleon_analysis.md)
for complete analysis visualization.

## Examples

``` r
if (FALSE) { # \dontrun{
result <- chameleon_analysis(x, verbose = TRUE)
plot(result$scale_metrics)
plot(result$scale_metrics, interactive = TRUE)
} # }
```
