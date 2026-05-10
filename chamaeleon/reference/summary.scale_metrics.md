# Summary method for scale_metrics objects

Compute detailed diagnostics for scale-dependent dynamical metrics,
including trend tests and chameleon behavior indicators. This summary
provides a rigorous assessment of whether the attractor exhibits
scale-dependent properties.

## Usage

``` r
# S3 method for class 'scale_metrics'
summary(object, ...)
```

## Arguments

- object:

  Object of class `scale_metrics`.

- ...:

  Additional arguments (ignored).

## Value

An object of class `summary.scale_metrics` containing:

- n_scales:

  Number of frequency scales analyzed.

- freq_range:

  Range of frequencies (min, max).

- D_trend:

  List with Spearman correlation, p-value, and direction.

- theta_trend:

  List with Spearman correlation, p-value, and direction.

- D_range:

  Range of mean dimension across scales.

- D_cv:

  Coefficient of variation of mean dimension.

- theta_range:

  Range of mean persistence across scales.

- theta_cv:

  Coefficient of variation of mean persistence.

- chameleon_assessment:

  Character string with overall assessment.

## Details

The summary performs Spearman rank correlation tests to assess whether
D(f) and theta(f) show significant trends with frequency. A chameleon
attractor is indicated by significant (p \< 0.05) trends in either
metric.

**Trend interpretation:**

- Positive D trend: dimension increases with frequency (more complex at
  small scales).

- Negative D trend: dimension decreases with frequency (more complex at
  large scales).

- Positive theta trend: dynamics become more stochastic at higher
  frequencies.

- Negative theta trend: dynamics become more persistent at higher
  frequencies.

## See also

[`scale_dependent_metrics`](https://robustecologies.github.io/chamaeleon/reference/scale_dependent_metrics.md)
for computing the metrics;
[`print.scale_metrics`](https://robustecologies.github.io/chamaeleon/reference/print.scale_metrics.md),
[`plot.scale_metrics`](https://robustecologies.github.io/chamaeleon/reference/plot.scale_metrics.md);
[`chameleon_test`](https://robustecologies.github.io/chamaeleon/reference/chameleon_test.md)
for rigorous statistical testing.
