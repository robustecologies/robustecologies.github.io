# Plot a single climate index

Creates a time series plot of a climate oscillation index with optional
smoothing and extreme value highlighting.

## Usage

``` r
plot_index(
  data,
  index = NULL,
  highlight_extremes = FALSE,
  smooth = FALSE,
  fill_anomaly = FALSE,
  title = NULL,
  subtitle = NULL,
  caption = NULL,
  dygraph = FALSE,
  ...
)
```

## Arguments

- data:

  A tibble containing climate index data. Must have columns `date` and
  `value`. Can be obtained from
  [`get_index()`](https://robustecologies.github.io/oscillatoRs/reference/get_index.md)
  or bundled data.

- index:

  Character. If data contains multiple indices, specify which one to
  plot. If NULL and data contains only one index, plots that index.

- highlight_extremes:

  Logical. If TRUE, highlights values beyond +/- 1.5 standard
  deviations.

- smooth:

  Logical or numeric. If TRUE, adds a LOESS smoother with default span.
  If numeric, uses that value as the span parameter.

- fill_anomaly:

  Logical. If TRUE, fills areas above/below zero.

- title:

  Character. Plot title. If NULL, uses index description.

- subtitle:

  Character. Plot subtitle.

- caption:

  Character. Plot caption. If NULL, uses default attribution.

- dygraph:

  Logical. If TRUE, returns an interactive dygraphs widget instead of a
  ggplot2 object. Requires the dygraphs and xts packages.

- ...:

  Additional arguments (currently unused).

## Value

A ggplot2 object, or a dygraphs widget if `dygraph = TRUE`.

## References

Wilke, C. O. (2019). *Fundamentals of data visualization*. O'Reilly
Media. ISBN 978-1492031086

## Examples

``` r
if (FALSE) { # \dontrun{
# Using bundled data
data(climate_monthly)
plot_index(climate_monthly, "NAO")

# With smoothing
plot_index(climate_monthly, "PDO", smooth = 0.1)

# With filled anomalies
plot_index(climate_monthly, "Nino34", fill_anomaly = TRUE)

# Interactive dygraphs version
plot_index(climate_monthly, "NAO", dygraph = TRUE)
} # }
```
