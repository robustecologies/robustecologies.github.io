# Plot multiple climate indices

Creates a faceted plot showing multiple climate indices, useful for
comparing temporal patterns across different oscillations.

## Usage

``` r
plot_indices(
  data,
  indices = NULL,
  facet = TRUE,
  ncol = 2,
  scales = "free_y",
  start_year = NULL,
  color_by = c("value", "category", "index"),
  title = NULL,
  dygraph = FALSE,
  ...
)
```

## Arguments

- data:

  A tibble containing climate index data in long format. Must have
  columns `date`, `value`, and `index`.

- indices:

  Character vector. Index names to include. If NULL, includes all
  indices in data.

- facet:

  Logical. If TRUE (default), creates faceted plot. If FALSE, overlays
  all indices on single panel.

- ncol:

  Integer. Number of columns for faceted plot.

- scales:

  Character. Scales argument for facet_wrap: "free_y" (default),
  "fixed", or "free".

- start_year:

  Integer. Filter to data starting from this year.

- color_by:

  Character. Variable to map to color: "value" (default), "category", or
  "index".

- title:

  Character. Plot title.

- dygraph:

  Logical. If TRUE, returns an interactive dygraphs widget displaying
  all indices overlaid on a single plot with synchronized range
  selector. Requires the dygraphs and xts packages.

- ...:

  Additional arguments passed to
  [`plot_index()`](https://robustecologies.github.io/oscillatoRs/reference/plot_index.md)
  for single-index plots.

## Value

A ggplot2 object, or a dygraphs widget if `dygraph = TRUE`.

## References

Wilke, C. O. (2019). *Fundamentals of data visualization*. O'Reilly
Media. ISBN 978-1492031086

## Examples

``` r
if (FALSE) { # \dontrun{
data(climate_monthly)

# ENSO family comparison
plot_indices(
  climate_monthly,
  indices = c("SOI", "Nino34", "ONI", "MEI_v2"),
  start_year = 1980
)

# All Atlantic indices
atlantic <- climate_monthly[climate_monthly$category == "Atlantic", ]
plot_indices(atlantic, color_by = "category")

# Interactive dygraphs version
plot_indices(
  climate_monthly,
  indices = c("SOI", "Nino34", "ONI"),
  dygraph = TRUE
)
} # }
```
