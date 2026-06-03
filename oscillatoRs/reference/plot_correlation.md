# Plot correlation matrix of climate indices

Computes pairwise correlations between climate indices and displays them
as a heatmap tile plot.

## Usage

``` r
plot_correlation(
  data,
  indices = NULL,
  method = c("pearson", "spearman", "kendall"),
  start_year = NULL,
  min_obs = 100,
  title = NULL
)
```

## Arguments

- data:

  A tibble in wide format with date column and one column per index. Can
  be obtained from
  [`get_indices()`](https://robustecologies.github.io/oscillatoRs/reference/get_indices.md)
  with format = "wide" or from bundled `climate_monthly_wide`.

- indices:

  Character vector. Subset of indices to include. If NULL, uses all
  numeric columns.

- method:

  Character. Correlation method: "pearson" (default), "spearman", or
  "kendall".

- start_year:

  Integer. Filter to data starting from this year.

- min_obs:

  Integer. Minimum pairwise observations required. Correlations with
  fewer observations are set to NA.

- title:

  Character. Plot title.

## Value

A ggplot2 object.

## References

Wilke, C. O. (2019). *Fundamentals of data visualization*. O'Reilly
Media. ISBN 978-1492031086

## Examples

``` r
if (FALSE) { # \dontrun{
data(climate_monthly_wide)
plot_correlation(climate_monthly_wide, start_year = 1980)

# ENSO indices only
enso_cols <- c("SOI", "Nino12", "Nino3", "Nino34", "Nino4", "ONI")
plot_correlation(climate_monthly_wide, indices = enso_cols)
} # }
```
