# Summary statistics for monthly climate indices

Pre-computed summary statistics for each monthly climate index,
including temporal coverage, completeness, and basic descriptive
statistics.

## Usage

``` r
climate_summary
```

## Format

A tibble with the following columns:

- index:

  Character. Index identifier.

- category:

  Character. Index category.

- description:

  Character. Short description.

- start_year:

  Integer. First year with data.

- end_year:

  Integer. Last year with data.

- n_obs:

  Integer. Number of non-missing observations.

- n_missing:

  Integer. Number of missing values.

- pct_complete:

  Numeric. Percentage of complete observations.

- mean:

  Numeric. Mean value.

- sd:

  Numeric. Standard deviation.

- min:

  Numeric. Minimum value.

- max:

  Numeric. Maximum value.

## Examples

``` r
if (FALSE) { # \dontrun{
data(climate_summary)

# Indices with longest records
climate_summary[order(climate_summary$start_year), c("index", "start_year", "end_year")]
} # }
```
