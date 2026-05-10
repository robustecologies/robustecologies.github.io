# Classify column names into species and covariates

Heuristic split based on the species-name regex `vt.species_regex`
(`^[A-Z][a-z]+_[a-z][a-z0-9-]*$` by default). The optional
`taxa_columns` and `covariate_columns` arguments allow per-dataset
overrides through the YAML sidecar.

## Usage

``` r
classify_columns(x, taxa_columns = NULL, covariate_columns = NULL)
```

## Arguments

- x:

  Character vector of column names from a wide CSV.

- taxa_columns:

  Optional character vector of column names to force into the taxa
  group.

- covariate_columns:

  Optional character vector of column names to force into the covariates
  group.

## Value

A list with components `taxa` and `covariates`, each a character vector.
The `year` column is always excluded from both.

## See also

[`vt_long()`](https://robustecologies.github.io/VerteTIME/reference/vt_long.md),
[`vt_ingest_dataset()`](https://robustecologies.github.io/VerteTIME/reference/vt_ingest_dataset.md)

## Examples

``` r
if (FALSE) { # \dontrun{
classify_columns(c("year","Geospiza_fortis","Geospiza_scandens","rain"))
} # }
```
