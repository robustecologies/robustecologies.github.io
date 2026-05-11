# Classify column names into species and covariates

Heuristic split of a wide community-CSV header into a species block and
a covariates block. The split is governed by the package option
`vt.species_regex`; per-dataset overrides flow through the
`taxa_columns` and `covariate_columns` arguments.

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

## Details

The default regex is `^[A-Z][a-z]+_[a-z][a-z0-9-]*$`, matching the
`Genus_species` underscore convention (initial capital genus, lowercase
epithet with optional hyphens and digits). When `taxa_columns` or
`covariate_columns` is supplied, the explicit assignment wins and any
leftover column is funneled through the regex. The `year` column is
excluded unconditionally because it is the temporal index, not a
covariate. The regex is configurable via
`options(vt.species_regex = ...)` to accommodate alternative naming
conventions in private forks.

## References

Wickham, H. (2014). *Tidy Data*. Journal of Statistical Software,
59(10).
[doi:10.18637/jss.v059.i10](https://doi.org/10.18637/jss.v059.i10) .

## See also

[`vt_long()`](https://robustecologies.github.io/VerteTIME/reference/vt_long.md),
[`vt_read()`](https://robustecologies.github.io/VerteTIME/reference/vt_read.md),
[`vt_read_csv()`](https://robustecologies.github.io/VerteTIME/reference/vt_read_csv.md),
[`vt_ingest_dataset()`](https://robustecologies.github.io/VerteTIME/reference/vt_ingest_dataset.md)

## Examples

``` r
classify_columns(c("year", "Geospiza_fortis", "Geospiza_scandens", "rain"))
#> $taxa
#> [1] "Geospiza_fortis"   "Geospiza_scandens"
#> 
#> $covariates
#> [1] "rain"
#> 
```
