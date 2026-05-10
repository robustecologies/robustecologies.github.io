# Pivot the wide observation block of a VerteTIME dataset into long form

Reshapes a wide-format community CSV (years as rows, species as columns)
into long form (one row per `(site_id, species_id, year)`). Covariate
columns identified by
[`classify_columns()`](https://robustecologies.github.io/VerteTIME/reference/classify_columns.md)
are pushed to a separate long frame returned via the `covariates`
attribute.

## Usage

``` r
vt_long(
  wide,
  dataset_id,
  site_id = dataset_id,
  taxa_columns = NULL,
  covariate_columns = NULL
)
```

## Arguments

- wide:

  A `data.table` produced by
  [`vt_read_csv()`](https://robustecologies.github.io/VerteTIME/reference/vt_read_csv.md)
  for a single dataset (or a single site of a multi-site dataset).

- dataset_id:

  Character scalar dataset identifier; written into the output rows.

- site_id:

  Character scalar site identifier; written into the output rows.

- taxa_columns:

  Optional character vector overriding the regex-based species/covariate
  split. Pass column names to be treated as taxa.

- covariate_columns:

  Optional character vector pairing with `taxa_columns`; remaining
  columns become covariates.

## Value

A `data.table` with columns `dataset_id`, `site_id`, `species_id`,
`year`, `value` (one row per non-NA observation, NAs preserved). The
`covariates` attribute carries a parallel `data.table` keyed by
`(dataset_id, site_id, year, covariate_name)`.

## References

Wickham, H. (2014). *Tidy Data*. Journal of Statistical Software,
59(10).
[doi:10.18637/jss.v059.i10](https://doi.org/10.18637/jss.v059.i10) .

## See also

[`vt_wide()`](https://robustecologies.github.io/VerteTIME/reference/vt_wide.md),
[`vt_read_csv()`](https://robustecologies.github.io/VerteTIME/reference/vt_read_csv.md),
[`vt_ingest_dataset()`](https://robustecologies.github.io/VerteTIME/reference/vt_ingest_dataset.md)

## Examples

``` r
if (FALSE) { # \dontrun{
x <- vt_read_csv(here::here("data-raw","VT_001","VT_001.csv"))
lo <- vt_long(x, dataset_id = "VT_001", site_id = "VT_001")
attr(lo, "covariates")
} # }
```
