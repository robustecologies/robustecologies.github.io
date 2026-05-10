# Pivot a long-format observation table back to wide

Inverse of
[`vt_long()`](https://robustecologies.github.io/VerteTIME/reference/vt_long.md).
Spreads species into columns with calendar year as the row index.
Covariate columns can optionally be appended to the right of the species
block.

## Usage

``` r
vt_wide(long, covariates = NULL)
```

## Arguments

- long:

  A `data.table` with columns `site_id`, `species_id`, `year`, `value`
  (covariates ignored unless `covariates` is supplied).

- covariates:

  Optional `data.table` to append as wide columns.

## Value

A `data.table` with `year` as the first column followed by one column
per species (and per covariate if supplied).

## References

Dowle, M., & Srinivasan, A. (2024). *data.table: Extension of
`data.frame`*. R package. <https://r-datatable.com>.

## See also

[`vt_long()`](https://robustecologies.github.io/VerteTIME/reference/vt_long.md),
[`vt_read_csv()`](https://robustecologies.github.io/VerteTIME/reference/vt_read_csv.md)

## Examples

``` r
if (FALSE) { # \dontrun{
lo <- vt_long(vt_read_csv(here::here("data-raw","VT_001","VT_001.csv")),
              dataset_id = "VT_001")
wd <- vt_wide(lo, covariates = attr(lo, "covariates"))
} # }
```
