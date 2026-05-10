# Read a VerteTIME dataset

Polymorphic reader. Accepts a `dataset_id` (in which case the function
looks for a staged Parquet file under `data-raw/_staging/<id>.parquet`
and falls back to the wide CSV under `data-raw/<id>/<id>.csv`), or an
absolute path to a CSV or Parquet file.

## Usage

``` r
vt_read(x, data_raw = here::here("data-raw"))
```

## Arguments

- x:

  Either a `dataset_id` like `"VT_001"`, or a path to a CSV or Parquet
  file.

- data_raw:

  Root of the source `data-raw/` directory. Defaults to
  `here::here("data-raw")`.

## Value

A `vt_dataset` if `x` is a `dataset_id`, or a `data.table` if `x` is a
raw file path.

## Details

This function is the user-friendly entry point; for batch ingestion of
the entire compilation call
[`vt_ingest_all()`](https://robustecologies.github.io/VerteTIME/reference/vt_ingest_all.md)
directly. The reader prefers the staged Parquet over the wide CSV when
both are available, since Parquet carries the canonical long-format
schema and the species/covariate classification already resolved.

## References

Apache Software Foundation. (2024). *Apache Parquet*.
<https://parquet.apache.org>.

## See also

[`vt_ingest_dataset()`](https://robustecologies.github.io/VerteTIME/reference/vt_ingest_dataset.md),
[`vt_ingest_all()`](https://robustecologies.github.io/VerteTIME/reference/vt_ingest_all.md),
[`vt_long()`](https://robustecologies.github.io/VerteTIME/reference/vt_long.md),
[`vt_wide()`](https://robustecologies.github.io/VerteTIME/reference/vt_wide.md)

## Examples

``` r
if (FALSE) { # \dontrun{
d <- vt_read("VT_001")
x <- vt_read(here::here("data-raw","VT_001","VT_001.csv"))
} # }
```
