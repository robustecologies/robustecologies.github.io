# Read a VerteTIME dataset

Polymorphic reader. Accepts a canonical `dataset_id` (e.g. `"VT_001"`),
in which case the function returns a single-dataset slice of the shipped
`vertetime` compilation; alternatively accepts an absolute path to a CSV
or Parquet file, in which case the function delegates to
[`vt_read_csv()`](https://robustecologies.github.io/VerteTIME/reference/vt_read_csv.md)
or to the `arrow` reader.

## Usage

``` r
vt_read(x, data_raw = NULL)
```

## Arguments

- x:

  Either a `dataset_id` (matching the canonical pattern
  `^[A-Z]+_[A-Za-z0-9]+$`) or a path to a CSV or Parquet file.

- data_raw:

  Optional path to a maintainer-side ingestion tree. When non-`NULL` and
  the directory exists, the function prefers a staged Parquet under
  `data_raw/_staging/<id>.parquet` and falls back to the wide CSV under
  `data_raw/<id>/<id>.csv`. When `NULL` (the default), the function
  reconstructs the dataset from the shipped
  [vertetime](https://robustecologies.github.io/VerteTIME/reference/vertetime.md)
  compilation, which is the user-facing code path. End users never need
  to supply `data_raw`.

## Value

A
[`vt_dataset()`](https://robustecologies.github.io/VerteTIME/reference/vt_dataset.md)
when `x` resolves to a `dataset_id`; a `data.table` when `x` resolves to
a raw file path.

## Details

On the installed package, the `dataset_id` branch loads the shipped
`vertetime` compilation into a private environment, filters every
relational table by `dataset_id == x`, and rebuilds a
[`vt_dataset()`](https://robustecologies.github.io/VerteTIME/reference/vt_dataset.md)
from the slice. The user's session-level objects are not touched. When
the identifier is absent from the compilation, the function aborts and
prints the first twenty available identifiers as a hint. The
parquet-and-CSV maintainer path is kept for forks that maintain a
private ingestion tree, gated by `data_raw` pointing to an existing
directory.

## References

Apache Software Foundation. (2024). *Apache Parquet*.
<https://parquet.apache.org>.

## See also

[`vt_ingest_dataset()`](https://robustecologies.github.io/VerteTIME/reference/vt_ingest_dataset.md),
[`vt_ingest_all()`](https://robustecologies.github.io/VerteTIME/reference/vt_ingest_all.md),
[`vt_long()`](https://robustecologies.github.io/VerteTIME/reference/vt_long.md),
[`vt_wide()`](https://robustecologies.github.io/VerteTIME/reference/vt_wide.md),
[vertetime](https://robustecologies.github.io/VerteTIME/reference/vertetime.md)

## Examples

``` r
if (FALSE) { # \dontrun{
d <- vt_read("VT_001")
print(d); summary(d); plot(d, type = "whittaker")
} # }
```
