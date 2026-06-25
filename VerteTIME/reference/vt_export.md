# Export a compilation to a single file format

Writes a
[`vt_compilation()`](https://robustecologies.github.io/VerteTIME/reference/vt_compilation.md)
to disk in one of four formats: CSV (one file per table plus a
`master.csv` long-form join), Apache Parquet (one file per table plus
`master.parquet`), SQLite single-file with foreign keys, or Frictionless
Data Package (a `datapackage.json` manifest with embedded table schemas
pointing at the CSVs in the same root).

## Usage

``` r
vt_export(x, format = c("csv", "parquet", "sqlite", "datapackage"), path)
```

## Arguments

- x:

  A `vt_compilation`.

- format:

  Character, one of `"csv"`, `"parquet"`, `"sqlite"`, `"datapackage"`.

- path:

  Character output path. For `"csv"` and `"parquet"` this is a directory
  (created if missing); for `"sqlite"` it is a single `.db` filename;
  for `"datapackage"` it is the directory whose root will hold
  `datapackage.json` (sibling to the CSV directory).

## Value

Invisibly, the value of `path`.

## Details

All exports preserve the column order and types of the canonical schema.
A `master.csv` (or `master.parquet`) is a denormalised long-form join of
`observations` with `datasets`, `sites` and `species` aimed at
downstream tools that prefer a single table. SQLite output materialises
the documented relational schema (`inst/schema/sqlite_schema.sql`):
primary keys, enforced foreign keys, surrogate auto-increment ids on
`observations` and `covariates`, and a `value_is_imputed` flag
defaulting to 0. Datapackage output writes a single `datapackage.json`
and points it at the CSV files alongside; if the CSV export has not been
run, this function calls it first.

## References

Frictionless Data. (2024). *Data Package specification*. Open Knowledge
Foundation. <https://specs.frictionlessdata.io/data-package/>.

## See also

[`vt_publish()`](https://robustecologies.github.io/VerteTIME/reference/vt_publish.md),
[`vt_compilation()`](https://robustecologies.github.io/VerteTIME/reference/vt_compilation.md),
[`vt_filter()`](https://robustecologies.github.io/VerteTIME/reference/vt_filter.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(vertetime); co <- vertetime
vt_export(co, "csv",         here::here("web-export","vertetime-v1.0","csv"))
vt_export(co, "parquet",     here::here("web-export","vertetime-v1.0","parquet"))
vt_export(co, "sqlite",      here::here("web-export","vertetime-v1.0","sqlite","vertetime.db"))
vt_export(co, "datapackage", here::here("web-export","vertetime-v1.0"))
} # }
```
