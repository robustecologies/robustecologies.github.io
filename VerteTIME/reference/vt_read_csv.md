# Read a wide-format community CSV from `data-raw/<id>/<id>.csv`

Loads a single dataset CSV into a `data.table` while transparently
absorbing the format quirks observed across the 145 source files:
variable header quoting, leading UTF-8 byte-order marks, the literal `;`
separator, and the `Genus_species` underscore convention for taxon
columns. The reader does not classify columns into species versus
covariates; that step happens during ingestion (see
[`vt_ingest_dataset()`](https://robustecologies.github.io/VerteTIME/reference/vt_ingest_dataset.md)).

## Usage

``` r
vt_read_csv(path, sep = ";", na_strings = c("", "NA", "NULL"))
```

## Arguments

- path:

  Character path to the CSV file. Absolute paths or paths resolvable
  through [`here::here()`](https://here.r-lib.org/reference/here.html)
  are both accepted.

- sep:

  Single-character field separator. Defaults to `";"`, the convention of
  the compilation.

- na_strings:

  Character vector of NA tokens. Defaults to `c("", "NA", "NULL")` to
  absorb both blank cells and the textual `NULL` sentinel used by some
  upstream sources.

## Value

A `data.table` with one row per calendar year, the first column coerced
to integer `year`, headers de-quoted and BOM-stripped, and trailing
all-empty columns dropped. Year values are validated against the
`vt.year_min` and `vt.year_max` options.

## Details

The reader is the single entry point for all VerteTIME ingestion paths.
It is exported so that users can inspect a CSV without going through the
full pipeline. Round-trip parity with
[`vt_export()`](https://robustecologies.github.io/VerteTIME/reference/vt_export.md)
CSV outputs is guaranteed because the writer uses the same separator and
NA conventions.

## References

Frictionless Data. (2024). *Data Resource specification*. Open Knowledge
Foundation. <https://specs.frictionlessdata.io/data-resource/>.

## See also

[`vt_ingest_dataset()`](https://robustecologies.github.io/VerteTIME/reference/vt_ingest_dataset.md),
[`vt_long()`](https://robustecologies.github.io/VerteTIME/reference/vt_long.md),
[`vt_wide()`](https://robustecologies.github.io/VerteTIME/reference/vt_wide.md),
[`vt_export()`](https://robustecologies.github.io/VerteTIME/reference/vt_export.md)

## Examples

``` r
if (FALSE) { # \dontrun{
p <- here::here("data-raw", "VT_001", "VT_001.csv")
x <- vt_read_csv(p)
head(x)
} # }
```
