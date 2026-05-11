# Publish a compilation as a complete `web-export/` tree

Orchestrates
[`vt_export()`](https://robustecologies.github.io/VerteTIME/reference/vt_export.md)
for all four formats under `<root>/csv`, `<root>/parquet`,
`<root>/sqlite/vertetime.db`, and the `<root>/datapackage.json`
manifest, then writes `README.md`, `LICENSE.md`, `CITATION.cff`,
`data_provenance.csv` and `CHECKSUMS.sha256` at the root. The output is
the canonical Zenodo deposit payload.

## Usage

``` r
vt_publish(
  x,
  root = here::here("web-export", "vertetime-v1.0"),
  overwrite = FALSE
)
```

## Arguments

- x:

  A `vt_compilation`.

- root:

  Character path to the publish root (will be created).

- overwrite:

  Logical, default `FALSE`. When `FALSE`, existing non-empty roots
  abort.

## Value

Invisibly, the path to `root`.

## Details

The function is deterministic: re-running on the same compilation
reproduces the same SHA-256 checksums (subject to file-mtime variation,
which is excluded from the hash). The `data_provenance.csv` file is
written at the root and is the audit reference for the licence and
per-dataset citation policy.

## References

Hudson, L. N., et al. (2017). *The database of the PREDICTS (Projecting
Responses of Ecological Diversity In Changing Terrestrial Systems)
project*. Ecology and Evolution, 7(1), 145-188.
[doi:10.1002/ece3.2579](https://doi.org/10.1002/ece3.2579) .

## See also

[`vt_export()`](https://robustecologies.github.io/VerteTIME/reference/vt_export.md),
[`vt_compilation()`](https://robustecologies.github.io/VerteTIME/reference/vt_compilation.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(vertetime); co <- vertetime
vt_publish(co, here::here("web-export","vertetime-v1.0"), overwrite = TRUE)
} # }
```
