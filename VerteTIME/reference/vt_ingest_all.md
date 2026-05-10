# Ingest the full compilation from `data-raw/`

Maintainer-side build tool. Walks every immediate subfolder of
`data_raw` whose name does not begin with an underscore (excluding
`_yaml/`, `_staging/` and similar control folders), runs
[`vt_ingest_dataset()`](https://robustecologies.github.io/VerteTIME/reference/vt_ingest_dataset.md)
on each, concatenates the resulting relational tables, drops
non-vertebrate species via the GBIF taxonomy cache, and remaps every
internal source identifier to the canonical `VT_NNN` namespace before
returning. The returned
[`vt_compilation()`](https://robustecologies.github.io/VerteTIME/reference/vt_compilation.md)
therefore never carries the internal source prefixes; the source-folder
names are private build metadata.

This function is a build-time tool used by `data-raw/build-data.R` to
assemble the shipped `vt_compilation` object. End users do not call it;
the public-facing entry point is `data(vertetime)`. The function is
exported for forks and audits, not for the regular analytical workflow.
Calling `vt_ingest_all()` from an installed copy of the package without
access to the original `data-raw/` tree errors out by design, because
that tree is excluded from the package build (`.Rbuildignore`).

## Usage

``` r
vt_ingest_all(
  data_raw = here::here("data-raw"),
  require_yaml = TRUE,
  verbose = FALSE
)
```

## Arguments

- data_raw:

  Root of the source `data-raw/` directory. Defaults to
  `here::here("data-raw")`. Absent in the installed package.

- require_yaml:

  Logical, default `TRUE`. Forwarded to
  [`vt_ingest_dataset()`](https://robustecologies.github.io/VerteTIME/reference/vt_ingest_dataset.md).

- verbose:

  Logical, default `FALSE`. When `TRUE`, emits aggregate progress to
  `stderr` (number of folders walked, number of datasets ingested,
  number of canonical `VT_NNN` identifiers issued). Per-dataset progress
  is never emitted with the source-folder name; the internal prefixes
  are private build metadata.

## Value

A `vt_compilation` object whose `dataset_id` and `site_id` columns live
entirely in the canonical `VT_NNN` namespace.

## Details

Datasets that fail validation contribute one row to a warning summary at
the end of the run; the rest of the compilation is built unaffected. Set
`vt.staging_dir` to a path to enable per-dataset Parquet caching of the
long-format observation block, which speeds up re-runs substantially
when only a handful of source CSVs have changed.

Internal-to-canonical mapping. The compilation is constructed first with
the raw source-folder identifiers (an internal index of incorporation
order, never a marker of provenance), and only after the relational
tables are assembled, validated and stripped of non-vertebrate species
are those identifiers replaced by the canonical `VT_NNN` namespace,
alphabetically over the surviving folder names. The maintainer-side
mapping is persisted as a TSV under `MEMORY/id_map.tsv` for audit; that
file is private and lives outside the package build.

## References

Magurran, A. E., & McGill, B. J. (2011). *Biological Diversity:
Frontiers in Measurement and Assessment*. Oxford University Press. ISBN
9780199580675.

## See also

[`vt_ingest_dataset()`](https://robustecologies.github.io/VerteTIME/reference/vt_ingest_dataset.md),
[`vt_register_dataset()`](https://robustecologies.github.io/VerteTIME/reference/vt_register_dataset.md),
[`vt_validate()`](https://robustecologies.github.io/VerteTIME/reference/vt_validate.md),
[`vt_publish()`](https://robustecologies.github.io/VerteTIME/reference/vt_publish.md),
[`vt_compilation()`](https://robustecologies.github.io/VerteTIME/reference/vt_compilation.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Maintainer-side, from the source repo:
vertetime <- vt_ingest_all(verbose = TRUE)
summary(vertetime)
usethis::use_data(vertetime, overwrite = TRUE)

# End-user side, from an installed copy of VerteTIME:
data(vertetime)
summary(vertetime)
} # }
```
