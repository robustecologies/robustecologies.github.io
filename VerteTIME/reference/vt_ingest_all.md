# Ingest the full compilation from a maintainer-side ingestion tree

Maintainer-side build tool. Walks every immediate subfolder of
`data_raw` whose name does not begin with an underscore (control folders
are skipped), runs
[`vt_ingest_dataset()`](https://robustecologies.github.io/VerteTIME/reference/vt_ingest_dataset.md)
on each, concatenates the resulting relational tables, drops
non-vertebrate species via the cached GBIF backbone, and remaps every
internal source identifier to the canonical `VT_NNN` namespace before
returning. The returned
[`vt_compilation()`](https://robustecologies.github.io/VerteTIME/reference/vt_compilation.md)
therefore never carries the internal source prefixes; the source-folder
names are private build metadata.

This function operates on a maintainer-side ingestion tree that is
excluded from the package build. End users do not call it; the
public-facing entry point is `data(vertetime)`. The function is exported
for forks and audits. Calling it from an installed package without
supplying a valid `data_raw` errors out by design.

## Usage

``` r
vt_ingest_all(data_raw = NULL, require_yaml = TRUE, verbose = FALSE)
```

## Arguments

- data_raw:

  Path to a maintainer-side ingestion tree. Absent in the installed
  package; the function aborts when the directory does not exist.

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
alphabetically over the surviving folder names. The
internal-to-canonical mapping is a private maintainer artefact and lives
outside the package build.

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
# Maintainer-side, from a private ingestion tree:
vertetime <- vt_ingest_all(data_raw = "<path-to-ingestion-tree>",
                           verbose  = TRUE)
summary(vertetime)

# End-user side, from an installed copy of VerteTIME:
data(vertetime)
summary(vertetime)
} # }
```
