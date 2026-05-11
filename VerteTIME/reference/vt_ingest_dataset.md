# Ingest a single dataset folder into the canonical relational schema

Maintainer-side build tool. Reads the wide community CSV under
`<data_raw>/<id>/<id>.csv` (or every `<id>*.csv` for multi-site studies)
and the per-dataset YAML sidecar under `<data_raw>/_yaml/<id>.yaml`,
classifies columns into species and covariates with
[`classify_columns()`](https://robustecologies.github.io/VerteTIME/reference/classify_columns.md),
pivots to long form with
[`vt_long()`](https://robustecologies.github.io/VerteTIME/reference/vt_long.md),
assembles the five canonical relational tables, and returns a
[`vt_dataset()`](https://robustecologies.github.io/VerteTIME/reference/vt_dataset.md)
object. The YAML sidecar is the single source of truth for site
coordinates, units, taxonomic focus and the primary reference citation.

This function operates on a maintainer-side ingestion tree that is
excluded from the package build. End users never call it; the
public-facing entry point is `data(vertetime)` followed by
[`vt_read()`](https://robustecologies.github.io/VerteTIME/reference/vt_read.md)
for single-dataset access. The function is exported for forks and
audits, not for the regular analytical workflow, and aborts with a
coherent message when `data_raw` does not exist on disk.

## Usage

``` r
vt_ingest_dataset(id, data_raw = NULL, yaml = NULL, require_yaml = TRUE)
```

## Arguments

- id:

  Character scalar dataset identifier matching a subfolder of
  `data_raw`.

- data_raw:

  Path to a maintainer-side ingestion tree. Absent in the installed
  package; the function aborts when the directory does not exist.

- yaml:

  Optional path to a YAML sidecar; defaults to
  `<data_raw>/_yaml/<id>.yaml`.

- require_yaml:

  Logical, default `TRUE`. When the sidecar is missing, error out
  (`TRUE`) or build a stub with placeholder fields (`FALSE`).

## Value

A `vt_dataset` carrying `datasets`, `sites`, `species`, `observations`,
`covariates`, `data_provenance` tables.

## Details

Multi-site studies (filenames `<id>-1.csv`, `<id>-2.csv`, ...) produce
one site per CSV; the sidecar YAML lists the per-site coordinates and
habitat under the `sites:` block. Coordinate plausibility is checked at
this stage. The `is_community_metric_eligible` flag is set to `TRUE`
when the dataset has at least `vt.community_min_spp` species and
`vt.community_min_yrs` observed years. The YAML schema is documented in
`system.file("templates", "dataset_template.yaml", package = "VerteTIME")`.

## References

McGill, B. J., et al. (2007). *Species abundance distributions: moving
beyond single prediction theories to integration within an ecological
framework*. Ecology Letters, 10(10), 995-1015.
[doi:10.1111/j.1461-0248.2007.01094.x](https://doi.org/10.1111/j.1461-0248.2007.01094.x)
.

## See also

[`vt_ingest_all()`](https://robustecologies.github.io/VerteTIME/reference/vt_ingest_all.md),
[`vt_register_dataset()`](https://robustecologies.github.io/VerteTIME/reference/vt_register_dataset.md),
[`vt_validate()`](https://robustecologies.github.io/VerteTIME/reference/vt_validate.md),
[`classify_columns()`](https://robustecologies.github.io/VerteTIME/reference/classify_columns.md),
[`vt_read()`](https://robustecologies.github.io/VerteTIME/reference/vt_read.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Maintainer-side, from a private ingestion tree:
d <- vt_ingest_dataset("VT_001", data_raw = "<path-to-ingestion-tree>")
summary(d)
} # }
```
