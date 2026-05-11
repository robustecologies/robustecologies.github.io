# Register a new dataset in a private VerteTIME extension tree

User-facing scaffolder. Given a private ingestion tree at `data_raw`,
registers a new dataset by first copying the YAML sidecar template when
absent, and then running the ingestion-and-validation pipeline once the
user has edited the template and dropped a wide-format CSV into the
expected location. Users who only consume the shipped compilation do not
need this function; it exists for users who maintain a private fork that
extends VerteTIME with additional series.

## Usage

``` r
vt_register_dataset(id, data_raw, dry_run = FALSE, strict = TRUE)
```

## Arguments

- id:

  Character scalar dataset identifier of the form `VT_NNN`. New datasets
  pick the next available number above the current maximum.

- data_raw:

  Required path to the user's private ingestion tree. The tree must
  contain (or will be initialised with) a `_yaml/` subfolder for
  sidecars and a per-dataset subfolder `<id>/` for the CSV. The path is
  created if missing.

- dry_run:

  Logical, default `FALSE`. When `TRUE`, and when the YAML sidecar at
  `<data_raw>/_yaml/<id>.yaml` is absent, the function initialises the
  tree, copies the package-shipped YAML template into the sidecar
  location, and returns a `list(scaffold = <yaml-path>, message = ...)`
  for the user to edit. When the sidecar already exists,
  `dry_run = TRUE` returns the would-be `vt_dataset` and its validation
  tibble without writing the audit log.

- strict:

  Logical, default `TRUE`. When `TRUE`, any non-empty
  [`vt_validate()`](https://robustecologies.github.io/VerteTIME/reference/vt_validate.md)
  output aborts.

## Value

When `dry_run = TRUE` and the sidecar is missing, a `list` with
`scaffold` (the written YAML path) and `message` (next-step
instruction). When `dry_run = TRUE` and the sidecar exists, a `list`
with `dataset` (the would-be `vt_dataset`) and `validation` (the issues
tibble). When `dry_run = FALSE`, the freshly built
[`vt_dataset()`](https://robustecologies.github.io/VerteTIME/reference/vt_dataset.md).

## Details

The YAML schema and the per-dataset checklist ship inside the installed
package and are reachable via
`system.file("templates", "dataset_template.yaml", package = "VerteTIME")`
and
`system.file("templates", "register-checklist.md", package = "VerteTIME")`.
Required YAML fields are `dataset_id`, `primary_reference_citation`,
`primary_reference_doi`, `primary_reference_kind`, the `sites:` block
with at least one site, and the year range. Optional fields are
`taxa_columns`, `covariate_columns`, `units`, `taxonomic_focus`,
`taxonomy_authority`, `notes` and `inherited_constraints`. Columns are
auto-classified into species and covariates by
[`classify_columns()`](https://robustecologies.github.io/VerteTIME/reference/classify_columns.md)
using the regex from `vt.species_regex`; manual overrides go through
`taxa_columns` and `covariate_columns`.

Validation is split between automatic checks (structure, types, year
range, regex name match, coordinate plausibility, DOI uniqueness, no
duplicated `(site, species, year)`) and human-in-the-loop checks
(taxonomic-backbone confirmation, coordinates transcription from the
primary reference, source verification). The audit log appended to
`<data_raw>/_yaml/_register_log.csv` records the timestamp, the SHA-256
of the YAML file and the validation outcome.

## References

Wilkinson, M. D., et al. (2016). *The FAIR Guiding Principles for
scientific data management and stewardship*. Scientific Data, 3, 160018.
[doi:10.1038/sdata.2016.18](https://doi.org/10.1038/sdata.2016.18) .

## See also

[`vt_ingest_dataset()`](https://robustecologies.github.io/VerteTIME/reference/vt_ingest_dataset.md),
[`vt_validate()`](https://robustecologies.github.io/VerteTIME/reference/vt_validate.md),
[`vt_ingest_all()`](https://robustecologies.github.io/VerteTIME/reference/vt_ingest_all.md),
[`classify_columns()`](https://robustecologies.github.io/VerteTIME/reference/classify_columns.md),
[`vt_read()`](https://robustecologies.github.io/VerteTIME/reference/vt_read.md)

## Examples

``` r
if (FALSE) { # \dontrun{
tree <- file.path(tempdir(), "my-vertetime-extension")
# First call: scaffold the YAML sidecar.
vt_register_dataset("VT_999", data_raw = tree, dry_run = TRUE)
# User edits the YAML and drops VT_999/VT_999.csv into the tree.
# Second call: dry-run validation.
vt_register_dataset("VT_999", data_raw = tree, dry_run = TRUE)
# Third call: commit the dataset.
d <- vt_register_dataset("VT_999", data_raw = tree)
summary(d)
} # }
```
