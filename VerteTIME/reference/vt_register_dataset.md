# Register a new dataset into the compilation

Validates the YAML sidecar of a new dataset, runs
[`vt_ingest_dataset()`](https://robustecologies.github.io/VerteTIME/reference/vt_ingest_dataset.md)
and
[`vt_validate()`](https://robustecologies.github.io/VerteTIME/reference/vt_validate.md),
and (if both succeed) returns the freshly built `vt_dataset`. Designed
to be called by users who add a dataset after VerteTIME v1.0 ships,
following the workflow documented in the `vt-ingestion-workflow`
vignette.

## Usage

``` r
vt_register_dataset(
  id,
  data_raw = here::here("data-raw"),
  dry_run = FALSE,
  strict = TRUE
)
```

## Arguments

- id:

  Character scalar dataset identifier of the form `VT_NNN` matching a
  subfolder of `data_raw` and a YAML sidecar `data-raw/_yaml/<id>.yaml`.
  New datasets get the next available number above the current maximum.

- data_raw:

  Root of the source `data-raw/` directory. Defaults to
  `here::here("data-raw")`.

- dry_run:

  Logical, default `FALSE`. When `TRUE`, returns the validation `tibble`
  and the would-be `vt_dataset` without writing anything.

- strict:

  Logical, default `TRUE`. When `TRUE`, any non-empty
  [`vt_validate()`](https://robustecologies.github.io/VerteTIME/reference/vt_validate.md)
  output aborts.

## Value

A `vt_dataset` (when `dry_run = FALSE`) or a list with components
`dataset` and `validation` (when `dry_run = TRUE`).

## Details

The YAML sidecar must follow the schema in
`inst/templates/dataset_template.yaml`. Required fields are
`dataset_id`, `primary_reference_citation`, `primary_reference_doi`,
`primary_reference_kind`, the `sites:` block with at least one site, and
the year range. Optional fields are `taxa_columns`, `covariate_columns`,
`units`, `taxonomic_focus`, `taxonomy_authority`, `notes` and
`inherited_constraints`. Columns are auto-classified into species and
covariates by
[`classify_columns()`](https://robustecologies.github.io/VerteTIME/reference/classify_columns.md)
using the regex `vt.species_regex`; manual overrides go through
`taxa_columns` and `covariate_columns`.

Validation is split between automatic checks (structure, types, year
range, regex name match, coordinate plausibility, DOI uniqueness, no
duplicated `(site, species, year)`) and human-in-the-loop checks
(taxonomic-backbone confirmation, coordinates transcription from the
primary reference, source verification). The
`inst/templates/register-checklist.md` file is a printable pre-flight
checklist.

## References

Wilkinson, M. D., et al. (2016). *The FAIR Guiding Principles for
scientific data management and stewardship*. Scientific Data, 3, 160018.
[doi:10.1038/sdata.2016.18](https://doi.org/10.1038/sdata.2016.18) .

## See also

[`vt_ingest_dataset()`](https://robustecologies.github.io/VerteTIME/reference/vt_ingest_dataset.md),
[`vt_validate()`](https://robustecologies.github.io/VerteTIME/reference/vt_validate.md),
[`vt_ingest_all()`](https://robustecologies.github.io/VerteTIME/reference/vt_ingest_all.md),
[`classify_columns()`](https://robustecologies.github.io/VerteTIME/reference/classify_columns.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# First, populate data-raw/VT_001/VT_001.csv and data-raw/_yaml/VT_001.yaml
result <- vt_register_dataset("VT_001", dry_run = TRUE)
if (nrow(result$validation) == 0L) {
  d <- vt_register_dataset("VT_001")
  summary(d)
}
} # }
```
