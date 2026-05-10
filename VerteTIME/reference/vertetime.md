# VerteTIME canonical compilation

The full VerteTIME compilation as a
[`vt_compilation()`](https://robustecologies.github.io/VerteTIME/reference/vt_compilation.md)
S3 object, with every dataset, site, species, observation, covariate row
and provenance entry already in the canonical `VT_NNN` namespace. This
is the entry point for the analytical workflow on the installed package:
end users do not call
[`vt_ingest_all()`](https://robustecologies.github.io/VerteTIME/reference/vt_ingest_all.md),
they `data(vertetime)`.

The symbol is named `vertetime` (not `vt_compilation`) because
[`vt_compilation()`](https://robustecologies.github.io/VerteTIME/reference/vt_compilation.md)
is the S3 constructor; with `LazyData: true`, a data object of the same
name would shadow the function in the package namespace and break every
constructor call. The dataset-level metadata and the provenance audit
table are exposed as `vertetime$datasets` and
`vertetime$data_provenance` respectively, so they are not duplicated as
separate data objects.

## Usage

``` r
vertetime
```

## Format

A `vt_compilation` carrying six `data.table` slots:

- `datasets`:

  One row per dataset: `dataset_id`, taxonomic focus, year range,
  `n_species`, `n_years`, eligibility flags.

- `sites`:

  One row per spatially-distinct site within a dataset: coordinates,
  elevation, habitat.

- `species`:

  Slim taxonomic backbone with class, order and family resolved against
  the GBIF backbone.

- `observations`:

  Long-format abundance, one row per (`dataset_id`, `site_id`,
  `species_id`, `year`).

- `covariates`:

  Long-format environmental covariates, one row per (`dataset_id`,
  `site_id`, `year`, `covariate_name`).

- `data_provenance`:

  Audit table linking every dataset to its primary reference (citation,
  DOI, kind, re-curation date).

## Source

Independent re-curation by Sergio Pico and Pablo Almaraz from the
primary peer-reviewed and grey literature cited in `data_provenance`.
Secondary compilations (LPI, BioTIME, GPDD, RivFishTIME) are referenced
in the manuscript only as positional comparators, never as a provenance
chain.

## Details

The compilation is rebuilt by the maintainer via
`Rscript data-raw/build-data.R`, which calls
[`vt_ingest_all()`](https://robustecologies.github.io/VerteTIME/reference/vt_ingest_all.md)
over the source `data-raw/` tree (excluded from the package build),
drops non-vertebrate species, remaps every internal source identifier to
the canonical `VT_NNN` namespace, and persists the resulting object
here. The internal source-folder names are private build metadata and
never surface in this object.

## See also

[`vt_compilation()`](https://robustecologies.github.io/VerteTIME/reference/vt_compilation.md),
[vt_demo](https://robustecologies.github.io/VerteTIME/reference/vt_demo.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(vertetime)
summary(vertetime)
vt_alpha_diversity(vertetime)
} # }
```
