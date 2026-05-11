# Registering a private dataset against VerteTIME

This vignette is for users who want to extend the shipped `vertetime`
compilation with one or more private datasets of their own and run the
full VerteTIME analysis surface against the combined corpus. The
workflow is fully end-user-callable; nothing here depends on the
maintainer-side build tree. The canonical schema and the printable
checklist ship inside the installed package and are reached through
`system.file("templates", "dataset_template.yaml", package = "VerteTIME")`
and
`system.file("templates", "register-checklist.md", package = "VerteTIME")`.

If you only consume the shipped compilation, you do not need this
vignette; `data(vertetime)` and `vt_read(id)` cover the analytical
workflow without writing anywhere.

  

## When to use the scaffolder

[`vt_register_dataset()`](https://robustecologies.github.io/VerteTIME/reference/vt_register_dataset.md)
is the user-facing scaffolder. It does two things. First, when the YAML
sidecar for a new identifier is absent, it copies the package template
into the user’s private tree so the user can edit it. Second, when both
the sidecar and the wide-format CSV are present, it runs the same
ingestion-and-validation pipeline that produced the shipped `vertetime`
compilation and returns a `vt_dataset` that the user can pass to every
analytical and visualisation function.

The function refuses to touch the maintainer-side ingestion tree; it
operates exclusively against the user-supplied `data_raw` argument. The
package writes inside that path only; it never writes inside the
installed library.

  

## Step 1. Pick a private ingestion path

Choose any folder on your machine to host the private extension tree.
The folder is created on demand:

``` r

tree <- file.path(Sys.getenv("HOME"), "vertetime-extension")
```

Pick a fresh dataset identifier of the form `VT_NNN` where `NNN` is
unique within your private tree. The canonical `VT_NNN` namespace of the
shipped compilation occupies a contiguous range; your private
identifiers can collide neither with that range nor with each other.

  

## Step 2. Scaffold the YAML sidecar

The first call writes the YAML template into the right location:

``` r

library(VerteTIME)
res <- vt_register_dataset("VT_999", data_raw = tree, dry_run = TRUE)
res$scaffold       # <tree>/_yaml/VT_999.yaml
res$csv_dir        # <tree>/VT_999/
res$message
```

The function returns a list naming the YAML path it wrote and the
per-dataset folder where you should place the CSV. Open the YAML in an
editor and fill in the required fields (marked `[REQUIRED]` in the
template): `dataset_id`, `primary_reference_citation`,
`primary_reference_doi`, `primary_reference_kind`, the `sites:` block
with at least one site (with non-placeholder coordinates), `year_min`
and `year_max`. The optional fields can stay empty.

  

## Step 3. Drop the wide-format CSV in place

VerteTIME’s wide CSV format is one row per calendar year. Column 1 is
`year` (integer); subsequent columns are species in `Genus_species` form
(hyphenated epithets allowed) and any environmental covariates. The
field separator is `;` and the encoding is UTF-8. Place the file at
`<tree>/VT_999/VT_999.csv`. Multi-site studies use the suffix convention
`VT_999-1.csv`, `VT_999-2.csv`, etc.

If your column names do not follow the underscore-binomial convention,
set `taxa_columns:` and `covariate_columns:` in the YAML to override the
regex classifier.

  

## Step 4. Dry-run the registration

Once both the YAML and the CSV are populated, repeat the dry run to
validate without writing the audit log:

``` r

res <- vt_register_dataset("VT_999", data_raw = tree, dry_run = TRUE)
res$validation
```

The returned `validation` tibble has zero rows on a clean dataset;
otherwise each row identifies a specific check that failed (year out of
range, negative abundance, duplicate `(site, species, year)`, coordinate
out of range, species name not matching the regex). Address every issue
before running with `dry_run = FALSE`.

  

## Step 5. Register the dataset

``` r

d <- vt_register_dataset("VT_999", data_raw = tree)
summary(d)
plot(d, type = "whittaker")
```

When validation succeeds, the function appends one row to
`<tree>/_yaml/_register_log.csv` with a timestamp, the dataset’s row
counts, and a SHA-256 of the YAML sidecar, and returns the freshly built
`vt_dataset` ready for analysis. You can call every analytical function
on `d` directly (`vt_alpha_diversity(d)`, `vt_temporal_turnover(d)`,
`vt_plot_community_spotlight(d)`, etc.).

  

## Step 6. Combine with the shipped compilation

To analyse your private datasets jointly with `vertetime`, concatenate
the relational tables of the new `vt_dataset` with the shipped
compilation:

``` r

data(vertetime)
co <- vt_compilation(
  datasets        = rbind(vertetime$datasets,        d$datasets,        fill = TRUE),
  sites           = rbind(vertetime$sites,           d$sites,           fill = TRUE),
  species         = unique(rbind(vertetime$species,  d$species,         fill = TRUE),
                           by = "species_id"),
  observations    = rbind(vertetime$observations,    d$observations,    fill = TRUE),
  covariates      = rbind(vertetime$covariates,      d$covariates,      fill = TRUE),
  data_provenance = rbind(vertetime$data_provenance, d$data_provenance, fill = TRUE)
)
summary(co)
```

The merged `vt_compilation` carries both namespaces side by side and
every plot or metric works against it without further ceremony. The
shipped object is not modified; the combined `co` lives in the current R
session only.

  

## What is automatic and what is human-in-the-loop

Automatic checks: structure (column 1 is `year`, separators, types,
range), regex name match, coordinate plausibility, DOI uniqueness, no
duplicated `(site, species, year)`, non-negative abundance.

Human-in-the-loop: taxonomic-backbone confirmation (some hyphenated
epithets and old combinations need manual mapping), coordinate
transcription from the primary reference, primary-reference verification
(open the DOI, confirm authors and year), and the licensing decision
when the source carries non-default terms. The pre-flight checklist at
`system.file("templates", "register-checklist.md", package = "VerteTIME")`
lists these explicitly.

  

## Limitations and mitigations

The auto-classifier of species versus covariate columns relies on the
`^[A-Z][a-z]+_[a-z][a-z0-9-]*$` regex. If your dataset uses non-binomial
taxa labels (genus-level, hybrid notation, or vernacular abbreviations),
populate `taxa_columns` and `covariate_columns` explicitly in the YAML.
The validation pipeline does not check whether the values *should* be
counts versus densities versus indices; that is recorded in `unit_class`
and is the user’s responsibility.

The GBIF backbone enrichment used by the shipped compilation is not
reproducible on a private tree without a cached taxonomy TSV.
[`enrich_taxonomy()`](https://robustecologies.github.io/VerteTIME/reference/enrich_taxonomy.md)
is offered as a no-op when the cache is absent, so the `class`, `order`
and `family` columns for private species stay `NA`. Joining against the
shipped `vertetime$species` table covers the species already in the
canonical compilation; the rest can be enriched manually or left
unresolved with no impact on the abundance-based metrics.

If you want to contribute the dataset back to VerteTIME proper, open an
issue or pull request at the project’s GitHub repository with the YAML,
the CSV and the primary-reference DOI. The maintainers run the same
[`vt_register_dataset()`](https://robustecologies.github.io/VerteTIME/reference/vt_register_dataset.md)
workflow on their side, against the maintainer-private ingestion tree
that produces the next shipped release.
