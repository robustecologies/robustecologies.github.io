# Adding a new dataset to VerteTIME

This vignette walks through the workflow for registering a new community
time-series dataset into the VerteTIME compilation. Treat it as the
canonical companion of `inst/templates/dataset_template.yaml` and
`inst/templates/register-checklist.md`.

  

## Where this vignette sits in the package

VerteTIME has a three-layer architecture. The first layer is the
**end-user workflow**, in which an analyst loads the shipped compilation
with `data(vertetime)` and runs the diversity, turnover and
visualisation functions on it. The second layer is the **audit and fork
workflow**, in which a third party clones the source repository, opens
the source tree at `data-raw/`, and re-executes
[`vt_ingest_all()`](https://robustecologies.github.io/VerteTIME/reference/vt_ingest_all.md)
to reproduce the compilation from the primary references. The third
layer, documented here, is the **maintainer workflow**, in which Sergio
Pico and Pablo Almaraz add a new dataset to the source tree and rebuild
the shipped objects via `Rscript data-raw/build-data.R`.

The source `data-raw/` tree is excluded from the package build
(`.Rbuildignore`); an installed copy of VerteTIME does not carry it.
[`vt_ingest_all()`](https://robustecologies.github.io/VerteTIME/reference/vt_ingest_all.md)
is exported for completeness and audit, but the function aborts with an
explanatory error when it cannot find the source tree, and it never
emits the internal source-folder names to the console: progress is
reported as aggregate counts, and the canonical `VT_NNN` namespace is
the only namespace that ever leaves the build script.

  

## Step 1. Identify the primary reference

VerteTIME is a compilation of independent re-curations from
peer-reviewed primary literature and grey literature. Before adding a
dataset to the compilation, locate the original peer-reviewed paper,
monograph, technical report or grey-literature source that publishes the
multi-species annual abundance time series, read it in full, and confirm
that the data are either published openly or that you have written
authorisation to redistribute them under CC-BY 4.0.

If the candidate dataset is also present in a secondary biodiversity
time-series compilation (LPI, BioTIME, GPDD, RivFishTIME, PREDICTS, BBS,
PECBMS, CBC, eBird Trends, TEAM/Wildlife Insights), do not extract from
that compilation. Go to the primary reference and re-extract. The fact
of the partial overlap is recorded comparatively in the dataset’s
`partial_overlap_with` field; it is never used to derive provenance.

  

## Step 2. Place the source files

Adopt a new dataset identifier of the form `VT_NNN`, where `NNN` is the
next available sequential number above the current compilation maximum.
Create the folder and place the CSV plus the primary-reference PDF:

    data-raw/VT_001/VT_001.csv
    data-raw/VT_001/<primary-reference>.pdf

The CSV is wide format: column 1 is `year` (integer); subsequent columns
are species in `Genus_species` form (hyphenated epithets allowed) and
any environmental covariates. The field separator is `;` and the
encoding is UTF-8.

  

## Step 3. Fill the YAML sidecar

Copy the template:

``` r

file.copy(
  from = system.file("templates", "dataset_template.yaml", package = "VerteTIME"),
  to   = here::here("data-raw","_yaml","VT_001.yaml")
)
```

Open the copy and edit the required fields (marked `[REQUIRED]` in the
template). The minimum complete set is `dataset_id`,
`primary_reference_citation`, `primary_reference_doi`,
`primary_reference_kind`, the `sites:` block with at least one site
(with non-placeholder coordinates), and `year_min` and `year_max`.
Optional fields can be left empty.

  

## Step 4. Dry-run the registration

``` r

res <- vt_register_dataset("VT_001", dry_run = TRUE)
res$validation
```

`dry_run = TRUE` runs the ingestion and validation pipeline without
writing anything. The returned `validation` tibble has zero rows on a
clean dataset; otherwise each row identifies a specific check that
failed (year out of range, negative abundance, duplicate
`(site, species, year)`, coordinate out of range, species name not
matching the regex). Address every issue before running with
`dry_run = FALSE`.

  

## Step 5. Register the dataset

``` r

d <- vt_register_dataset("VT_001")
summary(d)
```

When validation succeeds, the function appends one row to
`data-raw/_yaml/_register_log.csv` with a timestamp, the dataset’s row
counts, and a SHA-256 of the YAML sidecar, and returns the freshly built
`vt_dataset` ready for inclusion in the compilation.

  

## Step 6. Rebuild the shipped data objects and re-publish

The build script is the single legitimate entry point that ingests the
source tree and writes the two shipped data objects (`vertetime` and
`vt_demo`) to `data/*.rda`. The internal source-folder identifiers are
translated to the canonical `VT_NNN` namespace inside the script and
never persisted. The `vertetime` object holds the full canonical
compilation (its `$datasets` slot is the dataset-level metadata and its
`$data_provenance` slot is the audit table; neither is shipped as a
separate data object).

``` r

source(here::here("data-raw", "build-data.R"))
```

After the build, regenerate the multi-format public release:

``` r

data(vertetime)
vt_publish(vertetime, here::here("web-export","vertetime-v1.0"),
           overwrite = TRUE)
```

[`vt_publish()`](https://robustecologies.github.io/VerteTIME/reference/vt_publish.md)
regenerates the four-format public release (CSV, Parquet, SQLite,
datapackage) plus the `README.md`, `LICENSE.md`, `CITATION.cff`,
`data_provenance.csv` and `CHECKSUMS.sha256` files at the root. The
Zenodo deposit step is manual: upload the new
`web-export/vertetime-v1.0` tree to a new Zenodo version, copy the
resulting DOI back into `CITATION.cff`, and re-render the manuscript.

  

## What is automatic and what is human-in-the-loop

Automatic checks: structure (column 1 is `year`, separators, types,
range), regex name match, coordinate plausibility, DOI uniqueness, no
duplicated `(site, species, year)`, non-negative abundance.

Human-in-the-loop: taxonomic-backbone confirmation (some hyphenated
epithets and old combinations need manual mapping), coordinate
transcription from the primary reference, primary-reference verification
(open the DOI, confirm authors and year), and the licensing decision
when the source carries non-default terms. The pre-flight checklist in
`inst/templates/register-checklist.md` lists these explicitly.

  

## Limitations and mitigations

The auto-classifier of species versus covariate columns relies on the
`^[A-Z][a-z]+_[a-z][a-z0-9-]*$` regex. If your dataset uses non-binomial
taxa labels (genus-level, hybrid notation, or vernacular abbreviations),
populate `taxa_columns` and `covariate_columns` explicitly in the YAML.
The validation pipeline does not check whether the values *should* be
counts versus densities versus indices; that is recorded in `unit_class`
and is the user’s responsibility.
