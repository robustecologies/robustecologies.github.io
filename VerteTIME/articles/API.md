# Technical architecture of the VerteTIME compilation framework

  

## Design philosophy

VerteTIME is a curated compilation of world vertebrate community time
series. A single relational schema feeds a uniform analysis surface that
operates identically over a one-dataset slice and over the full
multi-hundred-site assembly. The architecture is designed so that a
curator writes the data once and every downstream tool, validator,
metric, plot and exporter consumes it through the same contract. This
section lays out the four architectural principles that make this
possible.

### Schema-first, function-second

Every entry point in the package consumes one of three types: a
`vt_dataset` (a single re-curated study), a `vt_compilation` (a
collection of datasets sharing the relational schema), or one of the
long observation tables that those objects expose (`observations`,
`covariates`, `data_provenance`). The functions never hold their own
copy of the data; they read the contract and return a typed result. A
new metric written tomorrow only has to know how to walk
`(dataset_id, site_id, species_id, year, value)` to slot in next to the
existing alpha, beta, turnover and nestedness routines. This separation
has two consequences: every analysis tool is agnostic to where the time
series came from, and a family of plotters can be compared on a fixed
compilation without rewriting anything.

### Public vs internal separation

The maintainer-private ingestion tree is excluded from the package
build. It carries the source CSVs, the per-dataset PDFs and the YAML
sidecars that record coordinates, primary-reference DOI and unit class.
Everything that any downstream user touches is materialised through the
public namespace `VT_NNN`, an opaque sequential identifier the package
assigns at ingestion time. The translation is invisible: the build-time
orchestrator returns a compilation already in the public namespace, and
the internal-to-canonical mapping is a private maintainer artefact that
never reaches the package build. The same separation keeps the
manuscript, the vignettes, the manual, the published export tree, the
SQLite database and the Frictionless Data Package consistent: every
public artefact reads from the same canonical surface.

### Provenance rule

Every entry began with the identification of a candidate primary
reference and the verbatim transcription of its abundance table. We did
not extract from secondary biodiversity time-series compilations;
partial overlaps, where they exist, are recorded in the per-dataset
`partial_overlap_with` field as comparator metadata only and are
surfaced as a positional figure
([`vt_plot_database_timeline()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_database_timeline.md)).
The compilation copyright belongs to the authors. The audit table
`data_provenance` carries the primary-reference citation, DOI, kind
(peer-reviewed, grey-literature, monograph, report), curation date and
any inherited constraint specific to a primary source.

### Progressive analysis

A trajectory through the API runs ingest, summary, metric, plot, export.
Each step returns an S3 object equipped with
[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html) and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods so that
interactive exploration is uniform. Plots emit the same theme
([`vt_theme()`](https://robustecologies.github.io/VerteTIME/reference/vt_theme.md))
and use one of two canonical palettes
([`vt_palette()`](https://robustecologies.github.io/VerteTIME/reference/vt_palette.md)
for taxonomic focus,
[`vt_class_palette()`](https://robustecologies.github.io/VerteTIME/reference/vt_class_palette.md)
for GBIF class). Filtering with
[`vt_filter()`](https://robustecologies.github.io/VerteTIME/reference/vt_filter.md)
cascades referential integrity automatically; the result is again a
`vt_compilation` that downstream tools accept without modification.
Export writes the same five tables in four formats (CSV, Apache Parquet,
SQLite single-file, Frictionless Data Package) plus the LPD-style
wide-master series table; cross-format byte-identical SHA-256 hashes are
guaranteed modulo file modification times.

Layered architecture of the VerteTIME framework. Each row is a
subsystem; arrows show the data flow from raw sources through the
ingestion engine and validation layer down to the public artefacts. The
five-table relational schema is the contract that every downstream tool
reads.

  

## Architecture overview

The package is organised around five subsystems. The **ingestion layer**
([`vt_ingest_dataset()`](https://robustecologies.github.io/VerteTIME/reference/vt_ingest_dataset.md),
[`vt_ingest_all()`](https://robustecologies.github.io/VerteTIME/reference/vt_ingest_all.md),
[`classify_columns()`](https://robustecologies.github.io/VerteTIME/reference/classify_columns.md),
[`enrich_taxonomy()`](https://robustecologies.github.io/VerteTIME/reference/enrich_taxonomy.md),
[`vt_register_dataset()`](https://robustecologies.github.io/VerteTIME/reference/vt_register_dataset.md))
walks the source folders and YAML sidecars in a maintainer-private
ingestion tree, transparently absorbs the format quirks of the CSVs
(BOM, separator, header-case, NA tokens), splits species columns from
environmental covariates with a regex heuristic, merges the GBIF
backbone for higher-rank taxonomy, and writes the result into the
relational schema. The **validation layer**
([`vt_validate()`](https://robustecologies.github.io/VerteTIME/reference/vt_validate.md))
runs structural and semantic checks against the schema and reports
failures one row at a time. The **compilation contract**
([`vt_compilation()`](https://robustecologies.github.io/VerteTIME/reference/vt_compilation.md),
[`vt_dataset()`](https://robustecologies.github.io/VerteTIME/reference/vt_dataset.md),
[`vt_wide_master()`](https://robustecologies.github.io/VerteTIME/reference/vt_wide_master.md))
is the typed surface every other tool reads. The **analysis toolkit**
(alpha, beta, turnover, composition, aggregators, filtering) consumes
the long-format observation table directly and returns either a tibble
or one of the typed S3 outputs (`vt_diversity`, `vt_turnover`,
`vt_provenance`). The **visualisation engine** (~33 plot functions
across descriptive, spatial, temporal, community-pattern,
biodiversity-portrait and meta families) reads the same compilation and
emits ggplot2 objects with a uniform theme. The **publication
subsystem**
([`vt_export()`](https://robustecologies.github.io/VerteTIME/reference/vt_export.md),
[`vt_publish()`](https://robustecologies.github.io/VerteTIME/reference/vt_publish.md),
[`vt_wide_master()`](https://robustecologies.github.io/VerteTIME/reference/vt_wide_master.md))
materialises the public release in four formats with a SHA-256 manifest.

### The compilation contract

A `vt_compilation` is a list with six slots that downstream tools read.
The core fields are `datasets` (one row per re-curated study with year
range, species and covariate counts, dataset-level coordinates, unit
class, taxonomic focus, eligibility flag), `sites` (one row per
spatially distinct site within a dataset), `species` (taxonomic backbone
via GBIF: `class`, `order`, `family`, `is_vertebrate`), `observations`
(long, the workhorse with
`(dataset_id, site_id, species_id, year, value)`), `covariates` (long,
environmental covariates measured on the same temporal index), and
`data_provenance` (the audit table with primary-reference DOI, citation,
kind, curation date, comparative `partial_overlap_with` list).

The contract is what makes the framework composable.
[`vt_alpha_diversity()`](https://robustecologies.github.io/VerteTIME/reference/vt_alpha_diversity.md)
reads `observations` and writes per-`(dataset, site, year)` index
values.
[`vt_plot_world_map()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_world_map.md)
reads `sites` and `datasets` and projects a Robinson world.
[`vt_publish()`](https://robustecologies.github.io/VerteTIME/reference/vt_publish.md)
writes every table into CSV, Parquet and SQLite formats keyed by `*_id`
foreign keys.
[`vt_filter()`](https://robustecologies.github.io/VerteTIME/reference/vt_filter.md)
cascades a tidy-eval predicate from `datasets` through the FK chain,
returning a smaller `vt_compilation` that every downstream tool accepts
without ceremony. None of these tools manipulate the internals directly;
they ask the contract.

Five-table relational schema with audit table. Foreign keys use the
\*\_id convention consistently. The same schema is materialised
identically across CSV, Parquet, SQLite and the Frictionless Data
Package manifest.

  

### Subsystem summary

| Subsystem | Entry points | Purpose |
|----|----|----|
| Ingestion | [`vt_ingest_dataset()`](https://robustecologies.github.io/VerteTIME/reference/vt_ingest_dataset.md), [`vt_ingest_all()`](https://robustecologies.github.io/VerteTIME/reference/vt_ingest_all.md), [`vt_register_dataset()`](https://robustecologies.github.io/VerteTIME/reference/vt_register_dataset.md) | Build a compilation from a maintainer-private ingestion tree; absorb format quirks; merge taxonomy |
| Reading | [`vt_read_csv()`](https://robustecologies.github.io/VerteTIME/reference/vt_read_csv.md), [`vt_read()`](https://robustecologies.github.io/VerteTIME/reference/vt_read.md), [`vt_long()`](https://robustecologies.github.io/VerteTIME/reference/vt_long.md), [`vt_wide()`](https://robustecologies.github.io/VerteTIME/reference/vt_wide.md), [`classify_columns()`](https://robustecologies.github.io/VerteTIME/reference/classify_columns.md) | Low-level IO, pivot between wide and long, regex-driven species/covariate split |
| Validation | [`vt_validate()`](https://robustecologies.github.io/VerteTIME/reference/vt_validate.md), [`enrich_taxonomy()`](https://robustecologies.github.io/VerteTIME/reference/enrich_taxonomy.md) | Structural and semantic checks; GBIF-backbone class/order/family enrichment |
| Construction | [`vt_compilation()`](https://robustecologies.github.io/VerteTIME/reference/vt_compilation.md), [`vt_dataset()`](https://robustecologies.github.io/VerteTIME/reference/vt_dataset.md), [`vt_filter()`](https://robustecologies.github.io/VerteTIME/reference/vt_filter.md) | S3 constructors and tidy-eval filter with cascading FK integrity |
| Diversity | [`vt_alpha_diversity()`](https://robustecologies.github.io/VerteTIME/reference/vt_alpha_diversity.md), [`vt_evenness()`](https://robustecologies.github.io/VerteTIME/reference/vt_evenness.md), [`vt_series_summary()`](https://robustecologies.github.io/VerteTIME/reference/vt_series_summary.md) | Richness, Shannon, Simpson, Hill numbers, Chao1, Pielou, Simpson evenness, per-series summaries |
| Beta/turnover | [`vt_beta_diversity()`](https://robustecologies.github.io/VerteTIME/reference/vt_beta_diversity.md), [`vt_temporal_turnover()`](https://robustecologies.github.io/VerteTIME/reference/vt_temporal_turnover.md), [`vt_beta_partition()`](https://robustecologies.github.io/VerteTIME/reference/vt_beta_partition.md) | Bray-Curtis, Jaccard, Hellinger, Whittaker multiplicative, Baselga (2010) decomposition |
| Composition | [`vt_dominance_curve()`](https://robustecologies.github.io/VerteTIME/reference/vt_dominance_curve.md), [`vt_rank_abundance()`](https://robustecologies.github.io/VerteTIME/reference/vt_rank_abundance.md), [`vt_rarefaction()`](https://robustecologies.github.io/VerteTIME/reference/vt_rarefaction.md), [`vt_nestedness()`](https://robustecologies.github.io/VerteTIME/reference/vt_nestedness.md) | Whittaker / RAC / Hurlbert-Sanders rarefaction / Almeida-Neto NODF |
| Aggregators | [`vt_community_summary()`](https://robustecologies.github.io/VerteTIME/reference/vt_community_summary.md), [`vt_compilation_summary()`](https://robustecologies.github.io/VerteTIME/reference/vt_compilation_summary.md) | Per-site-year community summary; per-dataset summary |
| Visualisation | [`vt_theme()`](https://robustecologies.github.io/VerteTIME/reference/vt_theme.md), [`vt_palette()`](https://robustecologies.github.io/VerteTIME/reference/vt_palette.md), [`vt_class_palette()`](https://robustecologies.github.io/VerteTIME/reference/vt_class_palette.md), ~33 `vt_plot_*` functions | Project-wide ggplot theme, two palettes, six plot families |
| Export | [`vt_export()`](https://robustecologies.github.io/VerteTIME/reference/vt_export.md), [`vt_publish()`](https://robustecologies.github.io/VerteTIME/reference/vt_publish.md), [`vt_wide_master()`](https://robustecologies.github.io/VerteTIME/reference/vt_wide_master.md) | Single-format and full-tree publication; LPD-style wide master series table |
| Self-test | [`vt_check()`](https://robustecologies.github.io/VerteTIME/reference/vt_check.md) | One-call smoke test that exercises every S3 family against the shipped `vt_demo` |

  

### Interoperability

Every analytical function returns either a tibble or one of the typed S3
outputs (`vt_diversity`, `vt_turnover`, `vt_provenance`). The tibble
outputs slot directly into `dplyr`, `tidyr`, `ggplot2` workflows. The S3
outputs carry [`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html) and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods so that
interactive use is uniform. The long observation table is
`tsibble`-ready when keyed by `(dataset_id, site_id, species_id)` with
`year` as the index. The wide-master output of
[`vt_wide_master()`](https://robustecologies.github.io/VerteTIME/reference/vt_wide_master.md)
is a drop-in for any downstream tool that expects the LPD wide layout.

  

## The ingestion pipeline

The ingestion pipeline operates on a maintainer-private tree that the
package build excludes. End users never run it on the shipped package;
the public-facing entry point is `data(vertetime)` followed by
`vt_read(id)` for any single-dataset slice. Users who maintain a private
fork of the schema can drive the same pipeline through
[`vt_register_dataset()`](https://robustecologies.github.io/VerteTIME/reference/vt_register_dataset.md)
(see the registration vignette).

### File discovery

`vt_ingest_dataset(id, data_raw)` walks the per-dataset subfolder of
`data_raw` and discovers CSVs that match `^<id>(-[0-9]+)?\\.csv$`. A
folder with a single `<id>.csv` is treated as a single-site study;
folders with `<id>-1.csv`, `<id>-2.csv`, … are multi-site studies in
which each CSV becomes its own site. The companion `<id>-N.pdf` (when
present) is parsed for per-site coordinates; the dataset-level
`<id>.pdf` (when present) supplies metadata that all sites share. The
discovery logic also handles auxiliary files (compressed archives,
helper scripts, supplementary PDFs) gracefully by ignoring anything that
does not match the CSV pattern.

### YAML sidecar contract

Each dataset has a YAML sidecar that carries the human-curated metadata:
`dataset_id`, `country_iso3`, `realm`, `biome`, `ecoregion`,
`primary_reference_citation`, `primary_reference_doi`,
`primary_reference_kind`, `re_curation_date`, `inherited_constraints`,
`partial_overlap_with`, `units`, `taxonomic_focus`,
`taxonomy_authority`, optional `taxa_columns` and `covariate_columns`
overrides, a `sites:` block with at least one site (each with `site_id`,
`latitude_dd`, `longitude_dd`, `coord_precision_km`, optional `habitat`
and `elevation_m`), `year_min`, `year_max`, and free-text `notes`. The
canonical schema ships inside the installed package and is reachable
through
`system.file("templates", "dataset_template.yaml", package = "VerteTIME")`;
a printable pre-flight checklist sits next to it at
`system.file("templates", "register-checklist.md", package = "VerteTIME")`.

### Column classification

[`classify_columns()`](https://robustecologies.github.io/VerteTIME/reference/classify_columns.md)
is the regex heuristic that splits species columns from covariate
columns. The default pattern is `^[A-Z][a-z]+_[a-z][a-z0-9-]*$`
(capitalised genus, underscore, lowercase species epithet with optional
hyphen and digits). Columns that match the regex are taxa; columns that
do not (`year`, `rain`, `wet_rain`, `temperature`, `ndvi`) are
covariates. The YAML sidecar can override the heuristic when needed via
`taxa_columns:` and `covariate_columns:` lists. The split is
deterministic and re-evaluated on every ingest.

### Taxonomy enrichment

[`enrich_taxonomy()`](https://robustecologies.github.io/VerteTIME/reference/enrich_taxonomy.md)
merges the GBIF backbone (cached in the maintainer-private ingestion
tree) into the species table. Every species name resolves to a `class`,
`order`, `family` and an `is_vertebrate` flag. Non-vertebrate species
(insects, gastropods, plants and algae that some primary references
include alongside the vertebrate community of interest) are filtered out
automatically by
[`vt_ingest_all()`](https://robustecologies.github.io/VerteTIME/reference/vt_ingest_all.md);
datasets that lose all observations after the filter are dropped from
the compilation with a console message naming them.

### Public namespace

After taxonomy enrichment and filtering,
[`vt_ingest_all()`](https://robustecologies.github.io/VerteTIME/reference/vt_ingest_all.md)
rewrites every `dataset_id` and `site_id` to the canonical `VT_NNN`
namespace. The translation is a sorted alphabetical mapping over the
source-folder names and is stable across re-ingestions of the same
compilation. The translation is invisible to user code:
[`vt_ingest_all()`](https://robustecologies.github.io/VerteTIME/reference/vt_ingest_all.md)
returns a compilation that already carries `VT_NNN` identifiers; the
internal-to-canonical mapping is a private maintainer artefact outside
the package build. Every public artefact (the package’s shipped
`vertetime` and `vt_demo` data objects, the manuscript appendix, the
published export tree, the SQLite database, the Frictionless Data
Package, the pkgdown website) reads exclusively the `VT_NNN` namespace.

End-to-end ingestion pipeline. Yellow nodes are private inputs; green
nodes are private operations; blue nodes are the public surface.

  

## Diversity catalogue

### Alpha diversity

`vt_alpha_diversity(x, indices)` returns one row per
`(dataset_id, site_id, year, index)` triple. Supported indices are `S`
(richness), `H` (Shannon entropy in nats), `D` (Simpson’s
`1 - sum(p^2)`), `q0` / `q1` / `q2` (Hill numbers of order 0, 1, 2), and
`Chao1` (bias-corrected Chao 1984 estimator). All indices are computed
in pure R, with `vegan` only as an optional cross-validation dependency.
The Hill-number profile across `q in [0, 3]` is exposed via
[`vt_plot_hill_profile()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_hill_profile.md)
for per-dataset diagnostics.

`vt_evenness(x, kind)` returns Pielou `J = H / log(S)` or Simpson
evenness `(1 / (1 - D)) / S`. Both are undefined when `S < 2`; the
function returns `NA` rather than `0` or `Inf` in that regime.

### Beta diversity and temporal turnover

`vt_beta_diversity(x, method)` computes pairwise dissimilarity between
sites within each dataset (Bray-Curtis, Jaccard, Hellinger) plus a
Whittaker multiplicative scalar `gamma / mean(alpha) - 1` summarised at
dataset level. Single-site datasets return `NA` for the pairwise
matrices and a defined Whittaker scalar.

`vt_temporal_turnover(x, method, pair)` computes pairwise year-by-year
dissimilarity within each site for one of three pair patterns:
`consecutive` (year `t` vs year `t+1`), `first_vs_each` (`year_min` vs
every other year), or `all_pairs` (every distinct `t_a < t_b`).

`vt_beta_partition(x, family)` implements the Baselga (2010)
decomposition of pairwise temporal Sorensen dissimilarity into the
Simpson turnover component and the nestedness-resultant component. The
identity `beta_sor = beta_sim + beta_nes` holds to numerical precision
across the compilation. The `vt_turnover` S3 output carries
[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html) and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods;
[`summary.vt_turnover()`](https://robustecologies.github.io/VerteTIME/reference/vt_turnover.md)
dispatches on the columns present and emits a tidy tibble keyed by
metric.

### Composition

[`vt_dominance_curve()`](https://robustecologies.github.io/VerteTIME/reference/vt_dominance_curve.md),
[`vt_rank_abundance()`](https://robustecologies.github.io/VerteTIME/reference/vt_rank_abundance.md),
[`vt_rarefaction()`](https://robustecologies.github.io/VerteTIME/reference/vt_rarefaction.md)
and
[`vt_nestedness()`](https://robustecologies.github.io/VerteTIME/reference/vt_nestedness.md)
complete the community-structure surface. NODF (Almeida-Neto et
al. 2008) is computed on year-by-species presence-absence matrices per
`(dataset_id, site_id)`. Hurlbert-Sanders rarefaction returns the
expected richness at a target sample size, with `NA` when the year-site
total falls below the requested sample.

  

## Visualisation engine

The package ships ~33 plot functions organised in six families. All emit
ggplot2 objects with a shared theme
([`vt_theme()`](https://robustecologies.github.io/VerteTIME/reference/vt_theme.md))
and one of two canonical palettes
([`vt_palette()`](https://robustecologies.github.io/VerteTIME/reference/vt_palette.md)
for taxonomic focus,
[`vt_class_palette()`](https://robustecologies.github.io/VerteTIME/reference/vt_class_palette.md)
for GBIF class).

| Family | Functions (selected) |
|----|----|
| Theme and palette | [`vt_theme()`](https://robustecologies.github.io/VerteTIME/reference/vt_theme.md), [`vt_palette()`](https://robustecologies.github.io/VerteTIME/reference/vt_palette.md), [`vt_class_palette()`](https://robustecologies.github.io/VerteTIME/reference/vt_class_palette.md) |
| Descriptive (compilation) | [`vt_plot_focus_bars()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_focus_bars.md), [`vt_plot_class_bars()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_class_bars.md), [`vt_plot_class_treemap()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_class_treemap.md), [`vt_plot_top_species()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_top_species.md), [`vt_plot_top_genera()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_top_genera.md) |
| Descriptive (coverage) | [`vt_plot_span_hist()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_span_hist.md), [`vt_plot_focus_span_box()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_focus_span_box.md), [`vt_plot_span_coverage_hex()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_span_coverage_hex.md), [`vt_plot_mean_cv()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_mean_cv.md), [`vt_plot_alpha_hist()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_alpha_hist.md), [`vt_plot_nodf_hist()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_nodf_hist.md) |
| Spatial | [`vt_plot_world_map()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_world_map.md), [`vt_plot_world_hex()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_world_hex.md), [`vt_plot_world_richness()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_world_richness.md), [`vt_plot_lat_hist()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_lat_hist.md), [`vt_plot_country_bar()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_country_bar.md), [`vt_plot_realm_bar()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_realm_bar.md), [`vt_plot_focus_realm_heatmap()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_focus_realm_heatmap.md), [`vt_plot_alluvial()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_alluvial.md) |
| Temporal | [`vt_plot_lasagna()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_lasagna.md), [`vt_plot_active_area()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_active_area.md), [`vt_plot_onset_cumulative()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_onset_cumulative.md), [`vt_plot_year_ridges()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_year_ridges.md), [`vt_plot_temporal_sankey()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_temporal_sankey.md), [`vt_plot_miss_by_year()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_miss_by_year.md), [`vt_plot_miss_class_heatmap()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_miss_class_heatmap.md) |
| Community pattern | [`vt_plot_whittaker()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_whittaker.md), [`vt_plot_rac()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_rac.md), [`vt_plot_alpha_timeline()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_alpha_timeline.md), [`vt_plot_beta_heatmap()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_beta_heatmap.md), [`vt_plot_metric_ridges()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_metric_ridges.md), [`vt_plot_hill_profile()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_hill_profile.md), [`vt_plot_community_series()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_community_series.md), [`vt_plot_community_spotlight()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_community_spotlight.md) |
| Meta (positional comparator) | [`vt_plot_database_timeline()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_database_timeline.md) (with [`vt_compilations()`](https://robustecologies.github.io/VerteTIME/reference/vt_compilations.md) lookup) |

The flagship per-dataset visualisation is
[`vt_plot_community_spotlight()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_community_spotlight.md),
which composes a community-time-series view (one line per species,
viridis-turbo contrasted colours, log-y by default) with the
Robinson-projected world map of the dataset’s site(s) using `patchwork`
for layout.

The two S3 plot methods of the compilation contract
([`plot.vt_dataset()`](https://robustecologies.github.io/VerteTIME/reference/vt_dataset.md),
[`plot.vt_compilation()`](https://robustecologies.github.io/VerteTIME/reference/vt_compilation.md))
accept a `type =` argument that dispatches to the corresponding
`vt_plot_*` family member. The legacy `kind =` alias is retained
indefinitely; new code should prefer `type =`.

  

## Publication subsystem

`vt_export(x, format, path)` writes a single format.
`vt_publish(x, root, overwrite)` orchestrates all four formats plus the
README, the licence, the citation file and a SHA-256 manifest at the
root. Re-running on an unchanged compilation reproduces byte-identical
SHA-256 hashes for every emitted file (modulo file modification times).

Publication subsystem. The same compilation feeds four format writers
and an integrity manifest. The Zenodo deposit step is manual.

The `master.csv` is the long-form denormalised join (`observations`
joined with `datasets`, `sites` and `species`); `master_wide.csv` is the
LPD-style wide table where each row is one series (`dataset_id`,
`site_id`, `species_id`) and the trailing columns are calendar years
from `year_min` to `year_max` populated with abundance values (or `NA`
where unobserved). The wide layout is a drop-in for any tool that
expects the historical wide format used by population-time-series
compilations.

  

## Validation contract

`vt_validate(x, level)` runs structural and semantic checks against the
schema. At the `dataset` level: column 1 is `year` (integer, in
`[vt.year_min, vt.year_max]`), `species_id` matches the regex, `value`
is non-negative or `NA`, no duplicated `(site_id, species_id, year)`
triples, coordinates are in valid ranges, `(0, 0)` is flagged as a
placeholder. At the `compilation` level: `dataset_id` uniqueness,
`primary_reference_doi` uniqueness (warning only because some studies
legitimately share a primary reference across multiple sites), and
cross-dataset `species_id` consistency. The function returns a tibble
with one row per failed check; an empty tibble means a clean bill of
health.

The dataset-registration workflow is
`vt_register_dataset(id, data_raw, dry_run)`. With `dry_run = TRUE` and
an absent sidecar the function scaffolds a YAML template into the user’s
private tree; with `dry_run = TRUE` and a complete sidecar the function
returns the would-be `vt_dataset` and the validation tibble without
writing anything; with `dry_run = FALSE` it appends the dataset to the
user’s private tree and writes a register-log row.

  

## Self-test

`vt_check(verbose, plots)` is the canonical smoke test for the installed
package. It loads `data(vt_demo)` into a private environment, exercises
one representative call per analytical family, exercises every S3
`print` / `summary` / `plot` on the resulting objects, exercises
[`vt_publish()`](https://robustecologies.github.io/VerteTIME/reference/vt_publish.md)
and
[`vt_export()`](https://robustecologies.github.io/VerteTIME/reference/vt_export.md)
against a [`tempfile()`](https://rdrr.io/r/base/tempfile.html), and
returns a tibble whose `status` column is `"ok"` for every row when the
package is healthy. The `attr(result, "ok")` flag is the single boolean
a reader can branch on.

  

## Changelog

See the
[Changelog](https://robustecologies.github.io/VerteTIME/news/index.md)
for release notes. `NEWS.md` is the canonical user-visible record of
every public release; per-release entries cover ingestion, validation,
metric, plot and export changes that downstream users would notice.
