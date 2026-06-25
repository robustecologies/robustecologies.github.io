# Changelog

## VerteTIME 1.0.1

### Bug fixes

- [`vt_plot_class_treemap()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_class_treemap.md)
  now includes reptile taxa. The previous version filtered species on a
  non-missing `order` field, but the GBIF backbone promotes the
  historical reptile orders `Squamata`, `Testudines` and `Crocodylia` to
  class rank and returns an empty `order` for them, so every reptile
  species was silently dropped from the figure. Species are now retained
  on `class` and `family`, and the missing `order` is coalesced to the
  `class` so the class -\> order -\> family nesting still renders.
- [`vt_rarefaction()`](https://robustecologies.github.io/VerteTIME/reference/vt_rarefaction.md)
  now returns the rarefied richness in `S_rarefied`; a key-separator
  mismatch previously left the column `NA` for every row.
- [`vt_series_summary()`](https://robustecologies.github.io/VerteTIME/reference/vt_series_summary.md)
  now computes `coverage_pct`, `longest_internal_gap` and
  `n_internal_gaps`; these were returned as `NA` regardless of the
  series.
- `vt_plot_alpha_timeline(colour_by = "realm")` now colours traces by
  realm instead of rendering every trace grey.
- [`vt_plot_focus_span_box()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_focus_span_box.md)
  no longer converts the input compilation’s `taxonomic_focus` column to
  a factor by reference; the caller’s object is left unchanged.
- [`vt_plot_lat_hist()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_lat_hist.md)
  and
  [`vt_plot_world_richness()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_world_richness.md)
  no longer drop sites lying exactly on the equator, and
  [`vt_plot_lat_hist()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_lat_hist.md)
  retains sites with an unrecorded realm under an `unknown` group
  instead of discarding them.
- [`vt_plot_temporal_sankey()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_temporal_sankey.md)
  now counts a dataset as active in the decade of its final year when
  that year falls on a decade boundary.
- [`vt_export()`](https://robustecologies.github.io/VerteTIME/reference/vt_export.md)
  and
  [`vt_publish()`](https://robustecologies.github.io/VerteTIME/reference/vt_publish.md)
  now store `re_curation_date` as ISO text in the SQLite output
  (previously a numeric days-since-epoch), write relative paths in
  `CHECKSUMS.sha256` (previously absolute when the root contained regex
  metacharacters), and list the `master_wide` resource in
  `datapackage.json`.
- [`vt_check()`](https://robustecologies.github.io/VerteTIME/reference/vt_check.md)
  now records an assertion that returns `FALSE` as an error rather than
  reporting it as `ok`.
- [`vt_beta_partition()`](https://robustecologies.github.io/VerteTIME/reference/vt_beta_partition.md)
  no longer errors on a site whose every year-pair is empty.
- The eligibility flag `is_community_metric_eligible` in the shipped
  `vertetime` compilation is corrected: it is now recomputed after
  non-vertebrate species are removed, so a dataset reduced to a single
  retained species is no longer flagged community-metric eligible.
- `vt_export(format = "sqlite")` now materialises the documented
  relational schema: primary keys, enforced foreign keys, surrogate
  auto-increment ids on `observations` and `covariates`, and a
  `value_is_imputed` flag. Previously the SQLite tables were
  unconstrained and lacked these columns.
- `vt_export(format = "datapackage")` now declares a field schema (with
  types) for the `data_provenance` resource, so an all-empty
  `partial_overlap_with` column round-trips as a string rather than
  being inferred as logical.
- [`vt_read()`](https://robustecologies.github.io/VerteTIME/reference/vt_read.md)
  now rejects a non-existent file path with a clear error instead of
  silently treating it as a dataset_id;
  [`vt_long()`](https://robustecologies.github.io/VerteTIME/reference/vt_long.md)
  now errors clearly when the wide table has no `year` column.
- [`vt_alpha_diversity()`](https://robustecologies.github.io/VerteTIME/reference/vt_alpha_diversity.md)
  now warns once when `Chao1` is requested on non-integer abundances,
  which the estimator floors and which can distort the
  singleton/doubleton counts.
- [`vt_plot_span_coverage_hex()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_span_coverage_hex.md)
  now averages the per-series coverage proportion within each dataset;
  the previous form divided per-series observation counts by the
  dataset-level span, mixing two aggregation scopes.
- Dataset registration now derives canonical site identifiers by
  stripping the full internal dataset-id prefix, so source identifiers
  that themselves contain a hyphen keep the correct site suffix.

  

## VerteTIME 1.0.0

Initial public release.

### New features

- Ingestion pipeline:
  [`vt_ingest_dataset()`](https://robustecologies.github.io/VerteTIME/reference/vt_ingest_dataset.md),
  [`vt_ingest_all()`](https://robustecologies.github.io/VerteTIME/reference/vt_ingest_all.md),
  [`vt_validate()`](https://robustecologies.github.io/VerteTIME/reference/vt_validate.md).
- User-facing dataset scaffolder:
  `vt_register_dataset(id, data_raw, dry_run)` copies the package YAML
  template into a user-supplied private tree, runs the same
  ingestion-and-validation pipeline that produced the shipped
  compilation against the user’s tree, and returns the freshly built
  `vt_dataset` together with the validation tibble. The function
  operates exclusively against the user-supplied `data_raw` argument and
  never touches the maintainer-side build tree.
- Reading and reshaping:
  [`vt_read()`](https://robustecologies.github.io/VerteTIME/reference/vt_read.md),
  [`vt_long()`](https://robustecologies.github.io/VerteTIME/reference/vt_long.md),
  [`vt_wide()`](https://robustecologies.github.io/VerteTIME/reference/vt_wide.md),
  [`vt_compilation()`](https://robustecologies.github.io/VerteTIME/reference/vt_compilation.md).
  `vt_read("VT_NNN")` reconstructs a `vt_dataset` from the shipped
  `vertetime` compilation when no maintainer-side ingestion tree is
  supplied, so end users can interrogate a single dataset slice with one
  call.
- Alpha-diversity metrics:
  [`vt_alpha_diversity()`](https://robustecologies.github.io/VerteTIME/reference/vt_alpha_diversity.md)
  (richness, Shannon, Simpson, Hill numbers q=0,1,2, Chao1) and
  [`vt_evenness()`](https://robustecologies.github.io/VerteTIME/reference/vt_evenness.md)
  (Pielou, Simpson).
- Beta-diversity and turnover:
  [`vt_beta_diversity()`](https://robustecologies.github.io/VerteTIME/reference/vt_beta_diversity.md)
  (Bray-Curtis, Jaccard, Hellinger, Whittaker multiplicative),
  [`vt_temporal_turnover()`](https://robustecologies.github.io/VerteTIME/reference/vt_temporal_turnover.md),
  [`vt_beta_partition()`](https://robustecologies.github.io/VerteTIME/reference/vt_beta_partition.md)
  (Baselga 2010 turnover/nestedness components).
- Composition:
  [`vt_dominance_curve()`](https://robustecologies.github.io/VerteTIME/reference/vt_dominance_curve.md),
  [`vt_rank_abundance()`](https://robustecologies.github.io/VerteTIME/reference/vt_rank_abundance.md),
  [`vt_rarefaction()`](https://robustecologies.github.io/VerteTIME/reference/vt_rarefaction.md),
  [`vt_nestedness()`](https://robustecologies.github.io/VerteTIME/reference/vt_nestedness.md)
  (NODF).
- Aggregators and filters:
  [`vt_community_summary()`](https://robustecologies.github.io/VerteTIME/reference/vt_community_summary.md),
  [`vt_compilation_summary()`](https://robustecologies.github.io/VerteTIME/reference/vt_compilation_summary.md),
  [`vt_filter()`](https://robustecologies.github.io/VerteTIME/reference/vt_filter.md).
- Multi-format export:
  [`vt_export()`](https://robustecologies.github.io/VerteTIME/reference/vt_export.md)
  writing CSV, Apache Parquet, SQLite single-file, and a Frictionless
  Data Package;
  [`vt_publish()`](https://robustecologies.github.io/VerteTIME/reference/vt_publish.md)
  orchestrator producing the full public-release tree into a
  user-supplied directory.
- Visualisation:
  [`vt_plot_whittaker()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_whittaker.md),
  [`vt_plot_rac()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_rac.md),
  [`vt_plot_alpha_timeline()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_alpha_timeline.md),
  [`vt_plot_beta_heatmap()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_beta_heatmap.md),
  [`vt_plot_world_map()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_world_map.md)
  (Robinson projection), and
  [`vt_plot_database_timeline()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_database_timeline.md)
  (positional comparator across vertebrate time-series compilations).
  The two compilation-level plot methods accept a `type =` argument that
  names the dispatched plotter; the `kind =` alias is retained
  indefinitely as a back-compatible synonym.
- S3 classes with bidirectional `print` / `summary` / `plot`:
  `vt_dataset`, `vt_compilation`, `vt_diversity`, `vt_turnover`,
  `vt_provenance`.
  [`summary.vt_turnover()`](https://robustecologies.github.io/VerteTIME/reference/vt_turnover.md)
  dispatches on the columns of the object (Baselga partition components
  or `metric`/`value` long form) and returns a tidy tibble keyed by
  metric with `n`, `median`, `q25`, `q75`, `mean`, `sd`.
- Installed-package self-test: `vt_check(verbose, plots)` loads
  `vt_demo` into a private environment, exercises one representative
  call per analytical family, runs every S3 `print` / `summary` / `plot`
  lifecycle method, exercises
  [`vt_publish()`](https://robustecologies.github.io/VerteTIME/reference/vt_publish.md)
  and
  [`vt_export()`](https://robustecologies.github.io/VerteTIME/reference/vt_export.md)
  against [`tempfile()`](https://rdrr.io/r/base/tempfile.html) paths,
  and returns a tibble whose `status` column is `"ok"` for every row
  when the package is healthy.

### Vignettes

- `vt-introduction`: a five-minute tour using a synthetic mid-sized
  compilation.
- `vt-community-metrics`: alpha/beta/turnover/nestedness deep dive.
- `vt-ingestion-workflow`: how to register a private dataset against the
  shipped compilation through
  [`vt_register_dataset()`](https://robustecologies.github.io/VerteTIME/reference/vt_register_dataset.md).
- `API`: technical architecture document covering the five-table
  relational schema, the ingestion pipeline, the diversity catalogue,
  the visualisation engine and the publication subsystem.
