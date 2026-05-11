# Changelog

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
