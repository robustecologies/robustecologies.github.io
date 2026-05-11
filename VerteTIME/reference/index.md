# Package index

## Ingestion, validation and registration

Build a compilation from a maintainer-private ingestion tree, validate
the relational schema, scaffold a user-side YAML sidecar, and run the
installed-package self-test.

- [`vt_ingest_dataset()`](https://robustecologies.github.io/VerteTIME/reference/vt_ingest_dataset.md)
  : Ingest a single dataset folder into the canonical relational schema

- [`vt_ingest_all()`](https://robustecologies.github.io/VerteTIME/reference/vt_ingest_all.md)
  : Ingest the full compilation from a maintainer-side ingestion tree

- [`vt_register_dataset()`](https://robustecologies.github.io/VerteTIME/reference/vt_register_dataset.md)
  : Register a new dataset in a private VerteTIME extension tree

- [`vt_validate()`](https://robustecologies.github.io/VerteTIME/reference/vt_validate.md)
  :

  Validate a `vt_dataset` or `vt_compilation`

- [`vt_check()`](https://robustecologies.github.io/VerteTIME/reference/vt_check.md)
  : One-call smoke test of the installed VerteTIME package

- [`enrich_taxonomy()`](https://robustecologies.github.io/VerteTIME/reference/enrich_taxonomy.md)
  : Enrich the species table with class/order/family from the taxonomy
  cache

## Reading and reshaping

Low-level CSV IO, long-wide pivot helpers and S3 constructors for
datasets, compilations and provenance objects.

- [`vt_read_csv()`](https://robustecologies.github.io/VerteTIME/reference/vt_read_csv.md)
  : Read a wide-format community CSV

- [`vt_read()`](https://robustecologies.github.io/VerteTIME/reference/vt_read.md)
  : Read a VerteTIME dataset

- [`vt_long()`](https://robustecologies.github.io/VerteTIME/reference/vt_long.md)
  : Pivot the wide observation block of a VerteTIME dataset into long
  form

- [`vt_wide()`](https://robustecologies.github.io/VerteTIME/reference/vt_wide.md)
  : Pivot a long-format observation table back to wide

- [`classify_columns()`](https://robustecologies.github.io/VerteTIME/reference/classify_columns.md)
  : Classify column names into species and covariates

- [`vt_compilation()`](https://robustecologies.github.io/VerteTIME/reference/vt_compilation.md)
  [`print(`*`<vt_compilation>`*`)`](https://robustecologies.github.io/VerteTIME/reference/vt_compilation.md)
  [`summary(`*`<vt_compilation>`*`)`](https://robustecologies.github.io/VerteTIME/reference/vt_compilation.md)
  [`plot(`*`<vt_compilation>`*`)`](https://robustecologies.github.io/VerteTIME/reference/vt_compilation.md)
  :

  Construct a `vt_compilation` object

- [`vt_dataset()`](https://robustecologies.github.io/VerteTIME/reference/vt_dataset.md)
  [`print(`*`<vt_dataset>`*`)`](https://robustecologies.github.io/VerteTIME/reference/vt_dataset.md)
  [`summary(`*`<vt_dataset>`*`)`](https://robustecologies.github.io/VerteTIME/reference/vt_dataset.md)
  [`plot(`*`<vt_dataset>`*`)`](https://robustecologies.github.io/VerteTIME/reference/vt_dataset.md)
  :

  Construct a `vt_dataset` object

- [`vt_provenance()`](https://robustecologies.github.io/VerteTIME/reference/vt_provenance.md)
  [`print(`*`<vt_provenance>`*`)`](https://robustecologies.github.io/VerteTIME/reference/vt_provenance.md)
  :

  Construct a `vt_provenance` object

## Diversity indices

Per-(dataset, site, year) alpha-diversity, Hill numbers and evenness
measures.

- [`vt_alpha_diversity()`](https://robustecologies.github.io/VerteTIME/reference/vt_alpha_diversity.md)
  : Compute alpha-diversity indices per site and year

- [`vt_evenness()`](https://robustecologies.github.io/VerteTIME/reference/vt_evenness.md)
  : Compute evenness indices per site and year

- [`vt_diversity()`](https://robustecologies.github.io/VerteTIME/reference/vt_diversity.md)
  [`print(`*`<vt_diversity>`*`)`](https://robustecologies.github.io/VerteTIME/reference/vt_diversity.md)
  [`summary(`*`<vt_diversity>`*`)`](https://robustecologies.github.io/VerteTIME/reference/vt_diversity.md)
  [`plot(`*`<vt_diversity>`*`)`](https://robustecologies.github.io/VerteTIME/reference/vt_diversity.md)
  :

  Construct a `vt_diversity` object

- [`vt_series_summary()`](https://robustecologies.github.io/VerteTIME/reference/vt_series_summary.md)
  : Per-series temporal-coverage and variability metrics for a long
  observation table

## Beta diversity and temporal turnover

Pairwise dissimilarity across sites and years, with the Baselga
turnover-nestedness decomposition.

- [`vt_beta_diversity()`](https://robustecologies.github.io/VerteTIME/reference/vt_beta_diversity.md)
  : Compute beta-diversity matrices across sites within each dataset

- [`vt_temporal_turnover()`](https://robustecologies.github.io/VerteTIME/reference/vt_temporal_turnover.md)
  : Compute temporal turnover (year-pair beta) per site

- [`vt_beta_partition()`](https://robustecologies.github.io/VerteTIME/reference/vt_beta_partition.md)
  : Partition Sorensen dissimilarity into turnover and nestedness
  components

- [`vt_turnover()`](https://robustecologies.github.io/VerteTIME/reference/vt_turnover.md)
  [`print(`*`<vt_turnover>`*`)`](https://robustecologies.github.io/VerteTIME/reference/vt_turnover.md)
  [`summary(`*`<vt_turnover>`*`)`](https://robustecologies.github.io/VerteTIME/reference/vt_turnover.md)
  [`plot(`*`<vt_turnover>`*`)`](https://robustecologies.github.io/VerteTIME/reference/vt_turnover.md)
  :

  Construct a `vt_turnover` object

## Composition

Rank-abundance, dominance, Hurlbert-Sanders rarefaction and NODF
nestedness diagnostics.

- [`vt_dominance_curve()`](https://robustecologies.github.io/VerteTIME/reference/vt_dominance_curve.md)
  : Dominance curves per site and year
- [`vt_rank_abundance()`](https://robustecologies.github.io/VerteTIME/reference/vt_rank_abundance.md)
  : Rank-abundance per site and year
- [`vt_rarefaction()`](https://robustecologies.github.io/VerteTIME/reference/vt_rarefaction.md)
  : Coverage-based rarefaction per site and year
- [`vt_nestedness()`](https://robustecologies.github.io/VerteTIME/reference/vt_nestedness.md)
  : Nestedness based on overlap and decreasing fill (NODF)

## Aggregators and filtering

Per-site-year summaries, compilation-level summaries and tidy-eval
filtering with FK-cascade integrity.

- [`vt_community_summary()`](https://robustecologies.github.io/VerteTIME/reference/vt_community_summary.md)
  : Per-site-year community summary

- [`vt_compilation_summary()`](https://robustecologies.github.io/VerteTIME/reference/vt_compilation_summary.md)
  : Compilation-level summary table

- [`vt_filter()`](https://robustecologies.github.io/VerteTIME/reference/vt_filter.md)
  :

  Filter a compilation by a tidy-eval predicate against the `datasets`
  table

## Export and publication

Single-format export, full publish tree and the LPD-style wide-master
series table.

- [`vt_export()`](https://robustecologies.github.io/VerteTIME/reference/vt_export.md)
  : Export a compilation to a single file format

- [`vt_publish()`](https://robustecologies.github.io/VerteTIME/reference/vt_publish.md)
  :

  Publish a compilation as a complete `web-export/` tree

- [`vt_wide_master()`](https://robustecologies.github.io/VerteTIME/reference/vt_wide_master.md)
  : Build the wide-format master series table for the compilation

## Visualisation theme and palettes

Project-wide ggplot theme; six-level focus palette and eleven-level
GBIF-class palette.

- [`vt_theme()`](https://robustecologies.github.io/VerteTIME/reference/vt_theme.md)
  : VerteTIME ggplot2 theme
- [`vt_palette()`](https://robustecologies.github.io/VerteTIME/reference/vt_palette.md)
  : Six-level taxonomic-focus palette
- [`vt_class_palette()`](https://robustecologies.github.io/VerteTIME/reference/vt_class_palette.md)
  : Class-level colour palette aligned with the GBIF backbone classes

## Plots - taxonomic structure

Class, focus, family and species composition of the compilation.

- [`vt_plot_class_bars()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_class_bars.md)
  : Two-panel bar chart of series and species per GBIF class
- [`vt_plot_class_treemap()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_class_treemap.md)
  : Hierarchical Class -\> Order -\> Family treemap
- [`vt_plot_focus_bars()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_focus_bars.md)
  : Bar charts of dataset count and species count by taxonomic focus
- [`vt_plot_top_species()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_top_species.md)
  : Top-N most monitored species across the compilation
- [`vt_plot_top_genera()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_top_genera.md)
  : Top-N genera by number of datasets they appear in

## Plots - geographic coverage

World map, hexbin density, latitudinal histogram, country and realm
breakdowns.

- [`vt_plot_world_map()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_world_map.md)
  : World map of compilation sites
- [`vt_plot_world_hex()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_world_hex.md)
  : Hexbin density of compilation sites on a Robinson world map
- [`vt_plot_world_richness()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_world_richness.md)
  : World richness heatmap on a coarse Robinson grid
- [`vt_plot_lat_hist()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_lat_hist.md)
  : Latitudinal histogram of compilation sites, faceted by realm
- [`vt_plot_country_bar()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_country_bar.md)
  : Top-N countries by number of datasets
- [`vt_plot_realm_bar()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_realm_bar.md)
  : Datasets per realm
- [`vt_plot_focus_realm_heatmap()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_focus_realm_heatmap.md)
  : Count heatmap of taxonomic focus by realm
- [`vt_plot_alluvial()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_alluvial.md)
  : Alluvial flow from taxonomic focus through realm to continent

## Plots - temporal coverage and missingness

Lasagna, active area, onset, year ridges, decadal Sankey and missingness
diagnostics.

- [`vt_plot_lasagna()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_lasagna.md)
  : Lasagna plot of dataset-by-year presence
- [`vt_plot_active_area()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_active_area.md)
  : Stacked-area count of datasets active per calendar year
- [`vt_plot_onset_cumulative()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_onset_cumulative.md)
  : Cumulative onset curve
- [`vt_plot_year_ridges()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_year_ridges.md)
  : Ridge plot of first-observed year by taxonomic focus
- [`vt_plot_temporal_sankey()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_temporal_sankey.md)
  : Temporal sankey of taxonomic focus across decades
- [`vt_plot_miss_by_year()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_miss_by_year.md)
  : Missingness rate per year by realm
- [`vt_plot_miss_class_heatmap()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_miss_class_heatmap.md)
  : Year-by-class missingness heatmap
- [`vt_plot_span_hist()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_span_hist.md)
  : Histograms of dataset span and observation count
- [`vt_plot_focus_span_box()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_focus_span_box.md)
  : Boxplot of dataset span by taxonomic focus
- [`vt_plot_span_coverage_hex()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_span_coverage_hex.md)
  : Span versus coverage scatter (hexbin density)
- [`vt_plot_mean_cv()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_mean_cv.md)
  : Mean abundance versus coefficient of variation per species

## Plots - community pattern

Whittaker, rank-abundance, alpha timeline and histogram, beta heatmap,
NODF distribution, Hill profile and per-community spotlights.

- [`vt_plot_whittaker()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_whittaker.md)
  : Whittaker rank-abundance plot
- [`vt_plot_rac()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_rac.md)
  : Rank-abundance curves with year facets
- [`vt_plot_alpha_timeline()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_alpha_timeline.md)
  : Alpha-diversity timeline
- [`vt_plot_alpha_hist()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_alpha_hist.md)
  : Distribution of an alpha-diversity index across the compilation
- [`vt_plot_beta_heatmap()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_beta_heatmap.md)
  : Beta-diversity / temporal-turnover heatmap
- [`vt_plot_nodf_hist()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_nodf_hist.md)
  : NODF nestedness across the compilation
- [`vt_plot_metric_ridges()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_metric_ridges.md)
  : Multi-metric ridge plot of community diagnostics
- [`vt_plot_hill_profile()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_hill_profile.md)
  : Hill-number profile per dataset across orders 0 to 3
- [`vt_plot_community_series()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_community_series.md)
  : Plot the abundance trajectory of every species in a community
- [`vt_plot_community_spotlight()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_community_spotlight.md)
  : Featured-community panel: time series + world location for one
  dataset

## Plots - positional comparator

Vertebrate community time-series compilations placed alongside VerteTIME
for orientation, never as a provenance chain.

- [`vt_plot_database_timeline()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_database_timeline.md)
  : Positional timeline of biodiversity time-series compilations
- [`vt_compilations()`](https://robustecologies.github.io/VerteTIME/reference/vt_compilations.md)
  : Built-in table of biodiversity time-series compilations cited by
  VerteTIME

## Shipped data objects

Canonical compilation and single-dataset slice loaded by data().

- [`vertetime`](https://robustecologies.github.io/VerteTIME/reference/vertetime.md)
  : VerteTIME canonical compilation
- [`vt_demo`](https://robustecologies.github.io/VerteTIME/reference/vt_demo.md)
  : VerteTIME single-dataset demo compilation

## S3 methods

Print, summary and plot methods attached to the exported S3 classes.
Users normally do not call these directly; they are listed here for
completeness and cross-reference.

- [`vt_compilation()`](https://robustecologies.github.io/VerteTIME/reference/vt_compilation.md)
  [`print(`*`<vt_compilation>`*`)`](https://robustecologies.github.io/VerteTIME/reference/vt_compilation.md)
  [`summary(`*`<vt_compilation>`*`)`](https://robustecologies.github.io/VerteTIME/reference/vt_compilation.md)
  [`plot(`*`<vt_compilation>`*`)`](https://robustecologies.github.io/VerteTIME/reference/vt_compilation.md)
  :

  Construct a `vt_compilation` object

- [`vt_dataset()`](https://robustecologies.github.io/VerteTIME/reference/vt_dataset.md)
  [`print(`*`<vt_dataset>`*`)`](https://robustecologies.github.io/VerteTIME/reference/vt_dataset.md)
  [`summary(`*`<vt_dataset>`*`)`](https://robustecologies.github.io/VerteTIME/reference/vt_dataset.md)
  [`plot(`*`<vt_dataset>`*`)`](https://robustecologies.github.io/VerteTIME/reference/vt_dataset.md)
  :

  Construct a `vt_dataset` object

- [`vt_diversity()`](https://robustecologies.github.io/VerteTIME/reference/vt_diversity.md)
  [`print(`*`<vt_diversity>`*`)`](https://robustecologies.github.io/VerteTIME/reference/vt_diversity.md)
  [`summary(`*`<vt_diversity>`*`)`](https://robustecologies.github.io/VerteTIME/reference/vt_diversity.md)
  [`plot(`*`<vt_diversity>`*`)`](https://robustecologies.github.io/VerteTIME/reference/vt_diversity.md)
  :

  Construct a `vt_diversity` object

- [`vt_provenance()`](https://robustecologies.github.io/VerteTIME/reference/vt_provenance.md)
  [`print(`*`<vt_provenance>`*`)`](https://robustecologies.github.io/VerteTIME/reference/vt_provenance.md)
  :

  Construct a `vt_provenance` object

- [`vt_turnover()`](https://robustecologies.github.io/VerteTIME/reference/vt_turnover.md)
  [`print(`*`<vt_turnover>`*`)`](https://robustecologies.github.io/VerteTIME/reference/vt_turnover.md)
  [`summary(`*`<vt_turnover>`*`)`](https://robustecologies.github.io/VerteTIME/reference/vt_turnover.md)
  [`plot(`*`<vt_turnover>`*`)`](https://robustecologies.github.io/VerteTIME/reference/vt_turnover.md)
  :

  Construct a `vt_turnover` object
