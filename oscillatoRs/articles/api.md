# API reference and architecture

The package surfaces 13 exported functions across four functional
layers: data access, index registry, visualization, and ggplot2
integration. This vignette is a navigational map. Section one states the
design philosophy. Section two diagrams each subsystem. Section three
catalogues the function families with cross-links to the reference
index.

  

## Design philosophy

oscillatoRs is a thin, dependency-light client for authoritative
climate-index time series. Three principles govern the architecture:

The first is **provider transparency**. Every index carries explicit
provenance: provider name (NOAA PSL, NOAA CPC, NCEP/NCAR), retrieval
URL, parser dispatch and last-update timestamp. The registry in
`R/index_registry.R` is the single source of truth; downstream code
never hard-codes URLs.

The second is **bundled defaults with optional refresh**. The package
ships pre-bundled datasets (`climate_monthly`, `climate_monthly_wide`,
`climate_metadata`) so a vignette renders or a reproducible analysis
runs without internet.
[`update_indices()`](https://robustecologies.github.io/oscillatoRs/reference/update_indices.md)
opts in to live retrieval and merges new observations into the bundled
tibbles.

The third is **ggplot2-native visualization with minimal opinion**. Plot
functions return `ggplot` objects with a tuned
[`theme_climate()`](https://robustecologies.github.io/oscillatoRs/reference/theme_climate.md)
and anomaly-aware
[`scale_color_anomaly()`](https://robustecologies.github.io/oscillatoRs/reference/scale_color_anomaly.md)
/
[`scale_fill_anomaly()`](https://robustecologies.github.io/oscillatoRs/reference/scale_fill_anomaly.md)
scales. Users compose, facet and annotate as in any other ggplot2
pipeline.

  

## Architecture overview

  

## Data-access pipeline

[`get_index()`](https://robustecologies.github.io/oscillatoRs/reference/get_index.md)
is the canonical accessor for a single named index;
[`get_indices()`](https://robustecologies.github.io/oscillatoRs/reference/get_indices.md)
is its vectorised counterpart. Both consult the registry, dispatch the
appropriate parser, hit the provider URL when called with
`download = TRUE`, and otherwise fall back to the bundled tibbles.
[`update_indices()`](https://robustecologies.github.io/oscillatoRs/reference/update_indices.md)
refreshes the bundled cache in place.

  

## Index registry

The registry organises 43 indices into nine semantic categories. The
diagram below names the categories and their canonical members.

  

## Visualization layer

The visualization stack is opinionated about a single thing: anomaly
plots should shade positive and negative excursions on either side of
zero, with a neutral midpoint. Everything else is a thin ggplot2
wrapper.

  

## Cache and update flow

[`update_indices()`](https://robustecologies.github.io/oscillatoRs/reference/update_indices.md)
walks the registry, calls each provider, parses, and merges new rows
into the bundled tibbles. Failed providers do not abort the refresh;
they are logged and the bundle keeps its prior state for the affected
indices.

  

## Function-family catalogue

The four reference families are mirrored in `_pkgdown.yml` under
`reference:`. Each family is exhaustively listed there and each function
carries a complete roxygen page.

The **data access** family (`get_index`, `get_indices`,
`update_indices`) is the only entry point for fetching observations.
Internally it dispatches to one of three parsers (PSL flat-file, MEI
multi-column, daily PSL) selected by the registry.

The **index registry** family (`list_indices`, `get_index_info`) exposes
the catalogue.
[`list_indices()`](https://robustecologies.github.io/oscillatoRs/reference/list_indices.md)
returns the full registry as a tibble;
[`get_index_info()`](https://robustecologies.github.io/oscillatoRs/reference/get_index_info.md)
returns the metadata block (provider, URL, parser, frequency, units) for
a single index.

The **visualization** family (`plot_index`, `plot_indices`,
`plot_correlation`, `plot_availability`) returns `ggplot` objects with
[`theme_climate()`](https://robustecologies.github.io/oscillatoRs/reference/theme_climate.md)
and the anomaly scales pre-applied. All four carry a mandatory subtitle
and caption per the package convention.

The **ggplot2 integration** family (`scale_color_anomaly`,
`scale_fill_anomaly`) is the diverging palette tuned for
centered-at-zero anomalies. Drop-in for any
[`ggplot()`](https://ggplot2.tidyverse.org/reference/ggplot.html) that
maps `value` (or its z-scored variant) to `colour` or `fill`.

  

## Datasets

`climate_monthly` is the long-format tibble (date, index, value)
covering every monthly index in the registry from the start of the
available record to the most recent bundle date. `climate_monthly_wide`
is the same data pivoted to wide format (one column per index).
`climate_metadata` carries the registry metadata as a tibble for joining
with analysis output.

  

## Reproducibility

All bundled tibbles are tagged with a `last_updated` attribute
reflecting the date of the bundle. Reproducible analyses pin to the
package version and use the bundled data; live analyses call
[`update_indices()`](https://robustecologies.github.io/oscillatoRs/reference/update_indices.md)
and record the new attribute.
