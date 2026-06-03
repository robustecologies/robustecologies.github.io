# Package index

## Data access

Download and refresh climate-index time series from authoritative
providers.

- [`get_index()`](https://robustecologies.github.io/oscillatoRs/reference/get_index.md)
  : Download a single climate index
- [`get_indices()`](https://robustecologies.github.io/oscillatoRs/reference/get_indices.md)
  : Download multiple climate indices
- [`update_indices()`](https://robustecologies.github.io/oscillatoRs/reference/update_indices.md)
  : Update bundled package data

## Index registry

Discover available indices and inspect their metadata.

- [`list_indices()`](https://robustecologies.github.io/oscillatoRs/reference/list_indices.md)
  : List available climate indices
- [`get_index_info()`](https://robustecologies.github.io/oscillatoRs/reference/get_index_info.md)
  : Get detailed information about a climate index

## Visualization

ggplot2 helpers tuned to anomaly time series, correlation matrices and
availability charts.

- [`plot_index()`](https://robustecologies.github.io/oscillatoRs/reference/plot_index.md)
  : Plot a single climate index
- [`plot_indices()`](https://robustecologies.github.io/oscillatoRs/reference/plot_indices.md)
  : Plot multiple climate indices
- [`plot_correlation()`](https://robustecologies.github.io/oscillatoRs/reference/plot_correlation.md)
  : Plot correlation matrix of climate indices
- [`plot_availability()`](https://robustecologies.github.io/oscillatoRs/reference/plot_availability.md)
  : Plot data availability timeline
- [`theme_climate()`](https://robustecologies.github.io/oscillatoRs/reference/theme_climate.md)
  : Climate visualization theme
- [`category_colors()`](https://robustecologies.github.io/oscillatoRs/reference/category_colors.md)
  : Category color palette

## ggplot2 integration

Drop-in colour and fill scales aligned with positive/negative anomaly
conventions.

- [`scale_color_anomaly()`](https://robustecologies.github.io/oscillatoRs/reference/scale_color_anomaly.md)
  : Anomaly color scale
- [`scale_fill_anomaly()`](https://robustecologies.github.io/oscillatoRs/reference/scale_fill_anomaly.md)
  : Anomaly fill scale

## Bundled datasets

Pre-bundled tibbles. Reproducible analyses pin to the package version;
live analyses call update_indices().

- [`climate_monthly`](https://robustecologies.github.io/oscillatoRs/reference/climate_monthly.md)
  : Monthly climate oscillation indices
- [`climate_monthly_wide`](https://robustecologies.github.io/oscillatoRs/reference/climate_monthly_wide.md)
  : Monthly climate indices in wide format
- [`climate_daily`](https://robustecologies.github.io/oscillatoRs/reference/climate_daily.md)
  : Daily climate oscillation indices
- [`climate_mei`](https://robustecologies.github.io/oscillatoRs/reference/climate_mei.md)
  : Multivariate ENSO index (MEI v2)
- [`climate_metadata`](https://robustecologies.github.io/oscillatoRs/reference/climate_metadata.md)
  : Climate index metadata
- [`climate_summary`](https://robustecologies.github.io/oscillatoRs/reference/climate_summary.md)
  : Summary statistics for monthly climate indices
