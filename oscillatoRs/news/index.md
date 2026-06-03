# Changelog

## oscillatoRs 0.1.1

### Changes

- [`plot_indices()`](https://robustecologies.github.io/oscillatoRs/reference/plot_indices.md)
  now carries the mandatory subtitle and caption that the other `plot_*`
  methods already produced. No change to returned data; visualization is
  now uniform across the visualization layer.

  

## oscillatoRs 0.1.0

Initial release of the oscillatoRs package. \## Features

- **43 climate oscillation indices** from NOAA PSL and CPC:

  - ENSO indices: SOI, Nino12, Nino3, Nino34, Nino4, ONI, BEST, TNI, MEI
    v2
  - Teleconnections: NAO (CPC and Jones), PNA, WP, EP/NP, EA/WR, NP, NOI
  - Pacific multidecadal: PDO, IPO (TPI)
  - Atlantic: AMO (smoothed and unsmoothed), TNA, TSA, WHWP, Atlantic
    Tripole, NTA, CAR
  - Polar: AO, AAO
  - Indian Ocean: DMI
  - Atmospheric: QBO
  - Pacific SST: Pacific Warmpool, Tropical Pacific EOF
  - Precipitation: Indian Monsoon, Sahel Rainfall
  - Other: Solar Flux, Global Temperature

- **Daily indices** (4): NAO, AO, PNA, AAO

- **Pre-bundled datasets**:

  - `climate_monthly`: Long-format monthly data
  - `climate_monthly_wide`: Wide-format monthly data
  - `climate_daily`: Daily teleconnection indices
  - `climate_mei`: MEI v2 bimonthly data
  - `climate_metadata`: Index descriptions and sources
  - `climate_summary`: Summary statistics

- **Download functions**:

  - [`get_index()`](https://robustecologies.github.io/oscillatoRs/reference/get_index.md):
    Download single index
  - [`get_indices()`](https://robustecologies.github.io/oscillatoRs/reference/get_indices.md):
    Download multiple indices
  - [`list_indices()`](https://robustecologies.github.io/oscillatoRs/reference/list_indices.md):
    List available indices
  - [`get_index_info()`](https://robustecologies.github.io/oscillatoRs/reference/get_index_info.md):
    Get index metadata
  - [`update_indices()`](https://robustecologies.github.io/oscillatoRs/reference/update_indices.md):
    Update bundled data

- **Visualization functions**:

  - [`plot_index()`](https://robustecologies.github.io/oscillatoRs/reference/plot_index.md):
    Time series plot of single index
  - [`plot_indices()`](https://robustecologies.github.io/oscillatoRs/reference/plot_indices.md):
    Faceted comparison of multiple indices
  - [`plot_correlation()`](https://robustecologies.github.io/oscillatoRs/reference/plot_correlation.md):
    Correlation heatmap
  - [`plot_availability()`](https://robustecologies.github.io/oscillatoRs/reference/plot_availability.md):
    Data coverage timeline
  - [`theme_climate()`](https://robustecologies.github.io/oscillatoRs/reference/theme_climate.md):
    Publication-ready ggplot2 theme
  - [`scale_fill_anomaly()`](https://robustecologies.github.io/oscillatoRs/reference/scale_fill_anomaly.md),
    [`scale_color_anomaly()`](https://robustecologies.github.io/oscillatoRs/reference/scale_color_anomaly.md):
    Diverging color scales
  - [`category_colors()`](https://robustecologies.github.io/oscillatoRs/reference/category_colors.md):
    Category color palette

- **Vignettes**:

  - User guide with examples
  - Comprehensive theory vignette covering physics of climate
    oscillations
