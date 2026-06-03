# Monthly climate oscillation indices

A dataset containing 38 monthly climate oscillation indices from NOAA
Physical Sciences Laboratory (PSL) and Climate Prediction Center (CPC).
Data are in long format with one row per index-month observation.

## Usage

``` r
climate_monthly
```

## Format

A tibble with the following columns:

- date:

  Date. First day of the month.

- year:

  Integer. Calendar year.

- month:

  Integer. Month (1-12).

- value:

  Numeric. Index value, typically standardized anomalies.

- index:

  Character. Index identifier (e.g., "NAO", "PDO", "Nino34").

- description:

  Character. Short description of the index.

- category:

  Character. Index category (e.g., "ENSO", "Teleconnection").

- source_url:

  Character. Data source URL.

## Source

NOAA Physical Sciences Laboratory: <https://psl.noaa.gov/>

NOAA Climate Prediction Center: <https://www.cpc.ncep.noaa.gov/>

## Details

The dataset includes indices from the following categories:

- **ENSO** (8 indices): SOI, Nino12, Nino3, Nino34, Nino4, ONI, BEST,
  TNI

- **Teleconnection** (8 indices): NAO, NAO_Jones, PNA, WP, EPNP, EAWR,
  NP, NOI

- **Pacific multidecadal** (2 indices): PDO, IPO_TPI

- **Atlantic** (8 indices): AMO_unsmoothed, AMO_smoothed, TNA, TSA,
  WHWP, Atlantic_Tripole, NTA, CAR

- **Polar** (2 indices): AO, AAO

- **Indian Ocean** (1 index): DMI

- **Atmospheric** (1 index): QBO

- **Pacific SST** (2 indices): Pacific_Warmpool, Tropical_Pacific_EOF

- **Precipitation** (2 indices): Indian_Monsoon, Sahel_Rainfall

- **Solar** (1 index): Solar_Flux

- **Temperature** (1 index): Global_Temp

Data coverage varies by index, with some extending back to the 1800s
(e.g., NAO_Jones from 1821, SOI from 1866) and others starting in the
satellite era (e.g., AAO from 1979).

## See also

[climate_monthly_wide](https://robustecologies.github.io/oscillatoRs/reference/climate_monthly_wide.md)
for wide format,
[climate_metadata](https://robustecologies.github.io/oscillatoRs/reference/climate_metadata.md)
for index descriptions,
[`get_index()`](https://robustecologies.github.io/oscillatoRs/reference/get_index.md)
to download updated data

## Examples

``` r
if (FALSE) { # \dontrun{
data(climate_monthly)
head(climate_monthly)

# Filter to ENSO indices
enso <- climate_monthly[climate_monthly$category == "ENSO", ]

# Get NAO time series
nao <- climate_monthly[climate_monthly$index == "NAO", ]
} # }
```
