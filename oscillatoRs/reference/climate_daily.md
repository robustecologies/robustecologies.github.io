# Daily climate oscillation indices

Daily values for atmospheric teleconnection indices from NOAA Climate
Prediction Center (CPC). Useful for high-frequency climate monitoring
and synoptic-scale analysis.

## Usage

``` r
climate_daily
```

## Format

A tibble with the following columns:

- date:

  Date. Observation date.

- year:

  Integer. Calendar year.

- month:

  Integer. Month (1-12).

- day:

  Integer. Day of month (1-31).

- value:

  Numeric. Index value (standardized).

- index:

  Character. Index identifier.

- description:

  Character. Short description.

- category:

  Character. Index category.

- source_url:

  Character. Data source URL.

## Source

NOAA Climate Prediction Center: <https://www.cpc.ncep.noaa.gov/>

## Details

Four daily indices are available:

- **NAO_daily**: North Atlantic Oscillation (1950-present)

- **AO_daily**: Arctic Oscillation (1950-present)

- **PNA_daily**: Pacific North American (1950-present)

- **AAO_daily**: Antarctic Oscillation (1979-present)

Daily indices are computed from EOF analysis of geopotential height
fields and represent projections onto the leading teleconnection
patterns.

## See also

[climate_monthly](https://robustecologies.github.io/oscillatoRs/reference/climate_monthly.md)
for monthly aggregates,
[`get_index()`](https://robustecologies.github.io/oscillatoRs/reference/get_index.md)
to download updated data

## Examples

``` r
if (FALSE) { # \dontrun{
data(climate_daily)

# Recent 30 days of NAO
nao_recent <- climate_daily[
  climate_daily$index == "NAO_daily" &
  climate_daily$date > Sys.Date() - 30,
]
} # }
```
