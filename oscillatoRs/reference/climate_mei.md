# Multivariate ENSO index (MEI v2)

The Multivariate ENSO Index version 2 (MEI.v2) combines multiple
atmospheric and oceanic variables into a single index that captures ENSO
variability more completely than single-variable indices.

## Usage

``` r
climate_mei
```

## Format

A tibble with the following columns:

- date:

  Date. Center of bimonthly period (15th of second month).

- year:

  Integer. Calendar year.

- month:

  Integer. Second month of bimonthly period.

- bimonth:

  Character. Bimonthly period label (e.g., "DJ", "JF").

- value:

  Numeric. MEI value (standardized).

- index:

  Character. "MEI_v2".

- description:

  Character. Index description.

- category:

  Character. "ENSO".

- source_url:

  Character. Data source URL.

## Source

NOAA PSL: <https://psl.noaa.gov/enso/mei/>

## Details

MEI.v2 combines six variables over the tropical Pacific:

1.  Sea level pressure

2.  Sea surface temperature

3.  Surface zonal wind

4.  Surface meridional wind

5.  Outgoing longwave radiation

The bimonthly periods are sliding two-month windows: DJ (Dec-Jan), JF
(Jan-Feb), FM (Feb-Mar), MA (Mar-Apr), AM (Apr-May), MJ (May-Jun), JJ
(Jun-Jul), JA (Jul-Aug), AS (Aug-Sep), SO (Sep-Oct), ON (Oct-Nov), ND
(Nov-Dec).

Positive values indicate El Nino conditions; negative indicate La Nina.

## References

Wolter, K. and M.S. Timlin (2011), El Nino/Southern Oscillation
behaviour since 1871 as diagnosed in an extended multivariate ENSO index
(MEI.ext). Int. J. Climatol., 31: 1074-1087.
[doi:10.1002/joc.2336](https://doi.org/10.1002/joc.2336)

## Examples

``` r
if (FALSE) { # \dontrun{
data(climate_mei)
head(climate_mei)
} # }
```
