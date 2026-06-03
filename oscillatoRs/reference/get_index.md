# Download a single climate index

Downloads and parses a climate oscillation index from its source URL.
Supports monthly, daily, and bimonthly indices from NOAA PSL and CPC.

## Usage

``` r
get_index(name, refresh = FALSE, timeout = 30, quiet = FALSE)
```

## Arguments

- name:

  Character. Index name (e.g., "NAO", "Nino34", "PDO"). Use
  [`list_indices()`](https://robustecologies.github.io/oscillatoRs/reference/list_indices.md)
  to see all available indices.

- refresh:

  Logical. If TRUE, forces download even if cached data exists.
  Currently not implemented (reserved for future caching).

- timeout:

  Numeric. Request timeout in seconds. Default is 30.

- quiet:

  Logical. If TRUE, suppresses progress messages.

## Value

A tibble containing the index data. Columns vary by resolution:

- Monthly:

  date, year, month, value, index, description, category, source_url

- Daily:

  date, year, month, day, value, index, description, category,
  source_url

- Bimonthly:

  date, year, month, bimonth, value, index, description, category,
  source_url

Returns NULL if download fails.

## Details

This function downloads data directly from authoritative sources:

- Monthly indices: NOAA Physical Sciences Laboratory (PSL)

- Daily indices: NOAA Climate Prediction Center (CPC) FTP

- MEI v2: NOAA PSL ENSO page

Missing values in source data (typically encoded as -99.9, -999, etc.)
are automatically converted to NA.

## References

Trenberth, K. E. (1997). The definition of El Nino. *Bulletin of the
American Meteorological Society*, 78(12), 2771-2778.
[doi:10.1175/1520-0477(1997)078\<2771:TDOENO\>2.0.CO;2](https://doi.org/10.1175/1520-0477%281997%29078%3C2771%3ATDOENO%3E2.0.CO%3B2)

Mantua, N. J., Hare, S. R., Zhang, Y., Wallace, J. M., & Francis, R. C.
(1997). A Pacific interdecadal climate oscillation with impacts on
salmon production. *Bulletin of the American Meteorological Society*,
78(6), 1069-1079.
[doi:10.1175/1520-0477(1997)078\<1069:APICOW\>2.0.CO;2](https://doi.org/10.1175/1520-0477%281997%29078%3C1069%3AAPICOW%3E2.0.CO%3B2)

Enfield, D. B., Mestas-Nunez, A. M., & Trimble, P. J. (2001). The
Atlantic multidecadal oscillation and its relation to rainfall and river
flows in the continental U.S. *Geophysical Research Letters*, 28(10),
2077-2080.
[doi:10.1029/2000GL012745](https://doi.org/10.1029/2000GL012745)

Hurrell, J. W. (1995). Decadal trends in the North Atlantic Oscillation:
regional temperatures and precipitation. *Science*, 269(5224), 676-679.
[doi:10.1126/science.269.5224.676](https://doi.org/10.1126/science.269.5224.676)

Thompson, D. W. J., & Wallace, J. M. (2000). Annular modes in the
extratropical circulation. Part I: month-to-month variability. *Journal
of Climate*, 13(5), 1000-1016.
[doi:10.1175/1520-0442(2000)013\<1000:AMITEC\>2.0.CO;2](https://doi.org/10.1175/1520-0442%282000%29013%3C1000%3AAMITEC%3E2.0.CO%3B2)

Barnston, A. G., & Livezey, R. E. (1987). Classification, seasonality
and persistence of low-frequency atmospheric circulation patterns.
*Monthly Weather Review*, 115(6), 1083-1126.
[doi:10.1175/1520-0493(1987)115\<1083:CSAPOL\>2.0.CO;2](https://doi.org/10.1175/1520-0493%281987%29115%3C1083%3ACSAPOL%3E2.0.CO%3B2)

## Examples

``` r
if (FALSE) { # \dontrun{
# Download the North Atlantic Oscillation
nao <- get_index("NAO")

# Download Nino 3.4 SST anomalies
nino34 <- get_index("Nino34")

# Download daily Arctic Oscillation
ao_daily <- get_index("AO_daily")
} # }
```
