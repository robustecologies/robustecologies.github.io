# Get detailed information about a climate index

Returns comprehensive metadata for a single climate index, including its
description, physical interpretation, data source, and temporal
coverage.

## Usage

``` r
get_index_info(name)
```

## Arguments

- name:

  Character. Index name (e.g., "NAO", "Nino34", "PDO").

## Value

A list containing:

- name:

  Index identifier

- description:

  Short description

- long_description:

  Detailed physical interpretation

- category:

  Index category

- resolution:

  Temporal resolution

- start_year:

  Approximate start year of data

- url:

  Data source URL

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
get_index_info("NAO")
#> $name
#> [1] "NAO"
#> 
#> $description
#> [1] "North Atlantic Oscillation (CPC)"
#> 
#> $long_description
#> [1] "Sea level pressure difference between the Icelandic Low and the Azores High. Positive NAO indicates stronger westerlies and mild, wet winters over northern Europe."
#> 
#> $category
#> [1] "Teleconnection"
#> 
#> $resolution
#> [1] "monthly"
#> 
#> $start_year
#> [1] 1950
#> 
#> $url
#> [1] "https://psl.noaa.gov/data/correlation/nao.data"
#> 
get_index_info("Nino34")
#> $name
#> [1] "Nino34"
#> 
#> $description
#> [1] "Nino 3.4 SST Anomaly (5N-5S, 170-120W)"
#> 
#> $long_description
#> [1] "SST anomaly in the central-eastern equatorial Pacific. Most commonly used ENSO monitoring region due to strong atmosphere-ocean coupling."
#> 
#> $category
#> [1] "ENSO"
#> 
#> $resolution
#> [1] "monthly"
#> 
#> $start_year
#> [1] 1950
#> 
#> $url
#> [1] "https://psl.noaa.gov/data/correlation/nina34.anom.data"
#> 
```
