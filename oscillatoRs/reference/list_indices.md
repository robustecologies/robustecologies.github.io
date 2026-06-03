# List available climate indices

Returns a tibble of all available climate oscillation indices with their
descriptions, categories, and temporal resolutions.

## Usage

``` r
list_indices(category = NULL, resolution = NULL)
```

## Arguments

- category:

  Character. Filter by category (e.g., "ENSO", "Teleconnection"). If
  NULL (default), returns all indices.

- resolution:

  Character. Filter by temporal resolution: "monthly", "daily", or
  "bimonthly". If NULL (default), returns all resolutions.

## Value

A tibble with columns: index, description, category, resolution.

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
# List all indices
list_indices()
#> # A tibble: 41 × 4
#>    index     description                                     category resolution
#>    <chr>     <chr>                                           <chr>    <chr>     
#>  1 NAO       North Atlantic Oscillation (CPC)                Telecon… monthly   
#>  2 NAO_Jones North Atlantic Oscillation (Jones/CRU)          Telecon… monthly   
#>  3 PNA       Pacific North American Index                    Telecon… monthly   
#>  4 WP        Western Pacific Index                           Telecon… monthly   
#>  5 EPNP      East Pacific/North Pacific Oscillation          Telecon… monthly   
#>  6 EAWR      Eastern Atlantic/Western Russia                 Telecon… monthly   
#>  7 NP        North Pacific Pattern                           Telecon… monthly   
#>  8 NOI       Northern Oscillation Index                      Telecon… monthly   
#>  9 PDO       Pacific Decadal Oscillation                     Pacific… monthly   
#> 10 IPO_TPI   Tripole Index for Interdecadal Pacific Oscilla… Pacific… monthly   
#> # ℹ 31 more rows

# List only ENSO indices
list_indices(category = "ENSO")
#> # A tibble: 9 × 4
#>   index  description                            category resolution
#>   <chr>  <chr>                                  <chr>    <chr>     
#> 1 SOI    Southern Oscillation Index             ENSO     monthly   
#> 2 Nino12 Nino 1+2 SST Anomaly (0-10S, 90W-80W)  ENSO     monthly   
#> 3 Nino3  Nino 3 SST Anomaly (5N-5S, 150W-90W)   ENSO     monthly   
#> 4 Nino34 Nino 3.4 SST Anomaly (5N-5S, 170-120W) ENSO     monthly   
#> 5 Nino4  Nino 4 SST Anomaly (5N-5S, 160E-150W)  ENSO     monthly   
#> 6 ONI    Oceanic Nino Index                     ENSO     monthly   
#> 7 BEST   Bivariate ENSO Timeseries              ENSO     monthly   
#> 8 TNI    Trans-Nino Index                       ENSO     monthly   
#> 9 MEI_v2 Multivariate ENSO Index v2 (bimonthly) ENSO     bimonthly 

# List daily indices
list_indices(resolution = "daily")
#> # A tibble: 4 × 4
#>   index     description                                      category resolution
#>   <chr>     <chr>                                            <chr>    <chr>     
#> 1 NAO_daily North Atlantic Oscillation (daily, 1950-present) Telecon… daily     
#> 2 AO_daily  Arctic Oscillation (daily, 1950-present)         Polar    daily     
#> 3 PNA_daily Pacific North American (daily, 1950-present)     Telecon… daily     
#> 4 AAO_daily Antarctic Oscillation (daily, 1979-present)      Polar    daily     
```
