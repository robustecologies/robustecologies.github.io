# Monthly climate indices in wide format

The same data as
[climate_monthly](https://robustecologies.github.io/oscillatoRs/reference/climate_monthly.md)
but pivoted to wide format with one column per index. Useful for
computing correlations and multivariate analyses.

## Usage

``` r
climate_monthly_wide
```

## Format

A tibble with the following structure:

- date:

  Date. First day of the month.

- NAO, PDO, ...:

  Numeric. One column per index containing monthly values.

## Source

NOAA Physical Sciences Laboratory: <https://psl.noaa.gov/>

## Details

Missing values (NA) occur when indices have different temporal coverage
or when source data contains gaps. Use
[`na.omit()`](https://rdrr.io/r/stats/na.fail.html) or
[`complete.cases()`](https://rdrr.io/r/stats/complete.cases.html) to
filter to common observation periods when computing correlations.

## See also

[climate_monthly](https://robustecologies.github.io/oscillatoRs/reference/climate_monthly.md)
for long format,
[`plot_correlation()`](https://robustecologies.github.io/oscillatoRs/reference/plot_correlation.md)
for correlation visualization

## Examples

``` r
data(climate_monthly_wide)
head(climate_monthly_wide)
#> # A tibble: 6 × 37
#>   date         AAO AMO_smoothed AMO_unsmoothed    AO Atlantic_Tripole  BEST
#>   <date>     <dbl>        <dbl>          <dbl> <dbl>            <dbl> <dbl>
#> 1 1870-01-01    NA           NA             NA    NA               NA    NA
#> 2 1870-02-01    NA           NA             NA    NA               NA    NA
#> 3 1870-03-01    NA           NA             NA    NA               NA    NA
#> 4 1870-04-01    NA           NA             NA    NA               NA    NA
#> 5 1870-05-01    NA           NA             NA    NA               NA    NA
#> 6 1870-06-01    NA           NA             NA    NA               NA    NA
#> # ℹ 30 more variables: CAR <dbl>, DMI <dbl>, EAWR <dbl>, EPNP <dbl>,
#> #   Global_Temp <dbl>, IPO_TPI <dbl>, Indian_Monsoon <dbl>, NAO <dbl>,
#> #   NAO_Jones <dbl>, NOI <dbl>, NP <dbl>, NTA <dbl>, Nino12 <dbl>, Nino3 <dbl>,
#> #   Nino34 <dbl>, Nino4 <dbl>, ONI <dbl>, PDO <dbl>, PNA <dbl>,
#> #   Pacific_Warmpool <dbl>, QBO <dbl>, SOI <dbl>, Sahel_Rainfall <dbl>,
#> #   Solar_Flux <dbl>, TNA <dbl>, TNI <dbl>, TSA <dbl>,
#> #   Tropical_Pacific_EOF <dbl>, WHWP <dbl>, WP <dbl>

if (FALSE) { # \dontrun{
# Compute correlation between NAO and AO (requires full data)
cor(climate_monthly_wide$NAO, climate_monthly_wide$AO, use = "complete.obs")
} # }
```
