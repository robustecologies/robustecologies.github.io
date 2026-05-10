# Built-in table of biodiversity time-series compilations cited by VerteTIME

Returns a curated `tibble` listing the major biodiversity time-series
compilations cited in the manuscript timeline figure. The list mixes
compilations of different taxonomic scope: vertebrate-only products
(LPD, RivFishTIME, BBS, PECBMS, CBC, eBird Trends, TEAM/Wildlife
Insights), multi-taxon products that include vertebrates among other
groups (GPDD, BioTIME, PREDICTS) and VerteTIME itself. The taxonomic
scope of every entry is recorded explicitly in the `taxonomic_scope`
column; the figure is positional only and never asserts that all listed
compilations are vertebrate-restricted. Inclusion here is comparative
only; it is not a statement of provenance.

## Usage

``` r
vt_compilations()
```

## Value

A `tibble` with one row per compilation and columns `compilation`,
`first_release_year`, `latest_release_year`, `community_data` (logical;
does the compilation publish multi-species community structure or is it
strictly population-level?), `geographic_scope`, `taxonomic_scope`,
`release_kind`, `notes`.

## Details

The list is intentionally non-exhaustive; users can extend it via the
`extra` argument of
[`vt_plot_database_timeline()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_database_timeline.md).

## References

Hudson, L. N., et al. (2017). *PREDICTS database*. Ecology and
Evolution, 7(1), 145-188.
[doi:10.1002/ece3.2579](https://doi.org/10.1002/ece3.2579) .

Dornelas, M., et al. (2018). *BioTIME*. Global Ecology and Biogeography,
27(7), 760-786.
[doi:10.1111/geb.12729](https://doi.org/10.1111/geb.12729) .

## See also

[`vt_plot_database_timeline()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_database_timeline.md)

## Examples

``` r
vt_compilations()
#> # A tibble: 12 × 8
#>    compilation             first_release_year latest_release_year community_data
#>    <chr>                                <int>               <int> <lgl>         
#>  1 GPDD                                  1994                2010 FALSE         
#>  2 Christmas Bird Count                  1900                2025 TRUE          
#>  3 North American BBS                    1966                2024 TRUE          
#>  4 PECBMS                                1980                2024 TRUE          
#>  5 TEAM / Wildlife Insigh…               2010                2024 TRUE          
#>  6 LPD                                   2013                2024 FALSE         
#>  7 BioTIME v1                            2018                2018 TRUE          
#>  8 BioTIME v2                            2025                2025 TRUE          
#>  9 RivFishTIME                           2021                2021 TRUE          
#> 10 eBird Trends                          2018                2025 TRUE          
#> 11 PREDICTS                              2014                2024 TRUE          
#> 12 VerteTIME                             2026                2026 TRUE          
#> # ℹ 4 more variables: geographic_scope <chr>, taxonomic_scope <chr>,
#> #   release_kind <chr>, notes <chr>
```
