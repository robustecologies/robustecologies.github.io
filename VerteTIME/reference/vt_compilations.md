# Built-in table of biodiversity time-series compilations cited by VerteTIME

Returns a curated `tibble` listing the major biodiversity time-series
compilations cited in the manuscript timeline figure. The list mixes
compilations of different taxonomic scope: vertebrate-only products
(LPI, RivFishTIME, BBS, PECBMS, CBC, eBird Trends, TEAM/Wildlife
Insights), multi-taxon products that include vertebrates among other
groups (GPDD, BioTIME, PREDICTS), and VerteTIME itself. The taxonomic
scope of every entry is recorded explicitly in `taxonomic_scope`; the
figure is positional only and never asserts that all listed compilations
are vertebrate-restricted. Inclusion here is comparative only; it is not
a statement of provenance.

## Usage

``` r
vt_compilations()
```

## Value

A `tibble` with one row per compilation and the following columns:
`compilation` (display label), `publication_year` (integer year of first
citable public release of the database), `citation_key` (Author et al.
YYYY rendered italic in the figure tab), `n_series` (numeric scale as
reported in the cited publication), `n_series_unit` (unit of the count:
`populations`, `studies`, `time series`, `species`, `circles`,
`site-species time series` or `images`), `data_year_min` and
`data_year_max` (first and last calendar year of any observation in the
compilation), `community_data` (logical: does the compilation publish
multi-species community structure or is it strictly population-level?),
`geographic_scope`, `taxonomic_scope`, `licence` (human-readable
string), `notes`.

## Details

Each row's numeric and licence cells are anchored to a primary source
(peer-reviewed paper or official institutional release) cited in
`@references`. Where a project has an annual rolling release,
`publication_year` points to the most-cited paper anchor and
`data_year_max` carries the latest data-year. The list is intentionally
non-exhaustive; users can extend it via the `extra` argument of
[`vt_plot_database_timeline()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_database_timeline.md).

## References

Dornelas, M., et al. (2018). *BioTIME: a database of biodiversity time
series for the Anthropocene*. Global Ecology and Biogeography, 27(7),
760-786. [doi:10.1111/geb.12729](https://doi.org/10.1111/geb.12729) .

Dornelas, M., et al. (2025). *BioTIME 2.0: Expanding and Improving a
Database of Biodiversity Time Series*. Global Ecology and Biogeography,
34(5), e70003.
[doi:10.1111/geb.70003](https://doi.org/10.1111/geb.70003) .

Comte, L., et al. (2021). *RivFishTIME: A global database of fish
time-series to study global change ecology in riverine systems*. Global
Ecology and Biogeography, 30(1), 38-50.
[doi:10.1111/geb.13210](https://doi.org/10.1111/geb.13210) .

Hudson, L. N., et al. (2017). *The database of the PREDICTS project*.
Ecology and Evolution, 7(1), 145-188.
[doi:10.1002/ece3.2579](https://doi.org/10.1002/ece3.2579) .

McRae, L., Deinet, S., & Freeman, R. (2017). *The Diversity-Weighted
Living Planet Index: Controlling for Taxonomic Bias in a Global
Biodiversity Indicator*. PLoS ONE, 12(1), e0169156.
[doi:10.1371/journal.pone.0169156](https://doi.org/10.1371/journal.pone.0169156)
.

Almond, R. E. A., Grooten, M., Juffe Bignoli, D., & Petersen, T. (Eds.)
(2024). *Living Planet Report 2024: A System in Peril*. WWF
International, Gland, Switzerland. <https://livingplanet.panda.org/>.

Brlík, V., et al. (2021). *Long-term and large-scale multispecies
dataset tracking population changes of common European breeding birds*.
Scientific Data, 8, 21.
[doi:10.1038/s41597-021-00804-2](https://doi.org/10.1038/s41597-021-00804-2)
.

## See also

[`vt_plot_database_timeline()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_database_timeline.md)

## Examples

``` r
vt_compilations()
#> # A tibble: 13 × 12
#>    compilation              publication_year citation_key n_series n_series_unit
#>    <chr>                               <int> <chr>           <int> <chr>        
#>  1 GPDD                                 1999 Inchausti &…     4975 populations  
#>  2 LPI v1                               2017 McRae et al…    14152 populations  
#>  3 LPI v2                               2024 Almond et a…    34836 populations  
#>  4 BioTIME v1                           2018 Dornelas et…      361 studies      
#>  5 BioTIME v2                           2025 Dornelas et…      708 studies      
#>  6 RivFishTIME                          2021 Comte et al…    11386 time series  
#>  7 PREDICTS                             2017 Hudson et a…      666 studies      
#>  8 BBS                                  2017 Sauer et al…      544 species      
#>  9 PECBMS                               2021 Brlik et al…      170 species      
#> 10 CBC                                  2017 LeBaron 2017     2500 circles      
#> 11 eBird Trends                         2023 Fink et al.…     2974 species      
#> 12 TEAM / Wildlife Insights             2020 Ahumada et …     2000 camera deplo…
#> 13 VerteTIME                            2026 Pico & Alma…     1802 site-species…
#> # ℹ 7 more variables: data_year_min <int>, data_year_max <int>,
#> #   community_data <lgl>, geographic_scope <chr>, taxonomic_scope <chr>,
#> #   licence <chr>, notes <chr>
```
