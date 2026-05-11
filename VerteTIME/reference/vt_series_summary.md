# Per-series temporal-coverage and variability metrics for a long observation table

Computes the standard set of derived per-series metrics used by the
package: number of observed years, observed window endpoints, span,
coverage, missingness, internal-gap structure, mean abundance and
coefficient of variation. The output is keyed by
`(dataset_id, site_id, species_id)` and joinable back to the compilation
tables.

## Usage

``` r
vt_series_summary(long)
```

## Arguments

- long:

  A `data.frame` or `data.table` in long format with at least columns
  `dataset_id`, `site_id`, `species_id`, `year`, `value`.

## Value

A `tibble` with one row per `(dataset_id, site_id, species_id)`
combination and the columns `n_obs`, `first_year`, `last_year`, `span`,
`coverage_pct`, `longest_internal_gap`, `n_internal_gaps`, `mean_value`,
`cv_value`, `log10_mean`. Coefficient of variation is set to `NA`
whenever the mean is non-positive, since it is not interpretable in that
regime.

## Details

The function distinguishes structural missingness (years outside
`[first_year, last_year]`) from internal missingness (NAs inside the
observed window). Use `coverage_pct` and `longest_internal_gap` rather
than a global missingness percentage for quality cuts.

## References

Dornelas, M., et al. (2018). BioTIME: a database of biodiversity time
series for the Anthropocene. *Global Ecology and Biogeography*, 27(7),
760-786. [doi:10.1111/geb.12729](https://doi.org/10.1111/geb.12729) .

## See also

[`vt_compilation_summary()`](https://robustecologies.github.io/VerteTIME/reference/vt_compilation_summary.md),
[`vt_validate()`](https://robustecologies.github.io/VerteTIME/reference/vt_validate.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(vertetime); co <- vertetime
ss <- vt_series_summary(co$observations)
summary(ss$span)
} # }
```
