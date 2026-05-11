# Per-site-year community summary

Returns one row per `(dataset_id, site_id, year)` with the most-used
descriptive metrics: richness `S`, Shannon `H`, Pielou `J`, total
abundance, dominance proportion (relative abundance of the most abundant
species).

## Usage

``` r
vt_community_summary(x, by = c("site", "year"))
```

## Arguments

- x:

  A `vt_compilation` or `vt_dataset`.

- by:

  Character, currently only `c("site","year")` is supported.

## Value

A `tibble` keyed by `(dataset_id, site_id, year)` with columns `S`, `H`,
`J`, `total_abundance`, `dominance`.

## References

Magurran, A. E. (2013). *Measuring Biological Diversity*. John Wiley &
Sons. ISBN 9780632056330.

## See also

[`vt_alpha_diversity()`](https://robustecologies.github.io/VerteTIME/reference/vt_alpha_diversity.md),
[`vt_evenness()`](https://robustecologies.github.io/VerteTIME/reference/vt_evenness.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(vertetime); co <- vertetime
s <- vt_community_summary(co)
} # }
```
