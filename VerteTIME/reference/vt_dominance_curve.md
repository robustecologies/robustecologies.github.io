# Dominance curves per site and year

Returns the cumulative-rank-abundance dominance curve per
`(dataset_id, site_id, year)` combination: species ordered from most to
least abundant, with rank, abundance, relative abundance, and cumulative
relative abundance.

## Usage

``` r
vt_dominance_curve(x)
```

## Arguments

- x:

  A `vt_compilation` or `vt_dataset`.

## Value

A `tibble` keyed by `(dataset_id, site_id, year, rank)` with columns
`species_id`, `value`, `relative`, `cumulative`.

## References

Whittaker, R. H. (1965). *Dominance and diversity in land plant
communities*. Science, 147(3655), 250-260.
[doi:10.1126/science.147.3655.250](https://doi.org/10.1126/science.147.3655.250)
.

## See also

[`vt_rank_abundance()`](https://robustecologies.github.io/VerteTIME/reference/vt_rank_abundance.md),
[`vt_plot_whittaker()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_whittaker.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(vertetime); co <- vertetime
dom <- vt_dominance_curve(co)
} # }
```
