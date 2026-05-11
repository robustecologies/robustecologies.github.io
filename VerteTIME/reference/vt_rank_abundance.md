# Rank-abundance per site and year

Convenience wrapper of
[`vt_dominance_curve()`](https://robustecologies.github.io/VerteTIME/reference/vt_dominance_curve.md)
that returns only `rank`, `species_id` and `value` columns (no relative
or cumulative).

## Usage

``` r
vt_rank_abundance(x)
```

## Arguments

- x:

  A `vt_compilation` or `vt_dataset`.

## Value

A `tibble` keyed by `(dataset_id, site_id, year, rank)` with columns
`species_id` and `value`.

## References

McGill, B. J. (2011). *Species abundance distributions*. In Magurran &
McGill (Eds.), *Biological Diversity*, Oxford University Press. ISBN
9780199580675.

## See also

[`vt_dominance_curve()`](https://robustecologies.github.io/VerteTIME/reference/vt_dominance_curve.md),
[`vt_plot_rac()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_rac.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(vertetime); co <- vertetime
rac <- vt_rank_abundance(co)
} # }
```
