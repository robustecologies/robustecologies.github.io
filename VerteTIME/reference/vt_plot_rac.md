# Rank-abundance curves with year facets

Same data as
[`vt_plot_whittaker()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_whittaker.md)
but on a linear y axis and with one panel per `(dataset, year)`
combination.

## Usage

``` r
vt_plot_rac(x, year = NULL)
```

## Arguments

- x:

  A `vt_compilation` or `vt_dataset`.

- year:

  Optional integer; when supplied, restricts to that calendar year (per
  site).

## Value

A `ggplot` object.

## References

McGill, B. J. (2011). *Species abundance distributions*. In Magurran &
McGill (Eds.), *Biological Diversity*, Oxford University Press. ISBN
9780199580675.

## See also

[`vt_plot_whittaker()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_whittaker.md),
[`vt_rank_abundance()`](https://robustecologies.github.io/VerteTIME/reference/vt_rank_abundance.md)

## Examples

``` r
if (FALSE) { # \dontrun{
co <- vt_ingest_all()
vt_plot_rac(vt_filter(co, dataset_id == "VT_001"))
} # }
```
