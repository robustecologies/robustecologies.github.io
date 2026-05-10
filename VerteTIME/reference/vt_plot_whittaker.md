# Whittaker rank-abundance plot

Log-abundance versus rank ordering, the canonical species-abundance
distribution diagnostic introduced by Whittaker (1965). One panel per
dataset by default; the user can subset prior to calling the function
for a single-dataset view.

## Usage

``` r
vt_plot_whittaker(x, year = NULL, facet_by = c("dataset", "dataset_year"))
```

## Arguments

- x:

  A `vt_compilation` or `vt_dataset`.

- year:

  Optional integer; when supplied, restricts to that calendar year (per
  site).

- facet_by:

  Character, one of `"dataset"` (default) or `"dataset_year"` for
  small-multiples by dataset and year.

## Value

A `ggplot` object.

## Details

Plots are produced with
[`vt_theme()`](https://robustecologies.github.io/VerteTIME/reference/vt_theme.md)
and use the
[`vt_palette()`](https://robustecologies.github.io/VerteTIME/reference/vt_palette.md)
six-level colour scheme. The y axis is log10-transformed; the x axis is
integer rank. A grey caption attributes the data to VerteTIME and notes
that abundance values are not normalised across datasets, so cross-panel
comparison should be qualitative.

## References

Whittaker, R. H. (1965). *Dominance and diversity in land plant
communities*. Science, 147(3655), 250-260.
[doi:10.1126/science.147.3655.250](https://doi.org/10.1126/science.147.3655.250)
.

## See also

[`vt_plot_rac()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_rac.md),
[`vt_dominance_curve()`](https://robustecologies.github.io/VerteTIME/reference/vt_dominance_curve.md)

## Examples

``` r
if (FALSE) { # \dontrun{
co <- vt_ingest_all()
p <- vt_plot_whittaker(vt_filter(co, dataset_id == "VT_001"))
print(p)
} # }
```
