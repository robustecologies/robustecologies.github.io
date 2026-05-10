# Histograms of dataset span and observation count

Histograms of dataset span and observation count

## Usage

``` r
vt_plot_span_hist(x)
```

## Arguments

- x:

  A `vt_compilation`.

## Value

A `ggplot` with two panels side-by-side.

## References

Magurran, A. E., et al. (2010). *Long-term datasets in biodiversity
research and monitoring*. Trends in Ecology & Evolution, 25(10),
574-582.
[doi:10.1016/j.tree.2010.06.016](https://doi.org/10.1016/j.tree.2010.06.016)
.

## See also

[`vt_plot_year_ridges()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_year_ridges.md),
[`vt_plot_lasagna()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_lasagna.md)

## Examples

``` r
if (FALSE) { # \dontrun{
vt_plot_span_hist(vt_ingest_all())
} # }
```
