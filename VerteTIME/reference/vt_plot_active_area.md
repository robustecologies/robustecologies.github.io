# Stacked-area count of datasets active per calendar year

A dataset is "active" in year `t` when `year_min <= t <= year_max`. The
stacked area is colour-encoded by `taxonomic_focus`.

## Usage

``` r
vt_plot_active_area(x)
```

## Arguments

- x:

  A `vt_compilation`.

## Value

A `ggplot` object.

## References

Dornelas, M., et al. (2014). *Assemblage time series reveal biodiversity
change but not systematic loss*. Science, 344(6181), 296-299.
[doi:10.1126/science.1248484](https://doi.org/10.1126/science.1248484) .

## See also

[`vt_plot_lasagna()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_lasagna.md),
[`vt_plot_onset_cumulative()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_onset_cumulative.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(vertetime); vt_plot_active_area(vertetime)
} # }
```
