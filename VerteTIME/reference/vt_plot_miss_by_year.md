# Missingness rate per year by realm

For every year, plots the proportion of (dataset, site, species) series
for which the annual abundance is `NA`, conditional on the series's
year_min:year_max window. One line per realm.

## Usage

``` r
vt_plot_miss_by_year(x)
```

## Arguments

- x:

  A `vt_compilation`.

## Value

A `ggplot`.

## References

Magurran et al. (2010).

## See also

[`vt_plot_miss_class_heatmap()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_miss_class_heatmap.md)

## Examples

``` r
if (FALSE)  data(vertetime); vt_plot_miss_by_year(vertetime)  # \dontrun{}
```
