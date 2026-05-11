# Year-by-class missingness heatmap

Tile heatmap with calendar year on the horizontal axis, GBIF class on
the vertical axis, and cell fill proportional to the missingness rate of
observations in that (class, year) cell. Reveals temporal gaps in
class-specific coverage that are hidden when missingness is averaged
across all classes.

## Usage

``` r
vt_plot_miss_class_heatmap(x)
```

## Arguments

- x:

  A `vt_compilation`.

## Value

A `ggplot`.

## References

Magurran, A. E., et al. (2010). *Long-term datasets in biodiversity
research and monitoring*. Trends in Ecology & Evolution, 25(10),
574-582.
[doi:10.1016/j.tree.2010.06.016](https://doi.org/10.1016/j.tree.2010.06.016)
.

## See also

[`vt_plot_miss_by_year()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_miss_by_year.md),
[`vt_plot_class_bars()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_class_bars.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(vertetime)
vt_plot_miss_class_heatmap(vertetime)
} # }
```
