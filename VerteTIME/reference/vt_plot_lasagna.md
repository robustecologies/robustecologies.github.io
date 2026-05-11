# Lasagna plot of dataset-by-year presence

Heat-style raster with one row per dataset and one cell per calendar
year. A filled cell indicates that the dataset had at least one observed
(non-`NA`) abundance value in that year. Datasets are ordered by
`year_min` so the cohort structure is visible.

## Usage

``` r
vt_plot_lasagna(x, sample_n = 200L)
```

## Arguments

- x:

  A `vt_compilation`.

- sample_n:

  Integer; if the compilation has more datasets than `sample_n`, a
  random sample of that size is drawn and plotted. `NULL` plots all.

## Value

A `ggplot` object.

## References

Bridge, E. S., et al. (2011). *Visualizing time-series data with the
lasagna plot*. Open Bioinformatics Journal, 5, 89-96.

## See also

[`vt_plot_year_ridges()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_year_ridges.md),
[`vt_plot_active_area()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_active_area.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(vertetime); vt_plot_lasagna(vertetime)
} # }
```
