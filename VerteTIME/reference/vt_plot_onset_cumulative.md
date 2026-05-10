# Cumulative onset curve

Cumulative count of datasets with `year_min <= t` plotted against `t`.
The slope is the rate at which new datasets enter the compilation over
time.

## Usage

``` r
vt_plot_onset_cumulative(x)
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

[`vt_plot_active_area()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_active_area.md)

## Examples

``` r
if (FALSE) { # \dontrun{
vt_plot_onset_cumulative(vt_ingest_all())
} # }
```
