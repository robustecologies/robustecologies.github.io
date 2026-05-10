# Multi-metric ridge plot of community diagnostics

Density ridge of richness `S`, Shannon `H` and Pielou `J` across every
(dataset, site, year) triple, faceted by metric. The portrait surfaces
the compilation's spread on three orthogonal axes of community structure
on a single canvas.

## Usage

``` r
vt_plot_metric_ridges(x)
```

## Arguments

- x:

  A `vt_compilation`.

## Value

A `ggplot`.

## References

Hill (1973); Pielou (1966).

## See also

[`vt_plot_alpha_hist()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_alpha_hist.md),
[`vt_plot_hill_profile()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_hill_profile.md)

## Examples

``` r
if (FALSE)  vt_plot_metric_ridges(vt_ingest_all())  # \dontrun{}
```
