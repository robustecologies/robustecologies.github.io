# Distribution of an alpha-diversity index across the compilation

Per `(dataset, site, year)` value of the chosen index, plotted as a
density and grouped by `taxonomic_focus`.

## Usage

``` r
vt_plot_alpha_hist(x, index = "H")
```

## Arguments

- x:

  A `vt_compilation`.

- index:

  Character index name passed to
  [`vt_alpha_diversity()`](https://robustecologies.github.io/VerteTIME/reference/vt_alpha_diversity.md).

## Value

A `ggplot` object.

## References

Hill, M. O. (1973). *Diversity and evenness*. Ecology, 54(2), 427-432.
[doi:10.2307/1934352](https://doi.org/10.2307/1934352) .

## See also

[`vt_alpha_diversity()`](https://robustecologies.github.io/VerteTIME/reference/vt_alpha_diversity.md),
[`vt_plot_alpha_timeline()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_alpha_timeline.md)

## Examples

``` r
if (FALSE) { # \dontrun{
vt_plot_alpha_hist(vt_ingest_all(), index = "H")
} # }
```
