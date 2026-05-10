# Top-N most monitored species across the compilation

Top-N most monitored species across the compilation

## Usage

``` r
vt_plot_top_species(x, n = 20)
```

## Arguments

- x:

  A `vt_compilation`.

- n:

  Integer; default 20. Number of species to list.

## Value

A `ggplot` object: horizontal bar chart with the species ordered from
most to least observed.

## References

Inger, R., et al. (2015). *Common European birds are declining rapidly
while less abundant species' numbers are rising*. Ecology Letters,
18(1), 28-36. [doi:10.1111/ele.12387](https://doi.org/10.1111/ele.12387)
.

## See also

[`vt_plot_focus_bars()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_focus_bars.md)

## Examples

``` r
if (FALSE) { # \dontrun{
vt_plot_top_species(vt_ingest_all(), n = 25)
} # }
```
