# Temporal sankey of taxonomic focus across decades

For every dataset that crosses a decade boundary, draws an alluvial flow
showing the persistence of taxonomic focus through decades. The figure
is a "biodiversity tide" portrait: where the compilation's coverage has
shifted, where it has held steady.

## Usage

``` r
vt_plot_temporal_sankey(x)
```

## Arguments

- x:

  A `vt_compilation`.

## Value

A `ggplot`. Requires `ggalluvial`.

## References

Brunson (2024).

## See also

[`vt_plot_active_area()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_active_area.md),
[`vt_plot_alluvial()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_alluvial.md)

## Examples

``` r
if (FALSE)  data(vertetime); vt_plot_temporal_sankey(vertetime)  # \dontrun{}
```
