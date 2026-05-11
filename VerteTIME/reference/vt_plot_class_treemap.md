# Hierarchical Class -\> Order -\> Family treemap

Treemap with three levels of nested rectangles, area proportional to the
number of distinct species. Requires `treemapify` (Suggests) and the
`class`/`order`/`family` columns populated by
[`enrich_taxonomy()`](https://robustecologies.github.io/VerteTIME/reference/enrich_taxonomy.md).

## Usage

``` r
vt_plot_class_treemap(x)
```

## Arguments

- x:

  A `vt_compilation`.

## Value

A `ggplot`.

## References

Wilkins (2017). Treemapify.

## See also

[`vt_plot_class_bars()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_class_bars.md)

## Examples

``` r
if (FALSE)  data(vertetime); vt_plot_class_treemap(vertetime)  # \dontrun{}
```
