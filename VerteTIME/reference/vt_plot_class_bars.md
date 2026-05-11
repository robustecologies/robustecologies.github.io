# Two-panel bar chart of series and species per GBIF class

Left panel: number of observation series (one per
`(dataset_id, site_id, species_id)`) per class. Right panel: number of
distinct species per class. Mirrors the LPD `tax-class-bars` figure
adapted to a community compilation where each dataset contributes
multiple species, hence multiple series.

## Usage

``` r
vt_plot_class_bars(x)
```

## Arguments

- x:

  A `vt_compilation` whose species table has `class` populated by
  [`enrich_taxonomy()`](https://robustecologies.github.io/VerteTIME/reference/enrich_taxonomy.md).

## Value

A `ggplot` (two-panel via `facet_wrap`).

## References

Magurran (2013).

## See also

[`vt_plot_class_treemap()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_class_treemap.md),
[`vt_plot_top_species()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_top_species.md)

## Examples

``` r
if (FALSE)  data(vertetime); vt_plot_class_bars(vertetime)  # \dontrun{}
```
