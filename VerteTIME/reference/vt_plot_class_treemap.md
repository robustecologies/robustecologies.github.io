# Hierarchical Class -\> Order -\> Family treemap

Treemap with three levels of nested rectangles, area proportional to the
number of distinct species. Requires `treemapify` (Suggests) and the
`class` and `family` columns populated by
[`enrich_taxonomy()`](https://robustecologies.github.io/VerteTIME/reference/enrich_taxonomy.md);
the intermediate `order` level is used where present.

## Usage

``` r
vt_plot_class_treemap(x)
```

## Arguments

- x:

  A `vt_compilation`.

## Value

A `ggplot`.

## Details

Species are retained when `class` and `family` are non-missing. The GBIF
backbone promotes the historical reptile orders `Squamata`, `Testudines`
and `Crocodylia` to class rank and returns an empty `order` for them; a
non-missing-order filter would therefore drop every reptile species
silently. The missing `order` is coalesced to the `class` so the class
-\> order -\> family nesting still renders and no class is omitted.

## References

Wilkins (2017). Treemapify.

## See also

[`vt_plot_class_bars()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_class_bars.md)

## Examples

``` r
if (FALSE)  data(vertetime); vt_plot_class_treemap(vertetime)  # \dontrun{}
```
