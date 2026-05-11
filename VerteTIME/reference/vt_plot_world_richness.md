# World richness heatmap on a coarse Robinson grid

For every 5° grid cell intersecting at least one compilation site, plots
the total number of distinct species observed in that cell. Provides a
single-glance richness portrait of the compilation's geography.

## Usage

``` r
vt_plot_world_richness(x, cell_deg = 5)
```

## Arguments

- x:

  A `vt_compilation`.

- cell_deg:

  Numeric grid cell size in degrees.

## Value

A `ggplot`.

## References

Pebesma (2018).

## See also

[`vt_plot_world_map()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_world_map.md),
[`vt_plot_world_hex()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_world_hex.md)

## Examples

``` r
if (FALSE)  data(vertetime); vt_plot_world_richness(vertetime)  # \dontrun{}
```
