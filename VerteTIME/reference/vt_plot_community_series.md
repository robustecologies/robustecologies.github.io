# Plot the abundance trajectory of every species in a community

Single-community time-series plot: one line per `species_id` over
calendar year, all species shown together with strongly contrasted
colours so that the user can read composition change at a glance. Unlike
the LPD-style single-population plot, this shows the full community in
one panel and is the canonical way to inspect a `vt_dataset`.

## Usage

``` r
vt_plot_community_series(
  x,
  dataset_id = NULL,
  log_y = TRUE,
  free_scales = TRUE,
  ncol = 3,
  max_species_in_legend = 15
)
```

## Arguments

- x:

  A `vt_dataset` or `vt_compilation`. When a compilation is passed and
  `dataset_id` is `NULL`, the function plots every dataset in the
  compilation as small multiples; pass a single `dataset_id` to focus on
  one.

- dataset_id:

  Optional character scalar to restrict to one dataset.

- log_y:

  Logical. When `TRUE` (default), the y axis is `log10(value + 1)` so
  that low-abundance species remain readable next to dominants.

- free_scales:

  Logical. When `TRUE` (default) and faceting by dataset, y scales are
  free per dataset. When `FALSE`, all datasets share a y scale.

- ncol:

  Integer; number of columns in the small-multiples grid when plotting
  many datasets. Default 3.

## Value

A `ggplot` object.

## Details

Colour assignment uses
[`viridis::turbo`](https://sjmgarnier.github.io/viridisLite/reference/viridis.html)
to maximise contrast across arbitrary species counts. When the dataset
has more species than `max_species_in_legend` (default 15) the legend is
dropped and species identity becomes a visual property without a key.
For multi-site datasets, linetype encodes the site so that overlapping
species traces remain distinguishable.

## References

Magurran, A. E., et al. (2010). *Long-term datasets in biodiversity
research and monitoring: assessing change in ecological communities
through time*. Trends in Ecology & Evolution, 25(10), 574-582.
[doi:10.1016/j.tree.2010.06.016](https://doi.org/10.1016/j.tree.2010.06.016)
.

## See also

[`vt_plot_alpha_timeline()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_alpha_timeline.md),
[`vt_plot_whittaker()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_whittaker.md),
[`vt_plot_world_map()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_world_map.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(vertetime); co <- vertetime
vt_plot_community_series(co, dataset_id = "VT_001")
vt_plot_community_series(vt_filter(co, taxonomic_focus == "birds"))
} # }
```
