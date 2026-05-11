# Positional timeline of biodiversity time-series compilations

Two-panel figure that places VerteTIME alongside the major biodiversity
time-series compilations cited in the field. The top panel is a
publication-event timeline: a horizontal calendar axis with one filled
marker per compilation at the year of its first citable public release,
with a labelled tab carrying the compilation name, the citation key, the
reported series count, and the licence. The bottom panel is a corrected
Gantt of underlying-data spans (first to last calendar year of any
observation in the compilation), coloured by `kind` (community,
population, or this work).

The compilations shown have heterogeneous taxonomic scope, ranging from
vertebrate-only to multi-taxon, and the figure does not claim a uniform
vertebrate scope across the rows; the `taxonomic_scope` column of
[`vt_compilations()`](https://robustecologies.github.io/VerteTIME/reference/vt_compilations.md)
disambiguates each entry. The figure is positional only: the
compilations shown are not the source of any VerteTIME series.

## Usage

``` r
vt_plot_database_timeline(x = NULL, extra = NULL)
```

## Arguments

- x:

  Optional `vt_compilation`. When supplied, VerteTIME's `data_year_min`
  and `data_year_max` in the bottom panel are repositioned to the actual
  `(year_min, year_max)` of the compilation; `n_series` is recomputed as
  the number of unique `(site_id, species_id)` time series in
  `x$observations`. Default `NULL` uses the canonical values shipped
  with
  [`vt_compilations()`](https://robustecologies.github.io/VerteTIME/reference/vt_compilations.md).

- extra:

  Optional `tibble` with the same shape as
  [`vt_compilations()`](https://robustecologies.github.io/VerteTIME/reference/vt_compilations.md);
  appended to the built-in list before plotting.

## Value

A `patchwork` composite of two `ggplot` panels.

## Details

The top panel places labelled tabs via
[ggrepel::geom_label_repel](https://ggrepel.slowkow.com/reference/geom_text_repel.html)
with a force-based non-overlap solver, seeded with deterministic
above/below nudges so VerteTIME's tab always sits below the axis and the
rest alternate. Each tab carries three lines: compilation name, citation
key, and a technical line of the form `N unit · licence`. The bottom
panel uses a standard `geom_segment` Gantt with rounded endpoints,
sorted by `data_year_min`. The two panels are composed via
[`patchwork::wrap_plots`](https://patchwork.data-imaginist.com/reference/wrap_plots.html)
with `heights = c(2.4, 1)` so the milestone panel dominates visually.

Both panels share the colour mapping: orange for "this work"
(VerteTIME), green for community-level compilations, purple for
population-level. VerteTIME's marker and tab are highlighted in orange.

When the optional packages `ggrepel` or `patchwork` are not available,
the function returns the bottom panel only with a warning.

## References

Slowikowski, K. (2024). *ggrepel: Automatically Position Non-Overlapping
Text Labels with ggplot2*. R package version 0.9.5.
<https://CRAN.R-project.org/package=ggrepel>.

Pedersen, T. L. (2024). *patchwork: The Composer of Plots*. R package
version 1.2.0. <https://CRAN.R-project.org/package=patchwork>.

## See also

[`vt_compilations()`](https://robustecologies.github.io/VerteTIME/reference/vt_compilations.md),
[`vt_plot_world_map()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_world_map.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(vertetime)
p <- vt_plot_database_timeline(vertetime)
print(p)
} # }
```
