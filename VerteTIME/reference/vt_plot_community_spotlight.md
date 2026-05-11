# Featured-community panel: time series + world location for one dataset

Side-by-side composite plot for a single dataset: left panel shows the
community time series via
[`vt_plot_community_series()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_community_series.md);
right panel shows the dataset's site(s) on a Robinson-projected world
map with a halo marker. Useful as a spotlight figure for the manuscript
and the package homepage.

## Usage

``` r
vt_plot_community_spotlight(x, dataset_id)
```

## Arguments

- x:

  A `vt_compilation`.

- dataset_id:

  Character scalar.

## Value

A `ggplot` (or `patchwork`) composite. When `patchwork` is missing,
returns the time-series panel only.

## References

Wickham, H. (2016). *ggplot2: Elegant Graphics for Data Analysis*.
Springer-Verlag.
[doi:10.1007/978-3-319-24277-4](https://doi.org/10.1007/978-3-319-24277-4)
.

## See also

[`vt_plot_community_series()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_community_series.md),
[`vt_plot_world_map()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_world_map.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(vertetime); co <- vertetime
vt_plot_community_spotlight(co, "VT_001")
} # }
```
