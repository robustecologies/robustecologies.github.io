# World map of compilation sites

Robinson-projected world map with one point per site. Point size is
mapped to the dataset's `n_species` and colour to `taxonomic_focus`. The
basemap is supplied by
`rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")` when
the suggested packages are available; otherwise the function falls back
to a plain Cartesian rectangle.

## Usage

``` r
vt_plot_world_map(x, projection = c("robinson", "wgs84"))
```

## Arguments

- x:

  A `vt_compilation`.

- projection:

  Character; one of `"robinson"` (default) or `"wgs84"`.

## Value

A `ggplot` object.

## References

Pebesma, E. (2018). *Simple Features for R: Standardized Support for
Spatial Vector Data*. The R Journal, 10(1), 439-446.
[doi:10.32614/RJ-2018-009](https://doi.org/10.32614/RJ-2018-009) .

## See also

[`vt_plot_database_timeline()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_database_timeline.md),
[`vt_compilation()`](https://robustecologies.github.io/VerteTIME/reference/vt_compilation.md)

## Examples

``` r
if (FALSE) { # \dontrun{
co <- vt_ingest_all()
vt_plot_world_map(co)
} # }
```
