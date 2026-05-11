# Hexbin density of compilation sites on a Robinson world map

Hexagonal binning of compilation site coordinates on a Robinson
projection world map, coloured by site count per cell with the magma
palette. Falls back to a Cartesian rectangle when `sf` or
`rnaturalearth` are not installed.

## Usage

``` r
vt_plot_world_hex(x)
```

## Arguments

- x:

  A `vt_compilation`.

## Value

A `ggplot`.

## References

Pebesma, E. (2018). *Simple Features for R: Standardized Support for
Spatial Vector Data*. The R Journal, 10(1), 439-446.
[doi:10.32614/RJ-2018-009](https://doi.org/10.32614/RJ-2018-009) .

## See also

[`vt_plot_world_map()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_world_map.md),
[`vt_plot_lat_hist()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_lat_hist.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(vertetime)
vt_plot_world_hex(vertetime)
} # }
```
