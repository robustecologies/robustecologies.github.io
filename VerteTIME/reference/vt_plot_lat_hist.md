# Latitudinal histogram of compilation sites, faceted by realm

Latitudinal histogram of compilation sites, faceted by realm

## Usage

``` r
vt_plot_lat_hist(x, bin_width = 10)
```

## Arguments

- x:

  A `vt_compilation`.

- bin_width:

  Numeric latitudinal bin width in degrees.

## Value

A `ggplot`.

## References

Magurran (2013).

## See also

[`vt_plot_world_map()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_world_map.md)

## Examples

``` r
if (FALSE)  vt_plot_lat_hist(vt_ingest_all())  # \dontrun{}
```
