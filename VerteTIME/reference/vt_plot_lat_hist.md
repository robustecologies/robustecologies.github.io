# Latitudinal histogram of compilation sites, faceted by realm

Histogram of compilation site latitudes binned by `bin_width` degrees,
faceted vertically by realm. Useful for diagnosing latitudinal sampling
gaps and confirming the tropical-temperate balance of the compilation.

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

Magurran, A. E. (2013). *Measuring Biological Diversity*. John Wiley &
Sons. ISBN 9780632056330.

## See also

[`vt_plot_world_map()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_world_map.md),
[`vt_plot_world_hex()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_world_hex.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(vertetime)
vt_plot_lat_hist(vertetime)
} # }
```
