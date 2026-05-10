# Top-N countries by number of datasets

Top-N countries by number of datasets

## Usage

``` r
vt_plot_country_bar(x, n = 15)
```

## Arguments

- x:

  A `vt_compilation`.

- n:

  Integer. Default 15.

## Value

A `ggplot` horizontal bar.

## References

Wilkinson et al. (2016) FAIR principles.

## See also

[`vt_plot_world_map()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_world_map.md)

## Examples

``` r
if (FALSE)  vt_plot_country_bar(vt_ingest_all())  # \dontrun{}
```
