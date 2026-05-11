# Top-N countries by number of datasets

Horizontal bar chart of the `n` countries (ISO3 codes) with the most
datasets in the compilation, descending. Counts are read from
`x$datasets$country_iso3`.

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

Wilkinson, M. D., et al. (2016). *The FAIR Guiding Principles for
scientific data management and stewardship*. Scientific Data, 3, 160018.
[doi:10.1038/sdata.2016.18](https://doi.org/10.1038/sdata.2016.18) .

## See also

[`vt_plot_world_map()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_world_map.md),
[`vt_plot_realm_bar()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_realm_bar.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(vertetime)
vt_plot_country_bar(vertetime)
} # }
```
