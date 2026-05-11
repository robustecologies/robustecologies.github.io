# Alluvial flow from taxonomic focus through realm to continent

Three-axis alluvial plot tracing the joint distribution of datasets
across taxonomic focus, realm and continent. Each ribbon is a flow of
datasets through the three axes; ribbon width is proportional to dataset
count. Requires the `ggalluvial` package; continent is derived from
`country_iso3`.

## Usage

``` r
vt_plot_alluvial(x)
```

## Arguments

- x:

  A `vt_compilation`.

## Value

A `ggplot`.

## References

Brunson, J. C. (2020). *ggalluvial: Layered Grammar for Alluvial Plots*.
Journal of Open Source Software, 5(49), 2017.
[doi:10.21105/joss.02017](https://doi.org/10.21105/joss.02017) .

## See also

[`vt_plot_focus_realm_heatmap()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_focus_realm_heatmap.md),
[`vt_plot_realm_bar()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_realm_bar.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(vertetime)
vt_plot_alluvial(vertetime)
} # }
```
