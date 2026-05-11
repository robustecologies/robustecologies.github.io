# Datasets per realm

Vertical bar chart of dataset count per realm (Terrestrial, Freshwater,
Marine, plus an `unknown` bucket for uncatalogued entries). Reads
`x$datasets$realm`.

## Usage

``` r
vt_plot_realm_bar(x)
```

## Arguments

- x:

  A `vt_compilation`.

## Value

A `ggplot`.

## References

Dornelas, M., et al. (2018). *BioTIME: A database of biodiversity time
series for the Anthropocene*. Global Ecology and Biogeography, 27(7),
760-786. [doi:10.1111/geb.12729](https://doi.org/10.1111/geb.12729) .

## See also

[`vt_plot_focus_realm_heatmap()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_focus_realm_heatmap.md),
[`vt_plot_country_bar()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_country_bar.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(vertetime)
vt_plot_realm_bar(vertetime)
} # }
```
