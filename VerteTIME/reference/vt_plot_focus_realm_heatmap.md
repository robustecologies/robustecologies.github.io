# Count heatmap of taxonomic focus by realm

Tile heatmap with `taxonomic_focus` on the vertical axis, `realm` on the
horizontal axis, and cell fill proportional to the number of datasets in
that focus-realm combination. The cell text repeats the count for
immediate scan-ability.

## Usage

``` r
vt_plot_focus_realm_heatmap(x)
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

[`vt_plot_realm_bar()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_realm_bar.md),
[`vt_plot_focus_bars()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_focus_bars.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(vertetime)
vt_plot_focus_realm_heatmap(vertetime)
} # }
```
