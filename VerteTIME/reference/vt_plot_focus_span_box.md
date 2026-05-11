# Boxplot of dataset span by taxonomic focus

Per-focus boxplot of `n_years` from `x$datasets`, providing a
non-parametric summary of how long each taxonomic focus has been
observed in the compilation. Useful as a sampling-effort diagnostic
across foci.

## Usage

``` r
vt_plot_focus_span_box(x)
```

## Arguments

- x:

  A `vt_compilation`.

## Value

A `ggplot` object.

## References

Dornelas, M., et al. (2018). *BioTIME: A database of biodiversity time
series for the Anthropocene*. Global Ecology and Biogeography, 27(7),
760-786. [doi:10.1111/geb.12729](https://doi.org/10.1111/geb.12729) .

## See also

[`vt_plot_span_hist()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_span_hist.md),
[`vt_plot_focus_bars()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_focus_bars.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(vertetime)
vt_plot_focus_span_box(vertetime)
} # }
```
