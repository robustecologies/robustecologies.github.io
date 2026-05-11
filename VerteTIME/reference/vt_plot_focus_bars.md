# Bar charts of dataset count and species count by taxonomic focus

Two stacked panels: number of datasets per `taxonomic_focus` and number
of distinct species per `taxonomic_focus`. The package recommendation is
to interpret these counts as a coverage diagnostic, not as an ecological
signal: skew towards a focus reflects monitoring infrastructure rather
than biological diversity.

## Usage

``` r
vt_plot_focus_bars(x)
```

## Arguments

- x:

  A `vt_compilation`.

## Value

A `ggplot` object.

## References

Magurran, A. E. (2013). *Measuring Biological Diversity*. John Wiley &
Sons. ISBN 9780632056330.

## See also

[`vt_plot_top_species()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_top_species.md),
[`vt_palette()`](https://robustecologies.github.io/VerteTIME/reference/vt_palette.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(vertetime); co <- vertetime
vt_plot_focus_bars(co)
} # }
```
