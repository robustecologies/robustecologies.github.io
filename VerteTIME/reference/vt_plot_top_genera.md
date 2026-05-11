# Top-N genera by number of datasets they appear in

Horizontal bar chart of the `n` genera (the leading token before the
underscore in `species_id`) appearing in the largest number of distinct
datasets, coloured by GBIF class. Complements
[`vt_plot_top_species()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_top_species.md)
at coarser taxonomic resolution.

## Usage

``` r
vt_plot_top_genera(x, n = 30)
```

## Arguments

- x:

  A `vt_compilation`.

- n:

  Integer. Default 30.

## Value

A `ggplot`.

## References

Magurran, A. E. (2013). *Measuring Biological Diversity*. John Wiley &
Sons. ISBN 9780632056330.

## See also

[`vt_plot_top_species()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_top_species.md),
[`vt_plot_class_bars()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_class_bars.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(vertetime)
vt_plot_top_genera(vertetime)
} # }
```
