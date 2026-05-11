# Mean abundance versus coefficient of variation per species

One point per `(dataset, species)` pair. The plot reveals the
compilation's unit-class heterogeneity: clusters of low-CV high-mean and
high-CV low-mean coexist because of the mix of counts, densities and
indices.

## Usage

``` r
vt_plot_mean_cv(x)
```

## Arguments

- x:

  A `vt_compilation`.

## Value

A `ggplot` object.

## References

McGill, B. J. (2011). *Species abundance distributions*. In Magurran &
McGill (Eds.), *Biological Diversity*, Oxford University Press. ISBN
9780199580675.

## See also

[`vt_plot_top_species()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_top_species.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(vertetime); vt_plot_mean_cv(vertetime)
} # }
```
