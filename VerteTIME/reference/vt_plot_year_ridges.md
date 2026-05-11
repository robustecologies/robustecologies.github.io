# Ridge plot of first-observed year by taxonomic focus

One smoothed density per focus showing the distribution of dataset onset
years. Useful as a cohort diagnostic.

## Usage

``` r
vt_plot_year_ridges(x, which = c("first", "last"))
```

## Arguments

- x:

  A `vt_compilation`.

- which:

  Character, one of `"first"` (default) or `"last"`.

## Value

A `ggplot` object.

## References

Dornelas, M., Gotelli, N. J., McGill, B., Shimadzu, H., Moyes, F.,
Sievers, C., & Magurran, A. E. (2014). *Assemblage time series reveal
biodiversity change but not systematic loss*. Science, 344(6181),
296-299.
[doi:10.1126/science.1248484](https://doi.org/10.1126/science.1248484) .

## See also

[`vt_plot_lasagna()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_lasagna.md),
[`vt_plot_active_area()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_active_area.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(vertetime); vt_plot_year_ridges(vertetime, which = "first")
} # }
```
