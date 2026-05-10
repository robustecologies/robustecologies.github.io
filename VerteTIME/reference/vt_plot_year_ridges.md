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

Joly, F.-X., et al. (2024). *Cohort effects in biodiversity time
series*. Methods in Ecology and Evolution, 15(2), 312-326.
[doi:10.1111/2041-210X.14245](https://doi.org/10.1111/2041-210X.14245) .

## See also

[`vt_plot_lasagna()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_lasagna.md),
[`vt_plot_active_area()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_active_area.md)

## Examples

``` r
if (FALSE) { # \dontrun{
vt_plot_year_ridges(vt_ingest_all(), which = "first")
} # }
```
