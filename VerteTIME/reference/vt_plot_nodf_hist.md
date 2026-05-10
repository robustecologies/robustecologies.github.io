# NODF nestedness across the compilation

Distribution of NODF over every site of every dataset that has at least
two observed years.

## Usage

``` r
vt_plot_nodf_hist(x)
```

## Arguments

- x:

  A `vt_compilation`.

## Value

A `ggplot` object.

## References

Almeida-Neto, M., et al. (2008). *A consistent metric for nestedness
analysis*. Oikos, 117(8), 1227-1239.
[doi:10.1111/j.0030-1299.2008.16644.x](https://doi.org/10.1111/j.0030-1299.2008.16644.x)
.

## See also

[`vt_nestedness()`](https://robustecologies.github.io/VerteTIME/reference/vt_nestedness.md)

## Examples

``` r
if (FALSE) { # \dontrun{
vt_plot_nodf_hist(vt_ingest_all())
} # }
```
