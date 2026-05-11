# Span versus coverage scatter (hexbin density)

Per-dataset span (years) against the same dataset's per-species mean
coverage proportion (observed years over span). High-quality datasets
sit towards the top-right.

## Usage

``` r
vt_plot_span_coverage_hex(x)
```

## Arguments

- x:

  A `vt_compilation`.

## Value

A `ggplot` object.

## References

Pebesma, E., & Bivand, R. (2023). *Spatial Data Science: With
Applications in R*. Chapman and Hall/CRC. ISBN 9781032538044.

## See also

[`vt_plot_span_hist()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_span_hist.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(vertetime); vt_plot_span_coverage_hex(vertetime)
} # }
```
