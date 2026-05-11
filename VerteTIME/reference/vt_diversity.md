# Construct a `vt_diversity` object

Output container of
[`vt_alpha_diversity()`](https://robustecologies.github.io/VerteTIME/reference/vt_alpha_diversity.md)
and
[`vt_evenness()`](https://robustecologies.github.io/VerteTIME/reference/vt_evenness.md);
carries the computed metrics as a long-format tibble keyed by
`(dataset_id, site_id, year)` and the call attributes used by
`print`/`summary`/`plot`.

## Usage

``` r
vt_diversity(values, call_args = list())

# S3 method for class 'vt_diversity'
print(x, ...)

# S3 method for class 'vt_diversity'
summary(object, ...)

# S3 method for class 'vt_diversity'
plot(x, ...)
```

## Arguments

- values:

  A `tibble` keyed by `(dataset_id, site_id, year)` with columns `index`
  and `value`.

- call_args:

  Named list with the function arguments used to build the object
  (echoed by the print method).

- x:

  A `vt_diversity`.

- ...:

  Forwarded to
  [`vt_plot_alpha_timeline()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_alpha_timeline.md).

- object:

  A `vt_diversity`.

## Value

Object of class `vt_diversity` (subclass of `tbl_df`).

A `ggplot` object (invisibly returned, drawn as side effect).

## References

Hill, M. O. (1973). *Diversity and evenness: a unifying notation and its
consequences*. Ecology, 54(2), 427-432.
[doi:10.2307/1934352](https://doi.org/10.2307/1934352) .

## See also

[`vt_alpha_diversity()`](https://robustecologies.github.io/VerteTIME/reference/vt_alpha_diversity.md),
[`vt_evenness()`](https://robustecologies.github.io/VerteTIME/reference/vt_evenness.md),
`print.vt_diversity()`, `summary.vt_diversity()`, `plot.vt_diversity()`,
[`vt_beta_diversity()`](https://robustecologies.github.io/VerteTIME/reference/vt_beta_diversity.md),
[`vt_turnover()`](https://robustecologies.github.io/VerteTIME/reference/vt_turnover.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(vertetime)
a <- vt_alpha_diversity(vertetime, indices = c("S", "H", "q1"))
print(a); summary(a); plot(a)
} # }
```
