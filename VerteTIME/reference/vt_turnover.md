# Construct a `vt_turnover` object

Output container of
[`vt_temporal_turnover()`](https://robustecologies.github.io/VerteTIME/reference/vt_temporal_turnover.md)
and
[`vt_beta_partition()`](https://robustecologies.github.io/VerteTIME/reference/vt_beta_partition.md);
a long-format tibble keyed by `(dataset_id, site_id, year_a, year_b)`
with columns for the chosen dissimilarity, the turnover and nestedness
components when applicable, and a `metric` label column.

## Usage

``` r
vt_turnover(values, call_args = list())

# S3 method for class 'vt_turnover'
print(x, ...)

# S3 method for class 'vt_turnover'
summary(object, ...)

# S3 method for class 'vt_turnover'
plot(x, ...)
```

## Arguments

- values:

  A `tibble` with the structure described above.

- call_args:

  Named list of function arguments used to build the object.

- x:

  A `vt_turnover`.

- ...:

  Forwarded to
  [`vt_plot_beta_heatmap()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_beta_heatmap.md).

- object:

  A `vt_turnover`.

## Value

Object of class `vt_turnover`.

For `summary`, a `tibble` with one row per turnover metric and the
columns `metric`, `n`, `median`, `q25`, `q75`, `mean`, `sd`. The metric
column reflects the columns present in the object: `beta_sor`,
`beta_sim`, `beta_nes` for a Baselga partition, or the value of the
`metric` column when the object comes from
[`vt_temporal_turnover()`](https://robustecologies.github.io/VerteTIME/reference/vt_temporal_turnover.md).

A `ggplot` object (invisibly returned, drawn as side effect).

## References

Baselga, A. (2010). *Partitioning the turnover and nestedness components
of beta diversity*. Global Ecology and Biogeography, 19(1), 134-143.
[doi:10.1111/j.1466-8238.2009.00490.x](https://doi.org/10.1111/j.1466-8238.2009.00490.x)
.

## See also

[`vt_temporal_turnover()`](https://robustecologies.github.io/VerteTIME/reference/vt_temporal_turnover.md),
[`vt_beta_partition()`](https://robustecologies.github.io/VerteTIME/reference/vt_beta_partition.md),
[`vt_beta_diversity()`](https://robustecologies.github.io/VerteTIME/reference/vt_beta_diversity.md),
`print.vt_turnover()`, `summary.vt_turnover()`, `plot.vt_turnover()`

## Examples

``` r
if (FALSE) { # \dontrun{
data(vertetime)
tt <- vt_temporal_turnover(vertetime, pair = "consecutive")
print(tt); summary(tt); plot(tt)
} # }
```
