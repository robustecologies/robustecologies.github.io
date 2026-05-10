# Beta-diversity / temporal-turnover heatmap

Renders a `(year_a, year_b)` heatmap of pairwise dissimilarity for every
site, faceted by `dataset_id`.

## Usage

``` r
vt_plot_beta_heatmap(x, value_col = NULL)
```

## Arguments

- x:

  A `vt_turnover` (output of
  [`vt_temporal_turnover()`](https://robustecologies.github.io/VerteTIME/reference/vt_temporal_turnover.md)
  or
  [`vt_beta_partition()`](https://robustecologies.github.io/VerteTIME/reference/vt_beta_partition.md))
  or a `vt_compilation` (which triggers a default
  `vt_temporal_turnover(method = "bray", pair = "all_pairs")`
  computation).

- value_col:

  Character name of the value column; defaults to `"value"` for
  [`vt_temporal_turnover()`](https://robustecologies.github.io/VerteTIME/reference/vt_temporal_turnover.md)
  output and `"beta_sor"` for
  [`vt_beta_partition()`](https://robustecologies.github.io/VerteTIME/reference/vt_beta_partition.md)
  output.

## Value

A `ggplot` object.

## References

Anderson, M. J., et al. (2011). *Navigating the multiple meanings of
beta diversity: a roadmap for the practicing ecologist*. Ecology
Letters, 14(1), 19-28.
[doi:10.1111/j.1461-0248.2010.01552.x](https://doi.org/10.1111/j.1461-0248.2010.01552.x)
.

## See also

[`vt_temporal_turnover()`](https://robustecologies.github.io/VerteTIME/reference/vt_temporal_turnover.md),
[`vt_beta_partition()`](https://robustecologies.github.io/VerteTIME/reference/vt_beta_partition.md)

## Examples

``` r
if (FALSE) { # \dontrun{
co <- vt_ingest_all()
tt <- vt_temporal_turnover(co, pair = "all_pairs")
vt_plot_beta_heatmap(tt)
} # }
```
