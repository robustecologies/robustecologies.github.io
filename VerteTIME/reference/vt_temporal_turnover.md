# Compute temporal turnover (year-pair beta) per site

For every site with at least two observed years, computes the chosen
pairwise dissimilarity between abundance vectors at distinct years.

## Usage

``` r
vt_temporal_turnover(
  x,
  method = c("bray", "jaccard", "hellinger"),
  pair = c("consecutive", "first_vs_each", "all_pairs")
)
```

## Arguments

- x:

  A `vt_compilation` or `vt_dataset`.

- method:

  Character; one of `"bray"`, `"jaccard"`, `"hellinger"`.

- pair:

  Character; one of `"consecutive"` (default; pairs `(t, t+1)`),
  `"first_vs_each"` (pairs `(year_min, t)` for every other observed
  year), `"all_pairs"` (every distinct pair `t_a < t_b`).

## Value

A
[`vt_turnover()`](https://robustecologies.github.io/VerteTIME/reference/vt_turnover.md)
tibble keyed by `(dataset_id, site_id, year_a, year_b)` with columns
`metric` and `value`.

## References

Magurran, A. E., et al. (2010). *Long-term datasets in biodiversity
research and monitoring: assessing change in ecological communities
through time*. Trends in Ecology & Evolution, 25(10), 574-582.
[doi:10.1016/j.tree.2010.06.016](https://doi.org/10.1016/j.tree.2010.06.016)
.

## See also

[`vt_beta_diversity()`](https://robustecologies.github.io/VerteTIME/reference/vt_beta_diversity.md),
[`vt_beta_partition()`](https://robustecologies.github.io/VerteTIME/reference/vt_beta_partition.md),
[`vt_turnover()`](https://robustecologies.github.io/VerteTIME/reference/vt_turnover.md)

## Examples

``` r
if (FALSE) { # \dontrun{
co <- vt_ingest_all()
tt <- vt_temporal_turnover(co, method = "bray", pair = "consecutive")
summary(tt$value)
} # }
```
