# Compute beta-diversity matrices across sites within each dataset

For every dataset with two or more sites, computes a pairwise
dissimilarity between site abundance vectors aggregated over the
dataset's full temporal window. Single-site datasets return `NA`. The
Whittaker multiplicative beta (`gamma / mean(alpha) - 1`) is reported as
a per-dataset scalar in addition to the pairwise matrices.

## Usage

``` r
vt_beta_diversity(x, method = "bray")
```

## Arguments

- x:

  A `vt_compilation`.

- method:

  Character, one or more of `"bray"` (Bray-Curtis), `"jaccard"`
  (presence-absence), `"hellinger"`, `"whittaker"`. Defaults to
  `"bray"`.

## Value

A `tibble` keyed by `(dataset_id, site_a, site_b, method)` with column
`value` (the pairwise dissimilarity for `bray`, `jaccard`, `hellinger`)
plus rows where `site_a == site_b == NA` and the column
`method == "whittaker"` carries the per-dataset scalar.

## Details

Bray-Curtis is the abundance-based dissimilarity \\1 - 2 \\sum
\\min(x_i, y_i) / \\sum (x_i + y_i)\\. Jaccard is the presence-absence
equivalent \\1 - \|A \\cap B\| / \|A \\cup B\|\\. Hellinger is the chord
distance after square-root transformation, \\\\sqrt{0.5 \\sum
(\\sqrt{p_i} - \\sqrt{q_i})^2}\\. Whittaker multiplicative beta is the
ratio of dataset-level richness (`gamma`) to mean per-site-year richness
(`alpha`), expressed as `gamma / mean(alpha) - 1`; values around zero
indicate spatial-temporal coherence in composition, larger values
indicate strong turnover.

## References

Bray, J. R., & Curtis, J. T. (1957). *An ordination of the upland forest
communities of southern Wisconsin*. Ecological Monographs, 27(4),
325-349. [doi:10.2307/1942268](https://doi.org/10.2307/1942268) .

Whittaker, R. H. (1960). *Vegetation of the Siskiyou Mountains, Oregon
and California*. Ecological Monographs, 30(3), 279-338.
[doi:10.2307/1943563](https://doi.org/10.2307/1943563) .

Legendre, P., & Gallagher, E. D. (2001). *Ecologically meaningful
transformations for ordination of species data*. Oecologia, 129(2),
271-280.
[doi:10.1007/s004420100716](https://doi.org/10.1007/s004420100716) .

## See also

[`vt_temporal_turnover()`](https://robustecologies.github.io/VerteTIME/reference/vt_temporal_turnover.md),
[`vt_beta_partition()`](https://robustecologies.github.io/VerteTIME/reference/vt_beta_partition.md),
[`vt_alpha_diversity()`](https://robustecologies.github.io/VerteTIME/reference/vt_alpha_diversity.md)

## Examples

``` r
if (FALSE) { # \dontrun{
co <- vt_ingest_all()
b  <- vt_beta_diversity(co, method = c("bray","whittaker"))
head(b)
} # }
```
