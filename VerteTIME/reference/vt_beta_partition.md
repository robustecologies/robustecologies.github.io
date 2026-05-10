# Partition Sorensen dissimilarity into turnover and nestedness components

Implements the Baselga (2010) decomposition of Sorensen pairwise
dissimilarity into a turnover (Simpson) component and a
nestedness-resultant component. Operates on presence-absence vectors per
year within each site.

## Usage

``` r
vt_beta_partition(x, family = "baselga")
```

## Arguments

- x:

  A `vt_compilation` or `vt_dataset`.

- family:

  Character, currently only `"baselga"` is supported.

## Value

A
[`vt_turnover()`](https://robustecologies.github.io/VerteTIME/reference/vt_turnover.md)
tibble keyed by `(dataset_id, site_id, year_a, year_b)` with columns
`beta_sor`, `beta_sim`, `beta_nes`. `beta_sim` is the Simpson turnover
component; `beta_nes = beta_sor - beta_sim` is the nestedness-resultant
component.

## Details

For two presence-absence vectors with `a` shared species, `b` species
unique to the first and `c` species unique to the second:

- \\\\beta\_{SOR} = (b + c) / (2a + b + c)\\ — total Sorensen
  dissimilarity;

- \\\\beta\_{SIM} = \\min(b, c) / (a + \\min(b, c))\\ — Simpson
  turnover;

- \\\\beta\_{NES} = \\beta\_{SOR} - \\beta\_{SIM}\\ —
  nestedness-resultant.

Pairs with `a + b + c == 0` (both years empty) return `NA` triplets.

## References

Baselga, A. (2010). *Partitioning the turnover and nestedness components
of beta diversity*. Global Ecology and Biogeography, 19(1), 134-143.
[doi:10.1111/j.1466-8238.2009.00490.x](https://doi.org/10.1111/j.1466-8238.2009.00490.x)
.

## See also

[`vt_temporal_turnover()`](https://robustecologies.github.io/VerteTIME/reference/vt_temporal_turnover.md),
[`vt_beta_diversity()`](https://robustecologies.github.io/VerteTIME/reference/vt_beta_diversity.md)

## Examples

``` r
if (FALSE) { # \dontrun{
co <- vt_ingest_all()
bp <- vt_beta_partition(co)
head(bp)
} # }
```
