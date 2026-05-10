# Coverage-based rarefaction per site and year

Computes a rarefied richness estimate for each
`(dataset_id, site_id, year)`: given a target sample size `m` (the
minimum non-NA total count across the compilation by default), returns
the expected number of species in a random sample of size `m` from each
year-site abundance vector.

## Usage

``` r
vt_rarefaction(x, sample = NULL)
```

## Arguments

- x:

  A `vt_compilation` or `vt_dataset`.

- sample:

  Optional integer sample size. When `NULL` (default), uses the minimum
  total count across all year-sites with at least one observation.

## Value

A `tibble` keyed by `(dataset_id, site_id, year)` with columns
`n_total`, `sample`, `S_obs` (observed richness), `S_rarefied` (expected
richness at the sample size).

## Details

Sanders' (1968) rarefaction formula is used: \\E\[S_m\] = \\sum_i
\\left\[ 1 - \\binom{N - x_i}{m} / \\binom{N}{m} \\right\]\\ where
\\x_i\\ is the count of species \\i\\ and \\N\\ the total count.
Year-sites with `n_total < sample` return `NA` for `S_rarefied`.

## References

Sanders, H. L. (1968). *Marine benthic diversity: a comparative study*.
The American Naturalist, 102(925), 243-282.
[doi:10.1086/282541](https://doi.org/10.1086/282541) .

Hurlbert, S. H. (1971). *The nonconcept of species diversity: a critique
and alternative parameters*. Ecology, 52(4), 577-586.
[doi:10.2307/1934145](https://doi.org/10.2307/1934145) .

## See also

[`vt_alpha_diversity()`](https://robustecologies.github.io/VerteTIME/reference/vt_alpha_diversity.md),
[`vt_dominance_curve()`](https://robustecologies.github.io/VerteTIME/reference/vt_dominance_curve.md)

## Examples

``` r
if (FALSE) { # \dontrun{
co <- vt_ingest_all()
r <- vt_rarefaction(co, sample = 50)
} # }
```
