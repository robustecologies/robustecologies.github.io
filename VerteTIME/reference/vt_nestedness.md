# Nestedness based on overlap and decreasing fill (NODF)

Computes NODF (Almeida-Neto et al. 2008) on a year-by-species
presence-absence matrix per `(dataset_id, site_id)`. NODF returns a
value in `[0, 100]` where higher values indicate stronger nestedness
(the species of poor years are subsets of the species of richer years).

## Usage

``` r
vt_nestedness(x, metric = "NODF")
```

## Arguments

- x:

  A `vt_compilation` or `vt_dataset`.

- metric:

  Character, currently only `"NODF"` is supported.

## Value

A `tibble` keyed by `(dataset_id, site_id)` with columns `metric`,
`value`, `n_pairs_rows`, `n_pairs_cols`. Sites with fewer than two
observed years return `NA`.

## Details

For each pair of rows ordered by decreasing marginal total, NODF
computes the proportion of species in the smaller-total row that also
occur in the larger-total row, conditional on the marginal totals being
strictly decreasing. The same calculation is performed on columns. The
final score averages row and column contributions, scaled to `[0, 100]`.

## References

Almeida-Neto, M., Guimaraes, P., Guimaraes, P. R., Loyola, R. D., &
Ulrich, W. (2008). *A consistent metric for nestedness analysis in
ecological systems: reconciling concept and measurement*. Oikos, 117(8),
1227-1239.
[doi:10.1111/j.0030-1299.2008.16644.x](https://doi.org/10.1111/j.0030-1299.2008.16644.x)
.

## See also

[`vt_beta_partition()`](https://robustecologies.github.io/VerteTIME/reference/vt_beta_partition.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(vertetime); co <- vertetime
n <- vt_nestedness(co)
} # }
```
